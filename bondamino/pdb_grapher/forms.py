"""Forms of Bondamino and the invitation email.

* Login and password reset: ``LoginForm``/``AdminLoginForm`` (login throttling via ``throttle.py``),
  ``PasswordResetForm`` (rate limited, also for invited users without a password).
* User administration: ``InvitationUserForm`` and ``send_invitation`` (used by ``admin.py``),
  ``UserSettingsForm`` (account page).
* Files and notes: ``AvailabilityForm``/``PDBUploadForm`` (who can see a file) and ``NoteForm`` (who can read a
  note). Their ``groups`` fields offer only the user's own groups plus the groups the object is already shared
  with.
"""
import logging
import unicodedata
import zoneinfo

from django import forms
from django.conf import settings
from django.contrib.admin.forms import AdminAuthenticationForm
from django.contrib.auth import forms as auth_forms
from django.contrib.auth import get_user_model
from django.contrib.auth.tokens import default_token_generator
from django.core.mail import send_mail
from django.template.loader import render_to_string
from django.urls import reverse
from django.utils.encoding import force_bytes
from django.utils.http import urlsafe_base64_encode
from django.core.files.uploadedfile import SimpleUploadedFile
from django.utils.translation import gettext_lazy as _

from . import rcsb, throttle
from .models import Language, Note, PDBFile, UserSettings, WorkGroup

logger = logging.getLogger(__name__)


def _same_email(a, b):
    """Compares two email addresses case-insensitively after Unicode normalization (as Django does)."""
    return unicodedata.normalize('NFKC', a).casefold() == unicodedata.normalize('NFKC', b).casefold()


class ThrottledLoginMixin:
    """Mixin for authentication forms that rejects logins after too many failed attempts.

    While locked, even the correct password is rejected, so an attacker cannot tell whether a guess was right
    (limits and duration in ``throttle.py``). Must precede the Django form class in the bases so its ``clean``
    wraps Django's authentication.
    """

    def clean(self):
        """Checks the lock, authenticates and updates the failure counters.

        Raises:
            ValidationError: With code ``'throttled'`` while the login is locked, or Django's own error for a
                wrong password or inactive user (which is counted as a failure).
        """
        username = self.cleaned_data.get('username') or ''
        if self.request is None or not username:
            return super().clean()
        if throttle.login_locked(self.request, username):
            raise forms.ValidationError(
                _('Zu viele fehlgeschlagene Anmeldeversuche. Bitte versuche es in %(minutes)d Minuten erneut.'),
                code='throttled', params={'minutes': throttle.LOGIN_LOCK_SECONDS // 60},
            )
        try:
            cleaned = super().clean()
        except forms.ValidationError as error:
            if error.code in ('invalid_login', 'inactive'):
                throttle.login_failed(self.request, username)
            raise
        throttle.login_succeeded(self.request, username)
        return cleaned


class LoginForm(ThrottledLoginMixin, auth_forms.AuthenticationForm):
    """Throttled login form of the site (wired up in ``bondamino/urls.py``)."""
    pass


class AdminLoginForm(ThrottledLoginMixin, AdminAuthenticationForm):
    """Throttled login form of the admin (wired up in ``admin.py``)."""
    pass


class PasswordResetForm(auth_forms.PasswordResetForm):
    """Forgotten password form.

    Unlike Django's form it also works for users without a usable password, i.e. after an invitation that was
    not yet accepted or has expired. Only a few emails per hour are sent per email address and per IP address;
    the response is the same either way so that nothing is revealed.
    """

    def get_users(self, email):
        """Returns the active users with this email address, including users without a usable password.

        Overrides Django's version, which skips users without a usable password.
        """
        user_model = get_user_model()
        email_field = user_model.get_email_field_name()
        users = user_model._default_manager.filter(**{f'{email_field}__iexact': email, 'is_active': True})
        return (u for u in users if _same_email(email, getattr(u, email_field)))

    def save(self, *args, request=None, **kwargs):
        """Sends the reset email unless the rate limit is reached (then only a warning is logged)."""
        if request is not None and not throttle.reset_allowed(request, self.cleaned_data['email']):
            logger.warning('Password reset rate limit reached for %s', throttle.client_ip(request))
            return
        super().save(*args, request=request, **kwargs)


def send_invitation(request, user):
    """Sends a newly created user an email with the link to set their password.

    The link leads to ``views.SetPasswordView`` and is valid for ``PASSWORD_RESET_TIMEOUT``.

    Args:
        request: The current request, used to build the absolute link.
        user: The invited user; must have an email address.

    Raises:
        Exception: Any error of the mail backend (e.g. SMTP errors) is passed on to the caller.
    """
    path = reverse('pdb_grapher:invitation', kwargs={
        'uidb64': urlsafe_base64_encode(force_bytes(user.pk)),
        'token': default_token_generator.make_token(user),
    })
    context = {
        'user': user,
        'link': request.build_absolute_uri(path),
        'valid_hours': settings.PASSWORD_RESET_TIMEOUT // 3600,
    }
    subject = ''.join(render_to_string('pdb_grapher/email/invitation_subject.txt', context).splitlines())
    body = render_to_string('pdb_grapher/email/invitation_email.txt', context)
    send_mail(subject, body, None, [user.email])


class InvitationUserForm(forms.ModelForm):
    """Admin form for creating a user without a password.

    The person sets the password via the link in the invitation email and can only log in afterwards. The email
    address is therefore required.
    """

    class Meta:
        model = get_user_model()
        fields = ['username', 'email', 'first_name', 'last_name']
        field_classes = {'username': auth_forms.UsernameField}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['email'].required = True

    def save(self, commit=True):
        user = super().save(commit=False)
        user.set_unusable_password()
        if commit:
            user.save()
        return user


class UserSettingsForm(forms.ModelForm):
    """Settings form on the account page: language, theme and time zone.

    Offers the languages created in the admin; without any, the language field is left out (the default
    language applies). The time zone list contains all IANA time zones; the empty choice means the server's
    ``TIME_ZONE``.

    Args:
        current_language: Active language code, preselected if the user has not stored a language yet.
    """
    language = forms.ChoiceField(label=_('Sprache'), widget=forms.Select(attrs={'class': 'form-select'}))
    time_zone = forms.ChoiceField(label=_('Zeitzone'), required=False,
                                  widget=forms.Select(attrs={'class': 'form-select'}))

    class Meta:
        model = UserSettings
        fields = ['language', 'theme', 'time_zone']
        widgets = {'theme': forms.Select(attrs={'class': 'form-select'})}

    def __init__(self, *args, current_language, **kwargs):
        super().__init__(*args, **kwargs)
        languages = [(lang.code, lang.name) for lang in Language.objects.all()]
        if languages:
            self.fields['language'].choices = languages
            if not self.instance.language:
                self.initial['language'] = current_language
        else:
            # A required field without choices could never be submitted
            del self.fields['language']
        self.fields['time_zone'].choices = [
            ('', _('Wie der Server (%(zone)s)') % {'zone': settings.TIME_ZONE}),
            *((zone, zone.replace('_', ' ')) for zone in sorted(zoneinfo.available_timezones())),
        ]


class AvailabilityForm(forms.ModelForm):
    """Availability of a PDB file: personal and/or for work groups.

    Selectable are the user's groups and the groups the file is already shared with (so the uploader can remove
    it from a group even after leaving that group). At least one of personal/groups is required.

    Args:
        user: The uploader, whose groups are offered.
    """
    personal = forms.BooleanField(
        label=_('Für mich persönlich'), required=False, initial=True,
        widget=forms.CheckboxInput(attrs={'class': 'form-check-input'}),
    )
    groups = forms.ModelMultipleChoiceField(
        label=_('Arbeitsgruppen'), queryset=WorkGroup.objects.none(), required=False,
        widget=forms.CheckboxSelectMultiple(attrs={'class': 'form-check-input'}),
    )

    class Meta:
        model = PDBFile
        fields = ['personal', 'groups']

    def __init__(self, *args, user, **kwargs):
        super().__init__(*args, **kwargs)
        groups = WorkGroup.objects.filter(members=user)
        if self.instance.pk:
            groups = groups | self.instance.groups.all()
        self.fields['groups'].queryset = groups.distinct()

    def clean(self):
        cleaned = super().clean()
        if not cleaned.get('personal') and not cleaned.get('groups'):
            raise forms.ValidationError(_(
                'Die Datei muss für dich persönlich oder für mindestens eine Arbeitsgruppe '
                'verfügbar sein. Um sie zu entfernen, lösche sie.'
            ))
        return cleaned


class PDBFileAdminForm(forms.ModelForm):
    """Admin form of a PDB file with the same availability rule as ``AvailabilityForm``.

    Without it the admin could leave a file available to no one; such a file would stay in the database
    although no user can see it (outside the admin) or delete it.
    """

    class Meta:
        model = PDBFile
        fields = '__all__'

    def clean(self):
        cleaned = super().clean()
        # "personal" only keeps a file available as long as it has an uploader (see delete_if_orphaned)
        personal = cleaned.get('personal') and cleaned.get('uploaded_by')
        if not personal and not cleaned.get('groups'):
            raise forms.ValidationError(_(
                'Die Datei muss für den Uploader persönlich oder für mindestens eine Arbeitsgruppe verfügbar '
                'sein. Um sie zu entfernen, lösche sie.'
            ))
        return cleaned


class PDBUploadForm(AvailabilityForm):
    """New PDB file: either uploaded or fetched from RCSB by PDB ID, plus its availability.

    After validation ``cleaned_data['pdb_file']`` always holds the file, also when it was fetched by ID.
    """
    pdb_file = forms.FileField(
        label=_('PDB-Datei'), required=False,
        help_text=_('Unterstützt .pdb Dateien'),
    )
    pdb_id = forms.CharField(
        label=_('PDB-ID'), required=False, max_length=10,
        widget=forms.TextInput(attrs={
            'class': 'form-control font-monospace', 'placeholder': _('z. B. 1LYZ'),
            'autocomplete': 'off', 'spellcheck': 'false',
        }),
    )

    class Meta(AvailabilityForm.Meta):
        fields = ['pdb_file', 'pdb_id', 'personal', 'groups']

    def clean_pdb_file(self):
        f = self.cleaned_data['pdb_file']
        if not f:
            return f
        if not f.name.lower().endswith('.pdb'):
            raise forms.ValidationError(_('Nur .pdb Dateien sind erlaubt.'))
        if f.size > rcsb.MAX_SIZE:
            raise forms.ValidationError(_('Datei darf maximal 20 MB groß sein.'))
        return f

    def clean_pdb_id(self):
        pdb_id = self.cleaned_data['pdb_id']
        if not pdb_id:
            return pdb_id
        try:
            return rcsb.normalize_pdb_id(pdb_id)
        except rcsb.RCSBError as e:
            raise forms.ValidationError(str(e))

    def clean(self):
        cleaned = super().clean()
        if self.errors:
            return cleaned  # no RCSB request while the form is invalid anyway
        pdb_file, pdb_id = cleaned.get('pdb_file'), cleaned.get('pdb_id')
        if pdb_file and pdb_id:
            self.add_error('pdb_file', _('Bitte entweder eine Datei auswählen oder eine PDB-ID eingeben, nicht beides.'))
        elif not pdb_file and not pdb_id:
            self.add_error('pdb_file', _('Bitte eine Datei auswählen oder eine PDB-ID eingeben.'))
        elif pdb_id:
            try:
                content = rcsb.fetch_pdb(pdb_id)
            except rcsb.RCSBError as e:
                self.add_error('pdb_id', str(e))
            else:
                cleaned['pdb_file'] = SimpleUploadedFile(f'{pdb_id}.pdb', content)
        return cleaned


class NoteForm(forms.ModelForm):
    """Note with visibility: only for the author or additionally for selected work groups.

    Selectable are the user's groups in which the PDB file is shared (members of other groups could never
    read the note, because they cannot see the file), plus the groups the note is already shared with.
    ``visibility`` is not a model field; with ``PERSONAL`` the selected groups are discarded.

    Args:
        user: The author, whose groups are offered.
        pdb: The PDBFile the note belongs to.
    """
    PERSONAL, GROUPS = 'personal', 'groups'

    visibility = forms.ChoiceField(
        label=_('Sichtbar für'), initial=PERSONAL,
        choices=[(PERSONAL, _('Nur für mich')), (GROUPS, _('Für Arbeitsgruppen'))],
        widget=forms.RadioSelect(attrs={'class': 'form-check-input'}),
    )
    groups = forms.ModelMultipleChoiceField(
        label=_('Arbeitsgruppen'), queryset=WorkGroup.objects.none(), required=False,
        widget=forms.CheckboxSelectMultiple(attrs={'class': 'form-check-input'}),
    )

    class Meta:
        model = Note
        fields = ['text', 'groups']
        widgets = {
            'text': forms.Textarea(attrs={
                'class': 'form-control', 'rows': 3,
                'placeholder': _('Anmerkung schreiben …'),
            }),
        }

    def __init__(self, *args, user, pdb, **kwargs):
        super().__init__(*args, **kwargs)
        groups = WorkGroup.objects.filter(members=user, pdb_files=pdb)
        if self.instance.pk:
            current = self.instance.groups.all()
            groups = groups | current
            self.initial['visibility'] = self.GROUPS if current else self.PERSONAL
        self.fields['groups'].queryset = groups.distinct()

    def clean(self):
        cleaned = super().clean()
        if cleaned.get('visibility') == self.GROUPS:
            if not cleaned.get('groups') and 'groups' not in self.errors:
                raise forms.ValidationError(
                    _('Wähle mindestens eine Arbeitsgruppe aus, die die Anmerkung lesen darf.')
                )
        else:
            cleaned['groups'] = WorkGroup.objects.none()
        return cleaned
