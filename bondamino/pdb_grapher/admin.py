"""Django admin configuration of Bondamino.

The admin is where languages, work groups (members and leader) and users are maintained. Users are created
without a password and receive an invitation email (``BondaminoUserAdmin``); the admin login is throttled like the
site login (``forms.AdminLoginForm``).
"""
from django.contrib import admin, messages
from django.contrib.auth import get_user_model
from django.contrib.auth.admin import UserAdmin
from django.db import transaction
from django.utils.translation import gettext, gettext_lazy as _, ngettext

from .forms import AdminLoginForm, InvitationUserForm, PDBFileAdminForm, send_invitation

from .models import Language, Note, PDBFile, UserSettings, WorkGroup


@admin.register(Language)
class LanguageAdmin(admin.ModelAdmin):
    """Languages offered in the navbar; only these can be activated."""
    list_display = ('name', 'code')


@admin.register(WorkGroup)
class WorkGroupAdmin(admin.ModelAdmin):
    """Work groups with leader and members."""
    list_display = ('name', 'leader', 'member_count')
    search_fields = ('name', 'leader__username', 'members__username')
    filter_horizontal = ('members',)

    @admin.display(description=_('Mitglieder'))
    def member_count(self, obj):
        """Number of members, shown in the change list."""
        return obj.members.count()

    def save_related(self, request, form, formsets, change):
        """Saves the members and makes sure the leader is one of them."""
        super().save_related(request, form, formsets, change)
        # The leader is always a member too; added after the M2M save, which would otherwise overwrite it
        form.instance.members.add(form.instance.leader)


@admin.register(PDBFile)
class PDBFileAdmin(admin.ModelAdmin):
    """PDB files; the admin can see and change all of them, independent of visibility."""
    form = PDBFileAdminForm
    list_display = ('name', 'uploaded_by', 'uploaded_at', 'personal')
    list_filter = ('groups', 'personal')
    search_fields = ('name', 'uploaded_by__username')
    filter_horizontal = ('groups',)
    readonly_fields = ('uploaded_at',)


@admin.register(Note)
class NoteAdmin(admin.ModelAdmin):
    """Notes; the admin can see and change all of them, independent of visibility."""
    list_display = ('__str__', 'pdb', 'author', 'created_at', 'edited_at')
    list_filter = ('groups',)
    search_fields = ('text', 'pdb__name', 'author__username')
    filter_horizontal = ('groups',)
    readonly_fields = ('created_at', 'edited_at')


class WorkGroupMembershipInline(admin.TabularInline):
    """Group memberships on the user page (edits the M2M table of ``WorkGroup.members``)."""
    model = WorkGroup.members.through
    verbose_name = _('Arbeitsgruppe')
    verbose_name_plural = _('Arbeitsgruppen')
    extra = 0


class UserSettingsInline(admin.StackedInline):
    """The user's personal settings on the user page."""
    model = UserSettings
    can_delete = False
    extra = 0


# Throttle login attempts in the admin as well
admin.site.login_form = AdminLoginForm


# Replace Django's default user admin with the invitation-based one below
User = get_user_model()
admin.site.unregister(User)


@admin.register(User)
class BondaminoUserAdmin(UserAdmin):
    """User admin with invitations.

    New users get no password but an invitation email with a link where they set it themselves. The action
    "resend invitation" sends a new link to users who have not set a password yet.
    """
    inlines = [*UserAdmin.inlines, WorkGroupMembershipInline, UserSettingsInline]
    add_form = InvitationUserForm
    add_form_template = 'admin/change_form.html'
    add_fieldsets = (
        (None, {'classes': ('wide',), 'fields': ('username', 'email', 'first_name', 'last_name')}),
    )
    list_display = (*UserAdmin.list_display, 'password_set')
    actions = ['resend_invitation']

    def get_inlines(self, request, obj):
        # When creating a user only the work groups; the user chooses the settings
        return [WorkGroupMembershipInline] if obj is None else self.inlines

    @admin.display(description=_('Passwort festgelegt'), boolean=True)
    def password_set(self, obj):
        """Whether the user has set a password, i.e. accepted the invitation."""
        return obj.has_usable_password()

    def _invite(self, request, user):
        """Sends the invitation and reports a failure as an admin message.

        Returns:
            True if the email was sent.
        """
        try:
            send_invitation(request, user)
        except Exception as error:  # SMTP errors etc.: the user is still created
            self.message_user(request, gettext('Die Einladung an %(email)s konnte nicht gesendet werden: %(error)s')
                              % {'email': user.email, 'error': error}, messages.ERROR)
            return False
        return True

    def save_model(self, request, obj, form, change):
        """Saves the user and, for a new user, sends the invitation."""
        super().save_model(request, obj, form, change)
        if not change:
            # Send only after the transaction (incl. inlines) is committed, so the user exists when the link is used
            # and a rollback does not leave an invitation for a user that was never created
            transaction.on_commit(lambda: self._invite(request, obj) and self.message_user(
                request, gettext('Die Einladung wurde an %(email)s gesendet.') % {'email': obj.email}, messages.SUCCESS))

    @admin.action(description=_('Einladung erneut senden'))
    def resend_invitation(self, request, queryset):
        """Admin action: resends the invitation to the selected active users without a password and with email.

        Users who do not qualify are skipped and counted in a warning.
        """
        invited = [u for u in queryset if not u.has_usable_password() and u.email and u.is_active]
        sent = sum(1 for user in invited if self._invite(request, user))
        if sent:
            self.message_user(request, ngettext('%(count)d Einladung wurde gesendet.',
                                                '%(count)d Einladungen wurden gesendet.', sent) % {'count': sent},
                              messages.SUCCESS)
        skipped = len(queryset) - len(invited)
        if skipped:
            self.message_user(request, ngettext(
                '%(count)d User hat bereits ein Passwort, keine E-Mail-Adresse oder ist inaktiv.',
                '%(count)d User haben bereits ein Passwort, keine E-Mail-Adresse oder sind inaktiv.', skipped,
            ) % {'count': skipped}, messages.WARNING)
