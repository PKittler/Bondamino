from django.contrib import admin, messages
from django.contrib.auth import get_user_model
from django.contrib.auth.admin import UserAdmin
from django.db import transaction
from django.utils.translation import gettext, gettext_lazy as _, ngettext

from .forms import AdminLoginForm, InvitationUserForm, send_invitation

from .models import Language, Note, PDBFile, UserSettings, WorkGroup


@admin.register(Language)
class LanguageAdmin(admin.ModelAdmin):
    list_display = ('name', 'code')


@admin.register(WorkGroup)
class WorkGroupAdmin(admin.ModelAdmin):
    list_display = ('name', 'leader', 'member_count')
    search_fields = ('name', 'leader__username', 'members__username')
    filter_horizontal = ('members',)

    @admin.display(description=_('Mitglieder'))
    def member_count(self, obj):
        return obj.members.count()

    def save_related(self, request, form, formsets, change):
        super().save_related(request, form, formsets, change)
        # Der Gruppenleiter ist immer auch Mitglied der Gruppe
        form.instance.members.add(form.instance.leader)


@admin.register(PDBFile)
class PDBFileAdmin(admin.ModelAdmin):
    list_display = ('name', 'uploaded_by', 'uploaded_at', 'personal')
    list_filter = ('groups', 'personal')
    search_fields = ('name', 'uploaded_by__username')
    filter_horizontal = ('groups',)
    readonly_fields = ('uploaded_at',)


@admin.register(Note)
class NoteAdmin(admin.ModelAdmin):
    list_display = ('__str__', 'pdb', 'author', 'created_at', 'edited_at')
    list_filter = ('groups',)
    search_fields = ('text', 'pdb__name', 'author__username')
    filter_horizontal = ('groups',)
    readonly_fields = ('created_at', 'edited_at')


class WorkGroupMembershipInline(admin.TabularInline):
    model = WorkGroup.members.through
    verbose_name = _('Arbeitsgruppe')
    verbose_name_plural = _('Arbeitsgruppen')
    extra = 0


class UserSettingsInline(admin.StackedInline):
    model = UserSettings
    can_delete = False
    extra = 0


# Begrenzung der Anmeldeversuche auch im Admin
admin.site.login_form = AdminLoginForm


User = get_user_model()
admin.site.unregister(User)


@admin.register(User)
class BondaminoUserAdmin(UserAdmin):
    """Neue User bekommen kein Passwort, sondern eine Einladung per E-Mail mit einem Link, unter
    dem sie es selbst festlegen."""
    inlines = [*UserAdmin.inlines, WorkGroupMembershipInline, UserSettingsInline]
    add_form = InvitationUserForm
    add_form_template = 'admin/change_form.html'
    add_fieldsets = (
        (None, {'classes': ('wide',), 'fields': ('username', 'email', 'first_name', 'last_name')}),
    )
    list_display = (*UserAdmin.list_display, 'password_set')
    actions = ['resend_invitation']

    def get_inlines(self, request, obj):
        # Beim Anlegen nur die Arbeitsgruppen; Einstellungen legt der User selbst fest
        return [WorkGroupMembershipInline] if obj is None else self.inlines

    @admin.display(description=_('Passwort festgelegt'), boolean=True)
    def password_set(self, obj):
        return obj.has_usable_password()

    def _invite(self, request, user):
        try:
            send_invitation(request, user)
        except Exception as error:  # SMTP-Fehler u. ä.: User bleibt angelegt
            self.message_user(request, gettext('Die Einladung an %(email)s konnte nicht gesendet werden: %(error)s')
                              % {'email': user.email, 'error': error}, messages.ERROR)
            return False
        return True

    def save_model(self, request, obj, form, change):
        super().save_model(request, obj, form, change)
        if not change:
            # Erst nach dem Speichern aller Daten senden, damit der Link gültig ist
            transaction.on_commit(lambda: self._invite(request, obj) and self.message_user(
                request, gettext('Die Einladung wurde an %(email)s gesendet.') % {'email': obj.email}, messages.SUCCESS))

    @admin.action(description=_('Einladung erneut senden'))
    def resend_invitation(self, request, queryset):
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
