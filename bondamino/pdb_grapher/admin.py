from django.contrib import admin
from django.contrib.auth import get_user_model
from django.contrib.auth.admin import UserAdmin
from django.utils.translation import gettext_lazy as _

from .models import Language, Note, PDBFile, WorkGroup


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


User = get_user_model()
admin.site.unregister(User)


@admin.register(User)
class BondaminoUserAdmin(UserAdmin):
    inlines = [*UserAdmin.inlines, WorkGroupMembershipInline]
