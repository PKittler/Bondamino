from django.conf import settings
from django.db import models
from django.db.models import Q
from django.utils.translation import gettext_lazy as _


class Language(models.Model):
    """Sprache, die User in der Navbar auswählen können. Wird nur im Admin gepflegt; die
    Texte stammen aus den Übersetzungsdateien in pdb_grapher/locale/."""
    code = models.CharField(_('Sprache'), max_length=15, unique=True, choices=settings.LANGUAGES)
    name = models.CharField(
        _('Anzeigename'), max_length=50,
        help_text=_('Name in der Sprachauswahl, in der Sprache selbst, z. B. „English“.'),
    )

    class Meta:
        verbose_name = _('Sprache')
        verbose_name_plural = _('Sprachen')
        ordering = ['name']

    def __str__(self):
        return self.name


class UserSettings(models.Model):
    """Persönliche Einstellungen eines Users. Wird beim ersten Speichern angelegt; bis dahin
    gelten Sprache aus Cookie bzw. Browser und die helle Darstellung."""
    THEME_LIGHT = 'light'
    THEME_DARK = 'dark'
    THEME_CHOICES = [(THEME_LIGHT, _('Hell')), (THEME_DARK, _('Dunkel'))]

    user = models.OneToOneField(
        settings.AUTH_USER_MODEL, verbose_name=_('User'), on_delete=models.CASCADE, related_name='user_settings',
    )
    # Code einer im Admin angelegten Sprache; leer = Cookie bzw. Browser
    language = models.CharField(_('Sprache'), max_length=15, blank=True)
    theme = models.CharField(_('Darstellung'), max_length=10, choices=THEME_CHOICES, default=THEME_LIGHT)

    class Meta:
        verbose_name = _('Einstellungen')
        verbose_name_plural = _('Einstellungen')

    def __str__(self):
        return str(self.user)

    @classmethod
    def of(cls, user):
        """Gespeicherte Einstellungen des Users oder ungespeicherte Standardwerte."""
        return cls.objects.filter(user=user).first() or cls(user=user)


class WorkGroup(models.Model):
    """Arbeitsgruppe eines Instituts/einer Universität. Mitglieder werden vom Admin
    zugewiesen; der Gruppenleiter darf PDB-Dateien aus der Gruppe entfernen."""
    name = models.CharField(_('Name'), max_length=150, unique=True)
    description = models.TextField(_('Beschreibung'), blank=True)
    leader = models.ForeignKey(
        settings.AUTH_USER_MODEL, verbose_name=_('Gruppenleiter'),
        on_delete=models.PROTECT, related_name='led_workgroups',
    )
    members = models.ManyToManyField(
        settings.AUTH_USER_MODEL, verbose_name=_('Mitglieder'),
        related_name='workgroups', blank=True,
    )

    class Meta:
        verbose_name = _('Arbeitsgruppe')
        verbose_name_plural = _('Arbeitsgruppen')
        ordering = ['name']

    def __str__(self):
        return self.name


class PDBFileQuerySet(models.QuerySet):
    def visible_to(self, user):
        """Dateien, die der User hochgeladen hat oder die für eine seiner Gruppen
        freigegeben sind."""
        return self.filter(Q(uploaded_by=user) | Q(groups__members=user)).distinct()


class PDBFile(models.Model):
    file = models.FileField(_('Datei'), upload_to='pdb_uploads/%Y/%m/')
    name = models.CharField(_('Dateiname'), max_length=255)
    uploaded_by = models.ForeignKey(
        settings.AUTH_USER_MODEL, verbose_name=_('Hochgeladen von'),
        on_delete=models.SET_NULL, null=True, related_name='pdb_uploads',
    )
    uploaded_at = models.DateTimeField(_('Hochgeladen am'), auto_now_add=True)
    # Verfügbarkeit: persönlich für den Uploader und/oder für Arbeitsgruppen
    personal = models.BooleanField(
        _('Persönlich'), default=False,
        help_text=_('Die Datei ist für den Uploader persönlich verfügbar.'),
    )
    groups = models.ManyToManyField(
        WorkGroup, verbose_name=_('Arbeitsgruppen'), related_name='pdb_files', blank=True,
    )

    objects = PDBFileQuerySet.as_manager()

    class Meta:
        verbose_name = _('PDB-Datei')
        verbose_name_plural = _('PDB-Dateien')
        ordering = ['-uploaded_at']

    def __str__(self):
        return self.name

    def can_manage(self, user):
        """Nur der Uploader darf die Datei löschen und ihre Verfügbarkeit bearbeiten."""
        return self.uploaded_by_id is not None and self.uploaded_by_id == user.pk

    def delete_if_orphaned(self):
        """Löscht die Datei, wenn sie weder persönlich noch für eine Gruppe verfügbar
        ist. Gibt True zurück, wenn gelöscht wurde."""
        if (self.personal and self.uploaded_by_id) or self.groups.exists():
            return False
        self.delete()
        return True


class NoteQuerySet(models.QuerySet):
    def visible_to(self, user):
        """Eigene Anmerkungen und solche, die für eine Gruppe des Users lesbar sind."""
        return self.filter(Q(author=user) | Q(groups__members=user)).distinct()


class Note(models.Model):
    """Anmerkung eines Users zu einer PDB-Datei. Ohne Arbeitsgruppen ist sie nur für
    den Verfasser sichtbar, sonst zusätzlich für die Mitglieder der gewählten Gruppen."""
    pdb = models.ForeignKey(
        PDBFile, verbose_name=_('PDB-Datei'), on_delete=models.CASCADE, related_name='notes',
    )
    author = models.ForeignKey(
        settings.AUTH_USER_MODEL, verbose_name=_('Verfasser'),
        on_delete=models.SET_NULL, null=True, related_name='notes',
    )
    text = models.TextField(_('Text'), max_length=10000)
    groups = models.ManyToManyField(
        WorkGroup, verbose_name=_('Lesbar für'), related_name='notes', blank=True,
    )
    created_at = models.DateTimeField(_('Erstellt am'), auto_now_add=True)
    edited_at = models.DateTimeField(_('Bearbeitet am'), null=True, blank=True)

    objects = NoteQuerySet.as_manager()

    class Meta:
        verbose_name = _('Anmerkung')
        verbose_name_plural = _('Anmerkungen')
        ordering = ['-created_at']

    def __str__(self):
        return self.text if len(self.text) <= 50 else f'{self.text[:49]}…'

    def editable_by(self, user):
        """Nur der Verfasser darf die Anmerkung bearbeiten und löschen."""
        return self.author_id is not None and self.author_id == user.pk

    def delete_if_orphaned(self):
        """Löscht die Anmerkung, wenn sie weder einen Verfasser hat noch für eine Gruppe
        lesbar ist. Gibt True zurück, wenn gelöscht wurde."""
        if self.author_id or self.groups.exists():
            return False
        self.delete()
        return True
