"""Data model of Bondamino.

Overview:

* ``Language`` and ``UserSettings`` hold the UI configuration. Only languages created in the admin can be
  selected; the actual texts come from the translation files in ``pdb_grapher/locale/``.
* ``WorkGroup`` groups users (members are assigned in the admin). Every group has a leader who may remove
  files and notes from the group.
* ``PDBFile`` is an uploaded (or RCSB-fetched) structure file. ``Note`` is a text annotation on a PDBFile.

Visibility and ownership rules:

* A PDBFile is visible to its uploader and to all members of the groups it is shared with
  (``PDBFile.objects.visible_to``). Only the uploader may change its availability or delete it
  (``PDBFile.can_manage``). The ``personal`` flag only controls whether the file stays available to the uploader
  once it is no longer shared with any group; a file that is neither personal nor shared is deleted
  (``PDBFile.delete_if_orphaned``).
* A Note is visible to its author and to the members of the groups it is shared with
  (``Note.objects.visible_to``). Only the author may edit or delete it (``Note.editable_by``).
* Group leaders may remove a file or note from their group; if it is then orphaned it is deleted
  (see ``views.remove_from_group`` and ``views.note_remove_from_group``).
* Deleting a user or a group triggers further clean-up, see ``signals.py``.

Every query in the views goes through ``visible_to`` so users never see files or notes they are not entitled to.
"""
from django.conf import settings
from django.db import models
from django.db.models import Q
from django.utils.translation import gettext_lazy as _


class Language(models.Model):
    """A UI language that users can select in the navbar.

    Maintained in the admin only. ``code`` must be one of ``settings.LANGUAGES``; the translated texts come
    from the translation files in ``pdb_grapher/locale/``. Languages not created here are never activated
    (see ``middleware.LanguageMiddleware``).

    Attributes:
        code: Language code, e.g. ``'en'``.
        name: Name shown in the language selector, written in the language itself.
    """
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
    """Personal settings of a user.

    The row is created the first time the user saves a setting. Until then the language is taken from the
    cookie or the browser and the light theme is used (see ``of``).

    Attributes:
        user: The user these settings belong to.
        language: Code of a ``Language``; empty means "use cookie or browser".
        theme: ``THEME_LIGHT`` or ``THEME_DARK``.
        time_zone: IANA time zone name (e.g. ``Europe/Berlin``) in which dates and times are shown to the
            user; empty means the server's ``TIME_ZONE``. Activated by ``middleware.LanguageMiddleware``.
    """
    THEME_LIGHT = 'light'
    THEME_DARK = 'dark'
    THEME_CHOICES = [(THEME_LIGHT, _('Hell')), (THEME_DARK, _('Dunkel'))]

    user = models.OneToOneField(
        settings.AUTH_USER_MODEL, verbose_name=_('User'), on_delete=models.CASCADE, related_name='user_settings',
    )
    # Code of a language created in the admin; empty = cookie or browser
    language = models.CharField(_('Sprache'), max_length=15, blank=True)
    theme = models.CharField(_('Darstellung'), max_length=10, choices=THEME_CHOICES, default=THEME_LIGHT)
    # IANA name; empty = settings.TIME_ZONE
    time_zone = models.CharField(_('Zeitzone'), max_length=64, blank=True)

    class Meta:
        verbose_name = _('Einstellungen')
        verbose_name_plural = _('Einstellungen')

    def __str__(self):
        return str(self.user)

    @classmethod
    def of(cls, user):
        """Returns the stored settings of a user, or unsaved defaults.

        Args:
            user: An authenticated user.

        Returns:
            The user's ``UserSettings``; a new, unsaved instance with default values if none is stored yet.
        """
        return cls.objects.filter(user=user).first() or cls(user=user)


class WorkGroup(models.Model):
    """Work group of an institute or university.

    Members are assigned in the admin; the admin also adds the leader as a member when the group is saved
    (``admin.WorkGroupAdmin.save_related``). Files and notes can be shared with groups; the leader may remove
    PDB files and notes from the group.

    Attributes:
        name: Unique group name.
        description: Optional free text.
        leader: Group leader; a user cannot be deleted while leading a group (``PROTECT``).
        members: Users who see the files and notes shared with the group.
    """
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
    """QuerySet (and manager) of ``PDBFile`` with the central visibility filter."""

    def visible_to(self, user):
        """Filters the files a user may see.

        Args:
            user: The requesting user.

        Returns:
            Files uploaded by the user or shared with one of the user's groups, without duplicates.
        """
        return self.filter(Q(uploaded_by=user) | Q(groups__members=user)).distinct()


class PDBFile(models.Model):
    """An uploaded PDB structure file.

    The file itself is stored under ``MEDIA_ROOT/pdb_uploads/`` and removed from disk when the row is deleted
    (``signals.delete_stored_file``). Parsing happens on demand in ``utils.py``; nothing derived from the file
    is stored in the database.

    Attributes:
        file: The stored file.
        name: Original file name, shown in the UI and used for downloads.
        uploaded_by: Uploader; ``None`` after the uploader was deleted (the file then lives on in its groups).
        uploaded_at: Upload timestamp.
        personal: Whether the file stays available to the uploader independently of any group.
        groups: Work groups the file is shared with.
    """
    file =models.FileField(_('Datei'), upload_to='pdb_uploads/%Y/%m/')
    name = models.CharField(_('Dateiname'), max_length=255)
    uploaded_by = models.ForeignKey(
        settings.AUTH_USER_MODEL, verbose_name=_('Hochgeladen von'),
        on_delete=models.SET_NULL, null=True, related_name='pdb_uploads',
    )
    uploaded_at = models.DateTimeField(_('Hochgeladen am'), auto_now_add=True)
    # Availability: personal for the uploader and/or for work groups
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
        """Checks whether a user may delete the file and edit its availability.

        Args:
            user: The requesting user.

        Returns:
            True only for the uploader.
        """
        return self.uploaded_by_id is not None and self.uploaded_by_id == user.pk

    def delete_if_orphaned(self):
        """Deletes the file if it is neither personal (with an existing uploader) nor shared with a group.

        Returns:
            True if the file was deleted.
        """
        if (self.personal and self.uploaded_by_id) or self.groups.exists():
            return False
        self.delete()
        return True


class NoteQuerySet(models.QuerySet):
    """QuerySet (and manager) of ``Note`` with the central visibility filter."""

    def visible_to(self, user):
        """Filters the notes a user may read.

        Args:
            user: The requesting user.

        Returns:
            Notes written by the user or shared with one of the user's groups, without duplicates.
        """
        return self.filter(Q(author=user) | Q(groups__members=user)).distinct()


class Note(models.Model):
    """A user's note on a PDB file.

    Without work groups it is visible only to its author, otherwise also to the members of the selected groups.

    Attributes:
        pdb: The annotated file; notes are deleted together with it.
        author: Author; ``None`` after the author was deleted (the note then lives on in its groups).
        text: The note text.
        groups: Work groups that may read the note.
        created_at: Creation timestamp.
        edited_at: Timestamp of the last edit, ``None`` if never edited.
    """
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
        """Checks whether a user may edit and delete the note.

        Args:
            user: The requesting user.

        Returns:
            True only for the author.
        """
        return self.author_id is not None and self.author_id == user.pk

    def delete_if_orphaned(self):
        """Deletes the note if it has no author and is not shared with any group.

        Returns:
            True if the note was deleted.
        """
        if self.author_id or self.groups.exists():
            return False
        self.delete()
        return True
