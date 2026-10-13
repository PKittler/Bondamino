from django.conf import settings
from django.db.models.signals import post_delete, pre_delete
from django.dispatch import receiver

from .models import PDBFile, WorkGroup


@receiver(post_delete, sender=PDBFile)
def delete_stored_file(sender, instance, **kwargs):
    """Entfernt die gespeicherte Datei vom Datenträger, sobald der Eintrag gelöscht ist."""
    if instance.file:
        instance.file.delete(save=False)


@receiver(pre_delete, sender=WorkGroup)
def remember_group_files(sender, instance, **kwargs):
    instance._orphan_candidates = list(instance.pdb_files.all())
    instance._orphan_note_candidates = list(instance.notes.all())


@receiver(post_delete, sender=WorkGroup)
def delete_files_orphaned_by_group(sender, instance, **kwargs):
    """Wird eine Gruppe gelöscht, verschwinden Dateien, die nur in ihr verfügbar waren,
    und Anmerkungen ohne Verfasser, die nur für sie lesbar waren."""
    for pdb in getattr(instance, '_orphan_candidates', []):
        pdb.delete_if_orphaned()
    for note in getattr(instance, '_orphan_note_candidates', []):
        note.delete_if_orphaned()


@receiver(pre_delete, sender=settings.AUTH_USER_MODEL)
def delete_personal_files_of_user(sender, instance, **kwargs):
    """Wird ein User gelöscht, entfallen seine persönlichen Dateien; in Gruppen
    freigegebene Dateien bleiben (ohne Uploader) erhalten."""
    for pdb in instance.pdb_uploads.all():
        if pdb.groups.exists():
            if pdb.personal:
                pdb.personal = False
                pdb.save(update_fields=['personal'])
        else:
            pdb.delete()


@receiver(pre_delete, sender=settings.AUTH_USER_MODEL)
def delete_personal_notes_of_user(sender, instance, **kwargs):
    """Wird ein User gelöscht, entfallen seine persönlichen Anmerkungen; für Gruppen
    lesbare Anmerkungen bleiben (ohne Verfasser) erhalten."""
    instance.notes.filter(groups__isnull=True).delete()
