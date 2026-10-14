"""Clean-up on deletion (connected in ``apps.PdbGrapherConfig.ready``).

These receivers keep the visibility rules of ``models.py`` consistent when objects are deleted, also in the admin
and in bulk deletes of querysets:

* A deleted PDBFile also loses its stored file on disk.
* A deleted work group takes along the files and notes that were available only through it.
* A deleted user takes along the personal files and notes; anything shared with a group stays, without uploader
  or author.
"""
from django.conf import settings
from django.db.models.signals import post_delete, pre_delete
from django.dispatch import receiver

from .models import PDBFile, WorkGroup


@receiver(post_delete, sender=PDBFile)
def delete_stored_file(sender, instance, **kwargs):
    """Removes the stored file from disk once the database row is deleted."""
    if instance.file:
        instance.file.delete(save=False)


@receiver(pre_delete, sender=WorkGroup)
def remember_group_files(sender, instance, **kwargs):
    """Remembers the group's files and notes before the group is deleted.

    Must happen in pre_delete: after the deletion the M2M rows are gone and the candidates can no longer be
    found. ``delete_files_orphaned_by_group`` checks them afterwards.
    """
    instance._orphan_candidates = list(instance.pdb_files.all())
    instance._orphan_note_candidates = list(instance.notes.all())


@receiver(post_delete, sender=WorkGroup)
def delete_files_orphaned_by_group(sender, instance, **kwargs):
    """Deletes the files and notes orphaned by a deleted group.

    These are files that were available only in the group and notes without author that were readable only by
    it (see ``delete_if_orphaned`` in ``models.py``).
    """
    for pdb in getattr(instance, '_orphan_candidates', []):
        pdb.delete_if_orphaned()
    for note in getattr(instance, '_orphan_note_candidates', []):
        note.delete_if_orphaned()


@receiver(pre_delete, sender=settings.AUTH_USER_MODEL)
def delete_personal_files_of_user(sender, instance, **kwargs):
    """When a user is deleted, deletes the user's files that are not shared with any group.

    Files shared with groups are kept; their ``personal`` flag is cleared and ``uploaded_by`` becomes ``None``
    through ``on_delete=SET_NULL``.
    """
    for pdb in instance.pdb_uploads.all():
        if pdb.groups.exists():
            if pdb.personal:
                pdb.personal = False
                pdb.save(update_fields=['personal'])
        else:
            pdb.delete()


@receiver(pre_delete, sender=settings.AUTH_USER_MODEL)
def delete_personal_notes_of_user(sender, instance, **kwargs):
    """When a user is deleted, deletes the user's personal notes.

    Notes readable by groups are kept; their ``author`` becomes ``None`` through ``on_delete=SET_NULL``.
    """
    instance.notes.filter(groups__isnull=True).delete()
