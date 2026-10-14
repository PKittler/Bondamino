"""Test suite of the pdb_grapher app.

Run it from the Django project directory:

    cd bondamino && python manage.py test pdb_grapher

A single class or test can be selected with its dotted path, e.g.
``python manage.py test pdb_grapher.tests.NoteTests.test_visibility``.

Organisation:

    Shared fixtures:
        BaseTestCase: creates four users and two work groups (see its docstring) and provides
            ``make_pdb()``. Uploaded files go to a temporary MEDIA_ROOT that is removed afterwards.
            Test classes that need users, groups or uploaded files derive from it; pure parser tests
            derive from ``django.test.TestCase`` directly and write PDB files to temporary files via
            their own ``write_pdb()`` helper.

    Access control, uploads and pages (views, using the Django test client):
        PermissionTests: who may see, download, edit, delete or remove PDB files from groups.
        NoteTests: personal and group notes and their permissions.
        IndexTests: file list on the start page (search, latest note, upload dialog).
        RCSBImportTests: import by PDB ID from RCSB (network access is mocked).
        LanguageTests: language selection, browser language, JavaScript catalog.
        AccountTests: account page, password change, stored language and theme.
        InvitationTests: users created in the admin without a password receive an invitation email.
        PasswordResetTests: "forgot password" flow.
        ThrottleTests: limits for failed logins and password reset requests (``throttle`` module).
        OutputSafetyTests: escaping of graph data and hiding of internal error details.
        ConformationViewTests, StructureInfoViewTests, AnglesTabTests: detail page elements.

    Configuration and installation:
        ConfigTests: helpers in ``bondamino.config`` (DATABASE_URL, .env file, environment variables).
        SetupCommandTests: the interactive ``manage.py setup`` command (input is mocked).

    PDB parsing and graph building (``utils`` module):
        ChainGapTests: no peptide bond across chain gaps.
        ConformationParsingTests: insertion codes and alternate locations (altlocs).
        ModifiedResidueTests: modified residues (MODRES, CONECT) and ligands.
        ReferenceValueTests: reference values for bond lengths and angles, header parsing.
        DisulfideTests: disulfide bridges by distance and from SSBOND records.
        RamachandranTests: phi/psi angles, residue types and Ramachandran regions.

PDB line helpers:
    PDB files are fixed-column text files (1-based columns: record name 1-6, atom serial 7-11,
    atom name 13-16, altLoc 17, residue name 18-20, chain 22, residue number 23-26,
    insertion code 27, x/y/z 31-38/39-46/47-54, occupancy 55-60, element 77-78). The helpers
    ``atom_line()``, ``pdb_line()``, ``xyz_line()`` and ``het_line()`` format single ATOM/HETATM lines;
    ``residue_lines()``, ``cys_lines()`` and ``peptide_lines()`` build whole residues or chains.
    Module-level constants (``PDB_CONTENT``, ``ALTLOC_LINES``, ``MSE_LINES``, ``HEADER_LINES``) are
    shared fixtures built from them.

Note: German strings asserted with ``assertContains`` are the UI source texts that are translated via
the locale files; they must match the templates and code exactly.
"""
import datetime
import math
import re
import os
import shutil
import tempfile
import urllib.error
from pathlib import Path
from unittest import mock

from django.conf import settings
from django.contrib.auth import get_user_model
from django.core import mail
from django.core.cache import cache
from django.core.management import call_command
from django.core.files.base import ContentFile
from django.core.files.uploadedfile import SimpleUploadedFile
from django.db import DatabaseError
from django.test import TestCase, override_settings
from django.urls import reverse

from bondamino.config import env, env_bool, env_list, parse_database_url, read_env_file

from . import throttle
from .forms import NoteForm, PDBFileAdminForm
from .models import Language, Note, PDBFile, UserSettings, WorkGroup
from .utils import (
    SIDECHAIN_BOND_REFERENCES, SIDECHAIN_BONDS, build_average_graph, build_bond_ensemble, build_graph, list_altlocs,
    list_models, parse_header, parse_pdb, ramachandran_region, ramachandran_value,
)

# Minimal PDB file: the backbone atoms N, CA, C, O of a single alanine (residue 1, chain A)
PDB_CONTENT = (
    b"ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N\n"
    b"ATOM      2  CA  ALA A   1      11.639   6.071  -5.147  1.00  0.00           C\n"
    b"ATOM      3  C   ALA A   1      13.161   6.195  -5.161  1.00  0.00           C\n"
    b"ATOM      4  O   ALA A   1      13.773   5.880  -6.191  1.00  0.00           O\n"
    b"END\n"
)

# Uploaded files are stored in a temporary directory instead of the real media directory
MEDIA_ROOT = tempfile.mkdtemp()
User = get_user_model()


@override_settings(MEDIA_ROOT=MEDIA_ROOT)
class BaseTestCase(TestCase):
    """Base class with users, work groups and helpers shared by the view tests.

    Fixtures created in ``setUp()`` (all users have the password ``'pw'``):
        uploader: member of groups A and B; uploads the files created by ``make_pdb()``.
        leader: leader and member of group A.
        member: plain member of group A.
        outsider: leader and member of group B, not in group A.
        group_a: work group "A" (leader, uploader, member).
        group_b: work group "B" (outsider, uploader).
    """

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()
        # Remove the temporary MEDIA_ROOT including all uploaded test files
        shutil.rmtree(MEDIA_ROOT, ignore_errors=True)

    def setUp(self):
        self.uploader = User.objects.create_user('uploader', password='pw')
        self.leader = User.objects.create_user('leader', password='pw')
        self.member = User.objects.create_user('member', password='pw')
        self.outsider = User.objects.create_user('outsider', password='pw')
        self.group_a = WorkGroup.objects.create(name='A', leader=self.leader)
        self.group_a.members.add(self.leader, self.uploader, self.member)
        self.group_b = WorkGroup.objects.create(name='B', leader=self.outsider)
        self.group_b.members.add(self.outsider, self.uploader)

    def make_pdb(self, personal=False, groups=()):
        """Creates a PDB file uploaded by ``self.uploader`` with the content ``PDB_CONTENT``.

        Args:
            personal: Whether the file is available in the uploader's personal area.
            groups: Work groups the file is shared with.

        Returns:
            The saved ``PDBFile`` instance.
        """
        pdb = PDBFile(name='test.pdb', uploaded_by=self.uploader, personal=personal)
        pdb.file.save('test.pdb', ContentFile(PDB_CONTENT), save=False)
        pdb.save()
        pdb.groups.set(groups)
        return pdb


class PermissionTests(BaseTestCase):
    """Upload, visibility, download, editing and deletion of PDB files and their group sharing."""

    def test_login_required(self):
        response = self.client.get(reverse('pdb_grapher:index'))
        self.assertRedirects(response, f"{reverse('login')}?next=/")

    def test_upload_stores_file_with_availability(self):
        self.client.login(username='uploader', password='pw')
        response = self.client.post(reverse('pdb_grapher:index'), {
            'pdb_file': SimpleUploadedFile('x.pdb', PDB_CONTENT),
            'groups': [self.group_a.pk],
        })
        pdb = PDBFile.objects.get()
        self.assertRedirects(response, reverse('pdb_grapher:detail', args=[pdb.pk]))
        self.assertEqual(pdb.uploaded_by, self.uploader)
        self.assertFalse(pdb.personal)
        self.assertEqual(list(pdb.groups.all()), [self.group_a])

    def test_upload_requires_availability(self):
        # Neither "personal" nor a group is selected, so the file would be visible to no one
        self.client.login(username='uploader', password='pw')
        self.client.post(reverse('pdb_grapher:index'), {'pdb_file': SimpleUploadedFile('x.pdb', PDB_CONTENT)})
        self.assertFalse(PDBFile.objects.exists())

    def test_upload_rejects_foreign_group(self):
        self.client.login(username='member', password='pw')
        self.client.post(reverse('pdb_grapher:index'), {
            'pdb_file': SimpleUploadedFile('x.pdb', PDB_CONTENT), 'groups': [self.group_b.pk],
        })
        self.assertFalse(PDBFile.objects.exists())

    def test_visibility(self):
        pdb = self.make_pdb(groups=[self.group_a])
        # Detail page, model data and download all answer 404 (not 403) for users who may not see the file
        for username, status in (('uploader', 200), ('member', 200), ('outsider', 404)):
            self.client.login(username=username, password='pw')
            self.assertEqual(self.client.get(reverse('pdb_grapher:detail', args=[pdb.pk])).status_code, status)
            self.assertEqual(self.client.get(reverse('pdb_grapher:model_data', args=[pdb.pk])).status_code, status)
            self.assertEqual(self.client.get(reverse('pdb_grapher:download', args=[pdb.pk])).status_code, status)

    def test_download_returns_file_as_attachment(self):
        pdb = self.make_pdb(groups=[self.group_a])
        self.client.login(username='member', password='pw')
        response = self.client.get(reverse('pdb_grapher:download', args=[pdb.pk]))
        self.assertEqual(response['Content-Disposition'], 'attachment; filename="test.pdb"')
        self.assertEqual(b''.join(response.streaming_content), PDB_CONTENT)

    def test_single_model_shows_disabled_model_select(self):
        pdb = self.make_pdb(personal=True)
        self.client.login(username='uploader', password='pw')
        response = self.client.get(reverse('pdb_grapher:detail', args=[pdb.pk]))
        self.assertContains(response, 'id="model-select"')
        self.assertContains(response, '<option>Model 1</option>', html=True)
        self.assertNotContains(response, 'model-select-wrap')
        self.assertContains(response, '<dd id="stat-models">1</dd>', html=True)

    def test_model_count_in_stats(self):
        pdb = self.make_pdb(personal=True)
        # NMR-style file with two MODEL ... ENDMDL blocks
        pdb.file.save('nmr.pdb', ContentFile(
            b'MODEL        1\n' + PDB_CONTENT.replace(b'END\n', b'ENDMDL\n') +
            b'MODEL        2\n' + PDB_CONTENT.replace(b'END\n', b'ENDMDL\n') + b'END\n'))
        self.client.login(username='uploader', password='pw')
        response = self.client.get(reverse('pdb_grapher:detail', args=[pdb.pk]))
        self.assertContains(response, '<dd id="stat-models">2</dd>', html=True)
        self.assertContains(response, 'model-select-wrap')

    def test_personal_file_only_visible_to_uploader(self):
        pdb = self.make_pdb(personal=True)
        self.client.login(username='member', password='pw')
        self.assertEqual(self.client.get(reverse('pdb_grapher:detail', args=[pdb.pk])).status_code, 404)

    def test_only_leader_removes_from_group(self):
        pdb = self.make_pdb(groups=[self.group_a, self.group_b])
        url = reverse('pdb_grapher:remove_from_group', args=[pdb.pk, self.group_a.pk])
        self.client.login(username='member', password='pw')
        self.assertEqual(self.client.post(url).status_code, 403)
        self.client.login(username='leader', password='pw')
        self.client.post(url)
        self.assertEqual(list(pdb.groups.all()), [self.group_b])

    def test_removal_from_last_group_deletes_file(self):
        pdb = self.make_pdb(groups=[self.group_a])
        path = pdb.file.path
        self.client.login(username='leader', password='pw')
        self.client.post(reverse('pdb_grapher:remove_from_group', args=[pdb.pk, self.group_a.pk]))
        self.assertFalse(PDBFile.objects.exists())
        self.assertFalse(os.path.exists(path))

    def test_removal_keeps_personal_file(self):
        pdb = self.make_pdb(personal=True, groups=[self.group_a])
        self.client.login(username='leader', password='pw')
        self.client.post(reverse('pdb_grapher:remove_from_group', args=[pdb.pk, self.group_a.pk]))
        self.assertTrue(PDBFile.objects.filter(pk=pdb.pk).exists())

    def test_only_uploader_edits_and_deletes(self):
        pdb = self.make_pdb(groups=[self.group_a])
        self.client.login(username='leader', password='pw')
        self.assertEqual(self.client.get(reverse('pdb_grapher:edit', args=[pdb.pk])).status_code, 403)
        self.assertEqual(self.client.post(reverse('pdb_grapher:delete', args=[pdb.pk])).status_code, 403)

        self.client.login(username='uploader', password='pw')
        self.client.post(reverse('pdb_grapher:edit', args=[pdb.pk]), {'personal': 'on', 'groups': [self.group_b.pk]})
        pdb.refresh_from_db()
        self.assertTrue(pdb.personal)
        self.assertEqual(list(pdb.groups.all()), [self.group_b])

        self.client.post(reverse('pdb_grapher:delete', args=[pdb.pk]))
        self.assertFalse(PDBFile.objects.exists())

    def test_uploader_manages_file_after_leaving_group(self):
        pdb = self.make_pdb(groups=[self.group_a])
        self.group_a.members.remove(self.uploader)
        self.client.login(username='uploader', password='pw')
        self.assertEqual(self.client.get(reverse('pdb_grapher:detail', args=[pdb.pk])).status_code, 200)
        self.client.post(reverse('pdb_grapher:edit', args=[pdb.pk]), {'personal': 'on'})
        pdb.refresh_from_db()
        self.assertFalse(pdb.groups.exists())

    def test_deleting_group_deletes_orphaned_files(self):
        orphan = self.make_pdb(groups=[self.group_a])
        shared = self.make_pdb(groups=[self.group_a, self.group_b])
        self.group_a.delete()
        self.assertFalse(PDBFile.objects.filter(pk=orphan.pk).exists())
        self.assertTrue(PDBFile.objects.filter(pk=shared.pk).exists())

    def test_deleting_uploader_keeps_group_files(self):
        """Deleting a user removes their personal-only files; shared files lose uploader and personal flag."""
        personal = self.make_pdb(personal=True)
        shared = self.make_pdb(personal=True, groups=[self.group_a])
        self.uploader.delete()
        self.assertFalse(PDBFile.objects.filter(pk=personal.pk).exists())
        shared.refresh_from_db()
        self.assertIsNone(shared.uploaded_by)
        self.assertFalse(shared.personal)

    def test_index_lists_files_with_actions(self):
        """The group leader sees "remove from group", the uploader sees "delete"."""
        pdb = self.make_pdb(groups=[self.group_a])
        remove_url = reverse('pdb_grapher:remove_from_group', args=[pdb.pk, self.group_a.pk])
        self.client.login(username='leader', password='pw')
        response = self.client.get(reverse('pdb_grapher:index'))
        self.assertContains(response, remove_url)
        self.assertNotContains(response, reverse('pdb_grapher:delete', args=[pdb.pk]))
        self.client.login(username='uploader', password='pw')
        response = self.client.get(reverse('pdb_grapher:index'), {'group': self.group_a.pk})
        self.assertContains(response, reverse('pdb_grapher:delete', args=[pdb.pk]))
        self.assertNotContains(response, remove_url)


class NoteTests(BaseTestCase):
    """Notes on a PDB file: creation, visibility (personal or per group), editing and deletion."""

    def setUp(self):
        super().setUp()
        self.pdb = self.make_pdb(groups=[self.group_a, self.group_b])

    def make_note(self, author, groups=()):
        """Creates a note on ``self.pdb``.

        Args:
            author: User who wrote the note.
            groups: Work groups the note is shared with; empty means personal.

        Returns:
            The saved ``Note`` instance.
        """
        note = Note.objects.create(pdb=self.pdb, author=author, text=f'Notiz von {author}')
        note.groups.set(groups)
        return note

    def post_note(self, **data):
        """Posts the form for a new note as the currently logged-in user.

        Args:
            **data: Form fields without prefix (e.g. ``text``, ``visibility``, ``groups``); the form
                prefix ``note-new-`` is added here.

        Returns:
            The response of the note_create view.
        """
        return self.client.post(reverse('pdb_grapher:note_create', args=[self.pdb.pk]),
                                {f'note-new-{key}': value for key, value in data.items()})

    def notes_on_detail_page(self, username):
        """Logs in as ``username`` and returns the primary keys of the notes shown on the detail page.

        Args:
            username: Name of the user to log in as.

        Returns:
            Set of ``Note`` primary keys in the ``notes`` context variable.
        """
        self.client.login(username=username, password='pw')
        response = self.client.get(reverse('pdb_grapher:detail', args=[self.pdb.pk]))
        return {note.pk for note in response.context['notes']}

    def test_create_personal_and_group_notes(self):
        self.client.login(username='member', password='pw')
        self.assertEqual(self.post_note(text='privat', visibility='personal', groups=[self.group_a.pk]).status_code, 200)
        self.assertEqual(self.post_note(text='geteilt', visibility='groups', groups=[self.group_a.pk]).status_code, 200)
        personal, shared = Note.objects.get(text='privat'), Note.objects.get(text='geteilt')
        self.assertEqual(personal.author, self.member)
        self.assertFalse(personal.groups.exists())
        self.assertEqual(list(shared.groups.all()), [self.group_a])

    def test_create_requires_group_and_membership(self):
        # Group visibility without a group, a group the member does not belong to, and empty text
        self.client.login(username='member', password='pw')
        self.assertEqual(self.post_note(text='x', visibility='groups').status_code, 400)
        self.assertEqual(self.post_note(text='x', visibility='groups', groups=[self.group_b.pk]).status_code, 400)
        self.assertEqual(self.post_note(text='', visibility='personal').status_code, 400)
        self.assertFalse(Note.objects.exists())

    def test_visibility(self):
        personal = self.make_note(self.member)
        shared_a = self.make_note(self.member, groups=[self.group_a])
        shared_b = self.make_note(self.outsider, groups=[self.group_b])
        # Personal notes only for their author; group notes for all members of the group
        self.assertEqual(self.notes_on_detail_page('member'), {personal.pk, shared_a.pk})
        self.assertEqual(self.notes_on_detail_page('leader'), {shared_a.pk})
        self.assertEqual(self.notes_on_detail_page('uploader'), {shared_a.pk, shared_b.pk})

    def test_only_author_edits(self):
        note = self.make_note(self.member, groups=[self.group_a])
        url = reverse('pdb_grapher:note_edit', args=[self.pdb.pk, note.pk])
        data = {f'note-{note.pk}-text': 'neu', f'note-{note.pk}-visibility': 'personal'}
        self.client.login(username='leader', password='pw')
        self.assertEqual(self.client.post(url, data).status_code, 403)
        self.client.login(username='member', password='pw')
        self.assertEqual(self.client.post(url, data).status_code, 200)
        note.refresh_from_db()
        self.assertEqual(note.text, 'neu')
        self.assertIsNotNone(note.edited_at)
        self.assertFalse(note.groups.exists())

    def test_only_author_deletes(self):
        note = self.make_note(self.member, groups=[self.group_a])
        url = reverse('pdb_grapher:note_delete', args=[self.pdb.pk, note.pk])
        for username in ('uploader', 'leader'):
            self.client.login(username=username, password='pw')
            self.assertEqual(self.client.post(url).status_code, 403)
        self.client.login(username='member', password='pw')
        self.client.post(url)
        self.assertFalse(Note.objects.filter(pk=note.pk).exists())

    def test_only_leader_removes_note_from_own_group(self):
        note = self.make_note(self.member, groups=[self.group_a, self.group_b])
        url = reverse('pdb_grapher:note_remove_from_group', args=[self.pdb.pk, note.pk, self.group_a.pk])
        for username in ('member', 'outsider'):  # author without leadership, or leader of the other group
            self.client.login(username=username, password='pw')
            self.assertEqual(self.client.post(url).status_code, 403)
        self.client.login(username='leader', password='pw')
        self.assertEqual(self.client.post(url).status_code, 200)
        self.assertEqual(list(note.groups.all()), [self.group_b])

    def test_removal_from_last_group_keeps_note_for_author(self):
        """A note removed from its last group stays as a personal note; without an author it is deleted."""
        note = self.make_note(self.member, groups=[self.group_a])
        orphan = self.make_note(self.member, groups=[self.group_a])
        # author=None simulates a note whose author account has been deleted
        Note.objects.filter(pk=orphan.pk).update(author=None)
        self.client.login(username='leader', password='pw')
        for n in (note, orphan):
            self.client.post(reverse('pdb_grapher:note_remove_from_group', args=[self.pdb.pk, n.pk, self.group_a.pk]))
        self.assertFalse(Note.objects.filter(pk=orphan.pk).exists())
        self.assertEqual(self.notes_on_detail_page('leader'), set())
        self.assertEqual(self.notes_on_detail_page('member'), {note.pk})

    def test_deleting_author_keeps_group_notes(self):
        personal = self.make_note(self.member)
        shared = self.make_note(self.member, groups=[self.group_a])
        self.member.delete()
        self.assertFalse(Note.objects.filter(pk=personal.pk).exists())
        shared.refresh_from_db()
        self.assertIsNone(shared.author)

    def test_deleting_group_deletes_orphaned_notes(self):
        """Notes without author and without remaining group are deleted with the group."""
        note = self.make_note(self.member, groups=[self.group_a])
        Note.objects.filter(pk=note.pk).update(author=None)
        kept = self.make_note(self.member, groups=[self.group_a])
        self.group_a.delete()
        self.assertFalse(Note.objects.filter(pk=note.pk).exists())
        self.assertTrue(Note.objects.filter(pk=kept.pk).exists())


class IndexTests(BaseTestCase):
    """File list on the start page: search, date of the latest visible note, upload dialog."""

    def index_files(self, **params):
        """Fetches the start page as ``uploader`` and returns the listed files.

        Args:
            **params: Query parameters, e.g. ``q`` (search term) or ``group``.

        Returns:
            List of ``PDBFile`` objects from the ``pdb_files`` context variable.
        """
        self.client.login(username='uploader', password='pw')
        return list(self.client.get(reverse('pdb_grapher:index'), params).context['pdb_files'])

    def test_search_by_name(self):
        match = self.make_pdb(personal=True)
        other = self.make_pdb(personal=True)
        PDBFile.objects.filter(pk=match.pk).update(name='Lysozym.pdb')
        self.assertEqual([p.pk for p in self.index_files(q='lyso')], [match.pk])
        self.assertEqual({p.pk for p in self.index_files()}, {match.pk, other.pk})

    def test_latest_note_only_counts_visible_notes(self):
        pdb = self.make_pdb(groups=[self.group_a])
        visible = Note.objects.create(pdb=pdb, author=self.member, text='a')
        visible.groups.set([self.group_a])
        hidden = Note.objects.create(pdb=pdb, author=self.member, text='b')  # personal, only for its author
        # The hidden note is newer, but must not count for the uploader who cannot see it
        Note.objects.filter(pk=hidden.pk).update(created_at=visible.created_at.replace(year=2099))
        [row] = self.index_files()
        self.assertEqual(row.latest_note_at, visible.created_at)

    def test_failed_upload_reopens_dialog(self):
        # The upload fails because no availability (personal or group) is selected
        self.client.login(username='uploader', password='pw')
        response = self.client.post(reverse('pdb_grapher:index'), {'pdb_file': SimpleUploadedFile('x.pdb', PDB_CONTENT)})
        self.assertTrue(response.context['show_upload'])
        self.assertFalse(self.client.get(reverse('pdb_grapher:index')).context['show_upload'])


class RCSBImportTests(BaseTestCase):
    """Import of a structure by its PDB ID from the RCSB server; the HTTP request is mocked."""

    def post_upload(self, **data):
        """Posts the upload form as ``uploader`` with "personal" availability.

        Args:
            **data: Additional form fields, e.g. ``pdb_id`` or ``pdb_file``.

        Returns:
            The response of the index view.
        """
        self.client.login(username='uploader', password='pw')
        return self.client.post(reverse('pdb_grapher:index'), {'personal': 'on', **data})

    def mock_rcsb(self, content=PDB_CONTENT, error=None):
        """Patches ``urlopen`` in ``pdb_grapher.rcsb`` so that no network access takes place.

        Args:
            content: Bytes returned as the downloaded file.
            error: Exception raised instead of returning a response, e.g. an ``HTTPError``.

        Returns:
            A ``mock.patch`` context manager; entering it yields the ``urlopen`` mock.
        """
        response = mock.MagicMock()
        response.__enter__.return_value.read.return_value = content
        return mock.patch('pdb_grapher.rcsb.urllib.request.urlopen', return_value=response, side_effect=error)

    def test_import_by_id(self):
        with self.mock_rcsb() as urlopen:
            response = self.post_upload(pdb_id=' 1lyz ')
        self.assertEqual(urlopen.call_args.args[0].full_url, 'https://files.rcsb.org/download/1LYZ.pdb')
        pdb = PDBFile.objects.get()
        self.assertRedirects(response, reverse('pdb_grapher:detail', args=[pdb.pk]))
        self.assertEqual(pdb.name, '1LYZ.pdb')
        with pdb.file.open('rb') as f:
            self.assertEqual(f.read(), PDB_CONTENT)

    def test_invalid_id_is_not_requested(self):
        # Malformed IDs are rejected by form validation before any request is made
        with self.mock_rcsb() as urlopen:
            response = self.post_upload(pdb_id='abc')
        urlopen.assert_not_called()
        self.assertTrue(response.context['form'].errors['pdb_id'])
        self.assertFalse(PDBFile.objects.exists())

    def test_unknown_id(self):
        error = urllib.error.HTTPError('url', 404, 'Not Found', {}, None)
        with self.mock_rcsb(error=error):
            response = self.post_upload(pdb_id='9ZZZ')
        self.assertIn('keine Datei im PDB-Format', response.context['form'].errors['pdb_id'][0])
        self.assertFalse(PDBFile.objects.exists())

    def test_requires_exactly_one_source(self):
        # Neither a PDB ID nor a file, then both at once
        with self.mock_rcsb() as urlopen:
            self.post_upload()
            self.post_upload(pdb_id='1LYZ', pdb_file=SimpleUploadedFile('x.pdb', PDB_CONTENT))
        urlopen.assert_not_called()
        self.assertFalse(PDBFile.objects.exists())

    def test_no_request_without_availability(self):
        with self.mock_rcsb() as urlopen:
            self.client.login(username='uploader', password='pw')
            self.client.post(reverse('pdb_grapher:index'), {'pdb_id': '1LYZ'})
        urlopen.assert_not_called()


class LanguageTests(BaseTestCase):
    """Interface language: default German, switching, browser language, admin-managed language list.

    The available languages are ``Language`` rows managed in the admin; the English entry (en-gb)
    is created by the data migration ``0003_language`` and therefore exists in the test database.
    """

    def setUp(self):
        super().setUp()
        self.client.login(username='uploader', password='pw')

    def index(self, **extra):
        """Fetches the start page.

        Args:
            **extra: Extra request headers, e.g. ``HTTP_ACCEPT_LANGUAGE``.

        Returns:
            The response of the index view.
        """
        return self.client.get(reverse('pdb_grapher:index'), **extra)

    def test_default_is_german(self):
        response = self.index()
        self.assertContains(response, '<html lang="de"')
        self.assertContains(response, 'Abmelden')

    def test_switch_to_english(self):
        response = self.client.post(reverse('pdb_grapher:set_language'), {'language': 'en-gb', 'next': '/?group=personal'})
        self.assertRedirects(response, '/?group=personal', fetch_redirect_response=False)
        response = self.index()
        self.assertContains(response, '<html lang="en-gb"')
        self.assertContains(response, 'Sign out')
        self.assertContains(response, 'All available')

    def test_only_languages_from_admin(self):
        # French is not enabled in the admin: neither selectable nor taken from the browser
        response = self.client.post(reverse('pdb_grapher:set_language'), {'language': 'fr'})
        self.assertNotIn(settings.LANGUAGE_COOKIE_NAME, response.cookies)
        self.assertContains(self.index(HTTP_ACCEPT_LANGUAGE='fr'), '<html lang="de"')
        # A language cookie for a language that was removed in the admin is ignored
        self.client.cookies[settings.LANGUAGE_COOKIE_NAME] = 'en-gb'
        Language.objects.filter(code='en-gb').delete()
        self.assertContains(self.index(), '<html lang="de"')

    def test_browser_language_matches_variant(self):
        # The browser asks for en-US/en, the closest enabled variant is en-gb
        self.assertContains(self.index(HTTP_ACCEPT_LANGUAGE='en-US,en;q=0.8'), '<html lang="en-gb"')

    def test_unsafe_next_url_is_ignored(self):
        response = self.client.post(reverse('pdb_grapher:set_language'),
                                    {'language': 'en-gb', 'next': 'https://example.com/'})
        self.assertRedirects(response, reverse('pdb_grapher:index'), fetch_redirect_response=False)

    def test_selector_lists_admin_languages(self):
        self.assertContains(self.index(), 'name="language" value="en-gb" lang="en-gb"')
        Language.objects.filter(code='en-gb').delete()
        self.assertNotContains(self.index(), 'name="language"')

    def test_anonymous_selector_on_login_page(self):
        self.client.logout()
        response = self.client.get(reverse('login'))
        self.assertContains(response, '<option value="en-gb" lang="en-gb" >English</option>', html=False)
        self.assertContains(response, 'id="theme-toggle"')
        self.assertNotContains(response, 'account-menu')

    def test_javascript_catalog(self):
        # Translations used by the JavaScript code (here an amino acid name)
        self.client.cookies[settings.LANGUAGE_COOKIE_NAME] = 'en-gb'
        response = self.client.get(reverse('javascript-catalog'))
        self.assertContains(response, 'Aspartic acid')

    def test_server_messages_are_translated(self):
        self.client.cookies[settings.LANGUAGE_COOKIE_NAME] = 'en-gb'
        pdb = self.make_pdb(personal=True)
        response = self.client.post(reverse('pdb_grapher:delete', args=[pdb.pk]), follow=True)
        self.assertContains(response, '“test.pdb” has been deleted.')


class AccountTests(BaseTestCase):
    """Account page and menu: user data, password change, stored language and theme settings."""

    def setUp(self):
        super().setUp()
        self.client.login(username='uploader', password='pw')
        self.url = reverse('pdb_grapher:account')

    def test_requires_login(self):
        self.client.logout()
        response = self.client.get(self.url)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url}", fetch_redirect_response=False)

    def test_page_shows_user_and_groups(self):
        self.uploader.first_name, self.uploader.last_name = 'Ada', 'Lovelace'
        self.uploader.save()
        response = self.client.get(self.url)
        self.assertContains(response, 'Ada Lovelace')
        self.assertContains(response, '<span class="badge text-bg-primary">A</span>', html=True)
        self.assertContains(response, '<span class="badge text-bg-primary">B</span>', html=True)

    def test_menu_in_navbar(self):
        # Logged-in users choose the theme in the account menu instead of the separate toggle button
        response = self.client.get(reverse('pdb_grapher:index'))
        self.assertContains(response, 'account-menu')
        self.assertContains(response, f'href="{self.url}"')
        self.assertContains(response, 'name="theme" value="dark"')
        self.assertContains(response, 'Abmelden')
        self.assertNotContains(response, 'id="theme-toggle"')

    def test_change_password(self):
        response = self.client.post(self.url, {
            'action': 'password', 'old_password': 'pw',
            'new_password1': 'Bondamino-2026!', 'new_password2': 'Bondamino-2026!',
        })
        self.assertRedirects(response, self.url)
        self.uploader.refresh_from_db()
        self.assertTrue(self.uploader.check_password('Bondamino-2026!'))
        # The session stays valid after the password change
        self.assertEqual(self.client.get(self.url).status_code, 200)

    def test_change_password_wrong_old_password(self):
        response = self.client.post(self.url, {
            'action': 'password', 'old_password': 'falsch',
            'new_password1': 'Bondamino-2026!', 'new_password2': 'Bondamino-2026!',
        })
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'is-invalid')
        self.uploader.refresh_from_db()
        self.assertTrue(self.uploader.check_password('pw'))

    def test_save_settings(self):
        response = self.client.post(self.url, {'action': 'settings', 'language': 'en-gb', 'theme': 'dark'})
        self.assertRedirects(response, self.url, fetch_redirect_response=False)
        saved = UserSettings.objects.get(user=self.uploader)
        self.assertEqual((saved.language, saved.theme), ('en-gb', 'dark'))
        response = self.client.get(reverse('pdb_grapher:index'))
        self.assertContains(response, '<html lang="en-gb" data-bs-theme="dark" data-theme-source="user">')

    def test_settings_reject_unknown_language(self):
        response = self.client.post(self.url, {'action': 'settings', 'language': 'fr', 'theme': 'dark'})
        self.assertEqual(response.status_code, 200)
        self.assertFalse(UserSettings.objects.filter(user=self.uploader).exists())

    def test_stored_language_wins_over_cookie(self):
        UserSettings.objects.create(user=self.uploader, language='en-gb')
        self.client.cookies[settings.LANGUAGE_COOKIE_NAME] = 'de'
        self.assertContains(self.client.get(reverse('pdb_grapher:index')), '<html lang="en-gb"')

    def test_language_menu_saves_setting(self):
        self.client.post(reverse('pdb_grapher:set_language'), {'language': 'en-gb', 'next': '/'})
        self.assertEqual(UserSettings.objects.get(user=self.uploader).language, 'en-gb')

    def test_set_theme(self):
        """Fetch requests get JSON (or 400 for unknown themes), form posts are redirected to ``next``."""
        url = reverse('pdb_grapher:set_theme')
        response = self.client.post(url, {'theme': 'dark'}, HTTP_X_REQUESTED_WITH='fetch')
        self.assertEqual(response.json(), {'theme': 'dark'})
        self.assertEqual(UserSettings.objects.get(user=self.uploader).theme, 'dark')
        response = self.client.post(url, {'theme': 'pink'}, HTTP_X_REQUESTED_WITH='fetch')
        self.assertEqual(response.status_code, 400)
        response = self.client.post(url, {'theme': 'light', 'next': '/?group=personal'})
        self.assertRedirects(response, '/?group=personal', fetch_redirect_response=False)
        self.assertEqual(UserSettings.objects.get(user=self.uploader).theme, 'light')

    def test_set_theme_requires_login(self):
        self.client.logout()
        self.client.post(reverse('pdb_grapher:set_theme'), {'theme': 'dark'})
        self.assertFalse(UserSettings.objects.exists())


def mail_link(message):
    """Extracts the path of the first link in an email.

    Args:
        message: An email from ``django.core.mail.outbox``.

    Returns:
        The URL path (without scheme and host), ready to be requested with the test client.
    """
    return re.search(r'https?://[^/\s]+(/\S+)', message.body).group(1)


class InvitationTests(BaseTestCase):
    """Users created in the admin without a password get an invitation email with a one-time link."""

    def setUp(self):
        super().setUp()
        self.admin_user = User.objects.create_superuser('admin', 'admin@example.org', 'pw')
        self.client.login(username='admin', password='pw')

    def add_user(self, **data):
        """Submits the admin "add user" form for the user ``neu``, who becomes a member of group A.

        The ``WorkGroup_members-*`` fields belong to the inline formset for work group membership.
        The invitation email is sent in an on-commit callback, which is executed here explicitly
        because a ``TestCase`` never commits its transaction.

        Args:
            **data: Form fields overriding the defaults, e.g. ``email=''``.

        Returns:
            The response of the admin view.
        """
        fields = {'username': 'neu', 'email': 'neu@example.org', 'first_name': 'Neue', 'last_name': 'Person',
                  'WorkGroup_members-TOTAL_FORMS': '1', 'WorkGroup_members-INITIAL_FORMS': '0',
                  'WorkGroup_members-0-workgroup': self.group_a.pk, **data}
        with self.captureOnCommitCallbacks(execute=True):
            return self.client.post(reverse('admin:auth_user_add'), fields)

    def test_admin_creates_user_without_password_and_sends_invitation(self):
        response = self.add_user()
        self.assertEqual(response.status_code, 302)
        user = User.objects.get(username='neu')
        self.assertFalse(user.has_usable_password())
        self.assertEqual(len(mail.outbox), 1)
        message = mail.outbox[0]
        self.assertEqual(message.to, ['neu@example.org'])
        self.assertIn('Neue Person', message.body)
        self.assertIn('24 Stunden', message.body)
        self.assertTrue(mail_link(message).startswith('/invitation/'))
        self.assertIn(user, self.group_a.members.all())

    def test_email_is_required(self):
        self.add_user(email='')
        self.assertFalse(User.objects.filter(username='neu').exists())
        self.assertEqual(mail.outbox, [])

    def test_no_login_before_password_is_set(self):
        self.add_user()
        self.client.logout()
        self.assertFalse(self.client.login(username='neu', password=''))

    def test_link_sets_password_and_logs_in(self):
        self.add_user()
        self.client.logout()
        response = self.client.get(mail_link(mail.outbox[0]), follow=True)
        self.assertContains(response, 'Passwort festlegen')
        set_url = response.redirect_chain[-1][0]
        response = self.client.post(set_url, {'new_password1': 'Bondamino-2026!', 'new_password2': 'Bondamino-2026!'})
        self.assertRedirects(response, reverse('pdb_grapher:index'))
        self.assertTrue(User.objects.get(username='neu').check_password('Bondamino-2026!'))
        self.assertEqual(self.client.get(reverse('pdb_grapher:account')).status_code, 200)
        # The link can only be used once
        self.client.logout()
        self.assertContains(self.client.get(mail_link(mail.outbox[0]), follow=True), 'ungültig oder abgelaufen')

    def test_link_expires_after_24_hours(self):
        self.add_user()
        self.client.logout()
        later = datetime.datetime.now() + datetime.timedelta(hours=24, minutes=1)
        # The token generator checks the token age against _now(); pretend the link is opened 24 h 1 min later
        with mock.patch('django.contrib.auth.tokens.PasswordResetTokenGenerator._now', return_value=later):
            response = self.client.get(mail_link(mail.outbox[0]), follow=True)
        self.assertContains(response, 'ungültig oder abgelaufen')

    def test_resend_invitation_only_without_password(self):
        # Admin action on the user list: "neu" has no password yet, "member" has one and gets no mail
        self.add_user()
        mail.outbox.clear()
        ids = [User.objects.get(username='neu').pk, self.member.pk]
        self.client.post(reverse('admin:auth_user_changelist'),
                         {'action': 'resend_invitation', '_selected_action': ids})
        self.assertEqual([m.to for m in mail.outbox], [['neu@example.org']])


class PasswordResetTests(BaseTestCase):
    """The "forgot password" flow: request a link by email address, then set a new password."""

    def setUp(self):
        super().setUp()
        self.member.email = 'member@example.org'
        self.member.save()

    def request_link(self, email):
        """Submits the password reset form.

        Args:
            email: Email address entered in the form.

        Returns:
            The response of the password_reset view.
        """
        return self.client.post(reverse('pdb_grapher:password_reset'), {'email': email})

    def test_login_page_links_to_reset(self):
        self.assertContains(self.client.get(reverse('login')), reverse('pdb_grapher:password_reset'))

    def test_reset_password(self):
        # The email address is matched case-insensitively; the mail contains the username
        response = self.request_link('MEMBER@example.org')
        self.assertRedirects(response, reverse('pdb_grapher:password_reset_done'))
        self.assertEqual(len(mail.outbox), 1)
        self.assertIn('member', mail.outbox[0].body)
        response = self.client.get(mail_link(mail.outbox[0]), follow=True)
        self.assertContains(response, 'Neues Passwort')
        response = self.client.post(response.redirect_chain[-1][0],
                                    {'new_password1': 'Bondamino-2026!', 'new_password2': 'Bondamino-2026!'})
        self.assertRedirects(response, reverse('pdb_grapher:index'))
        self.member.refresh_from_db()
        self.assertTrue(self.member.check_password('Bondamino-2026!'))

    def test_unknown_email_sends_nothing_but_same_page(self):
        response = self.request_link('niemand@example.org')
        self.assertRedirects(response, reverse('pdb_grapher:password_reset_done'))
        self.assertEqual(mail.outbox, [])

    def test_works_for_user_without_password(self):
        # E.g. an invited user whose invitation link has expired
        self.member.set_unusable_password()
        self.member.save()
        self.request_link('member@example.org')
        self.assertEqual(len(mail.outbox), 1)

    def test_inactive_user_gets_no_mail(self):
        self.member.is_active = False
        self.member.save()
        self.request_link('member@example.org')
        self.assertEqual(mail.outbox, [])

    def test_mismatched_passwords(self):
        self.request_link('member@example.org')
        response = self.client.get(mail_link(mail.outbox[0]), follow=True)
        response = self.client.post(response.redirect_chain[-1][0],
                                    {'new_password1': 'Bondamino-2026!', 'new_password2': 'anders'})
        self.assertContains(response, 'is-invalid')
        self.member.refresh_from_db()
        self.assertTrue(self.member.check_password('pw'))


class ThrottleTests(BaseTestCase):
    """Limits on failed logins and password reset requests (see the ``throttle`` module).

    The counters live in the database cache (CACHES in settings.py), i.e. in a table of the test
    database, so every test starts with fresh counters when its transaction is rolled back.
    """

    def login(self, username='member', password='falsch', ip='10.0.0.1', url=None, **extra):
        """Posts the login form; by default with a wrong password ("falsch").

        Args:
            username: Username entered in the form.
            password: Password entered in the form.
            ip: Client address (``REMOTE_ADDR``).
            url: Login URL; defaults to the normal login page (the admin login is tested as well).
            **extra: Further request options, e.g. ``HTTP_X_FORWARDED_FOR`` or ``secure``.

        Returns:
            The response of the login view.
        """
        return self.client.post(url or reverse('login'), {'username': username, 'password': password},
                                REMOTE_ADDR=ip, **extra)

    def test_login_locked_after_failures(self):
        for _ in range(throttle.LOGIN_FAILURES_PER_USER):
            self.assertContains(self.login(), 'Bitte Benutzername und Passwort eingeben')
        # Afterwards locked even with the correct password, but not from another address
        # (the per-user counter is keyed by IP address and username)
        self.assertContains(self.login(password='pw'), 'Zu viele fehlgeschlagene Anmeldeversuche')
        self.assertNotIn('_auth_user_id', self.client.session)
        self.assertRedirects(self.login(password='pw', ip='10.0.0.2'), reverse('pdb_grapher:index'))

    def test_success_resets_counter(self):
        # One failure short of the lock, then a success; one more failure afterwards must not lock
        for _ in range(throttle.LOGIN_FAILURES_PER_USER - 1):
            self.login()
        self.assertRedirects(self.login(password='pw'), reverse('pdb_grapher:index'))
        self.client.logout()
        self.login()
        self.assertRedirects(self.login(password='pw'), reverse('pdb_grapher:index'))

    def test_lock_expires(self):
        for _ in range(throttle.LOGIN_FAILURES_PER_USER):
            self.login()
        cache.clear()  # as if the lock period had expired
        self.assertRedirects(self.login(password='pw'), reverse('pdb_grapher:index'))

    def test_failures_per_ip_across_usernames(self):
        # Many usernames from one address lock that address for all users, including "member"
        for i in range(throttle.LOGIN_FAILURES_PER_IP):
            self.login(username=f'user{i}')
        self.assertContains(self.login(password='pw'), 'Zu viele fehlgeschlagene Anmeldeversuche')

    def test_admin_login_throttled(self):
        self.member.is_staff = True
        self.member.save()
        url = reverse('admin:login')
        for _ in range(throttle.LOGIN_FAILURES_PER_USER):
            self.login(url=url)
        self.assertContains(self.login(password='pw', url=url), 'Zu viele fehlgeschlagene Anmeldeversuche')

    def test_forwarded_for_only_behind_proxy(self):
        for i in range(throttle.LOGIN_FAILURES_PER_USER):
            # Without HTTPS=1, X-Forwarded-For is ignored (otherwise anyone could fake the address)
            self.login(HTTP_X_FORWARDED_FOR=f'192.0.2.{i}')
        self.assertContains(self.login(password='pw'), 'Zu viele fehlgeschlagene Anmeldeversuche')
        with override_settings(HTTPS=True):
            # Behind the proxy the last entry counts, which the proxy itself appends; 192.0.2.7 is not
            # locked although REMOTE_ADDR (the proxy) and the client-supplied first entry are 10.0.0.1
            self.assertRedirects(self.login(password='pw', HTTP_X_FORWARDED_FOR='10.0.0.1, 192.0.2.7', secure=True),
                                 reverse('pdb_grapher:index'))

    def test_missing_cache_table_does_not_block_login(self):
        # Simulates a missing cache table (createcachetable not run): no throttling, but an error is logged
        with mock.patch('pdb_grapher.throttle.cache.get', side_effect=DatabaseError), \
                mock.patch('pdb_grapher.throttle.cache.add', side_effect=DatabaseError), \
                self.assertLogs('pdb_grapher.throttle', 'ERROR'):
            self.assertRedirects(self.login(password='pw'), reverse('pdb_grapher:index'))

    def test_password_reset_limited_per_email(self):
        self.member.email = 'member@example.org'
        self.member.save()
        with self.assertLogs('pdb_grapher.forms', 'WARNING'):
            for _ in range(throttle.RESET_REQUESTS_PER_EMAIL + 2):
                response = self.client.post(reverse('pdb_grapher:password_reset'), {'email': 'member@example.org'})
                # Always the same response, so nothing is revealed about the account
                self.assertRedirects(response, reverse('pdb_grapher:password_reset_done'))
        self.assertEqual(len(mail.outbox), throttle.RESET_REQUESTS_PER_EMAIL)

    def test_password_reset_limited_per_ip(self):
        self.member.email = 'member@example.org'
        self.member.save()
        # Use up the per-IP limit with unknown addresses, then a known address gets no mail either
        for i in range(throttle.RESET_REQUESTS_PER_IP):
            self.client.post(reverse('pdb_grapher:password_reset'), {'email': f'x{i}@example.org'})
        with self.assertLogs('pdb_grapher.forms', 'WARNING'):
            self.client.post(reverse('pdb_grapher:password_reset'), {'email': 'member@example.org'})
        self.assertEqual(mail.outbox, [])


class OutputSafetyTests(BaseTestCase):
    """Uploaded content cannot inject HTML, and internal error details are not shown to users."""

    def test_graph_data_cannot_break_out_of_script(self):
        # The atom name "<!--" (four characters, columns 13-16) would confuse the parsing of a <script> block;
        # json_script must escape "<" as <
        content = (atom_line(1, 'N', 1, 0.0) + atom_line(2, 'CA', 1, 1.46)
                   + f"ATOM  {3:5d} <!-- ALA A   1    {2.0:8.3f}{0:8.3f}{0:8.3f}  1.00  0.00           C\n"
                   + atom_line(4, 'C', 1, 2.98)).encode()
        pdb = PDBFile(name='x.pdb', uploaded_by=self.uploader, personal=True)
        pdb.file.save('x.pdb', ContentFile(content), save=False)
        pdb.save()
        self.client.force_login(self.uploader)
        html = self.client.get(reverse('pdb_grapher:detail', args=[pdb.pk])).content.decode()
        script = re.search(r'<script id="graph-data" type="application/json">(.*?)</script>', html, re.S).group(1)
        self.assertNotIn('<', script)
        self.assertIn('\\u003C!--', script)

    def test_processing_error_details_not_shown(self):
        # The exception message contains a server path that must only be logged, not shown
        pdb = self.make_pdb(personal=True)
        self.client.force_login(self.uploader)
        with mock.patch('pdb_grapher.views.parse_pdb', side_effect=OSError('/srv/geheim/pfad')), \
                self.assertLogs('pdb_grapher.views', 'ERROR'):
            detail = self.client.get(reverse('pdb_grapher:detail', args=[pdb.pk]), follow=True)
            data = self.client.get(reverse('pdb_grapher:model_data', args=[pdb.pk]))
        self.assertContains(detail, 'konnte nicht verarbeitet werden')
        self.assertNotContains(detail, '/srv/geheim')
        self.assertEqual(data.status_code, 400)
        self.assertNotIn('/srv/geheim', data.json()['error'])


class NoteGroupChoiceTests(BaseTestCase):
    """Notes can only be shared with groups in which the file is available."""

    def setUp(self):
        super().setUp()
        # The uploader is a member of A and B; the file is shared with A only
        self.pdb = self.make_pdb(personal=True, groups=[self.group_a])
        self.client.force_login(self.uploader)

    def post_note(self, groups):
        """Posts a new note shared with ``groups`` and returns the response."""
        return self.client.post(reverse('pdb_grapher:note_create', args=[self.pdb.pk]), {
            'note-new-text': 'Notiz', 'note-new-visibility': 'groups', 'note-new-groups': [g.pk for g in groups],
        })

    def test_only_groups_with_access_are_offered(self):
        response = self.client.get(reverse('pdb_grapher:detail', args=[self.pdb.pk]))
        self.assertEqual(list(response.context['note_form'].fields['groups'].queryset), [self.group_a])

    def test_group_without_access_is_rejected(self):
        self.assertEqual(self.post_note([self.group_b]).status_code, 400)
        self.assertEqual(self.post_note([self.group_a]).status_code, 200)

    def test_note_keeps_its_groups_when_the_file_leaves_a_group(self):
        note = Note.objects.create(pdb=self.pdb, author=self.uploader, text='Notiz')
        note.groups.set([self.group_a])
        self.pdb.groups.remove(self.group_a)
        self.assertEqual(list(note.groups.all()), [self.group_a])
        # The author can still edit it without losing the group
        form = NoteForm(instance=note, user=self.uploader, pdb=self.pdb)
        self.assertIn(self.group_a, form.fields['groups'].queryset)


class AdminAvailabilityTests(BaseTestCase):
    """The admin form of a PDB file enforces the same availability rule as the file list."""

    def form(self, pdb, **data):
        """Returns the bound admin form for ``pdb`` with its current values changed by ``data``."""
        values = {'name': pdb.name, 'uploaded_by': pdb.uploaded_by_id, 'personal': pdb.personal,
                  'groups': [g.pk for g in pdb.groups.all()], **data}
        return PDBFileAdminForm(values, instance=pdb)

    def test_file_must_stay_available(self):
        pdb = self.make_pdb(personal=True, groups=[self.group_a])
        self.assertTrue(self.form(pdb).is_valid())
        self.assertFalse(self.form(pdb, personal=False, groups=[]).is_valid())
        self.assertTrue(self.form(pdb, personal=False).is_valid())
        # "Personal" does not count without an uploader
        self.assertFalse(self.form(pdb, uploaded_by='', groups=[]).is_valid())


class TimeZoneTests(BaseTestCase):
    """Dates and times are shown in the user's time zone; the server default is UTC."""

    def setUp(self):
        super().setUp()
        pdb = self.make_pdb(personal=True)
        PDBFile.objects.filter(pk=pdb.pk).update(
            uploaded_at=datetime.datetime(2026, 1, 15, 12, 0, tzinfo=datetime.timezone.utc))
        self.client.force_login(self.uploader)

    def shown(self):
        """Returns the index page as text."""
        return self.client.get(reverse('pdb_grapher:index')).content.decode()

    def test_default_and_user_time_zone(self):
        self.assertIn('15.01.2026 12:00', self.shown())
        UserSettings.objects.create(user=self.uploader, time_zone='Europe/Berlin')
        self.assertIn('15.01.2026 13:00', self.shown())

    def test_unknown_time_zone_falls_back_to_default(self):
        UserSettings.objects.create(user=self.uploader, time_zone='Mars/Olympus_Mons')
        self.assertIn('15.01.2026 12:00', self.shown())

    def test_save_time_zone_in_settings(self):
        url = reverse('pdb_grapher:account')
        response = self.client.post(url, {'action': 'settings', 'language': 'de', 'theme': 'light',
                                          'time_zone': 'Asia/Tokyo'})
        self.assertRedirects(response, url, fetch_redirect_response=False)
        self.assertEqual(UserSettings.objects.get(user=self.uploader).time_zone, 'Asia/Tokyo')
        response = self.client.post(url, {'action': 'settings', 'language': 'de', 'theme': 'light',
                                          'time_zone': 'Mars/Olympus_Mons'})
        self.assertEqual(response.status_code, 200)

    def test_settings_without_languages(self):
        # Without languages in the admin the language field is left out and the form can still be saved
        Language.objects.all().delete()
        url = reverse('pdb_grapher:account')
        self.assertNotContains(self.client.get(url), 'name="language"')
        response = self.client.post(url, {'action': 'settings', 'theme': 'dark', 'time_zone': ''})
        self.assertRedirects(response, url, fetch_redirect_response=False)
        self.assertEqual(UserSettings.objects.get(user=self.uploader).theme, 'dark')


class ModelDataTests(BaseTestCase):
    """The model parameter of model_data must name a model of the file."""

    def test_unknown_model_is_rejected(self):
        pdb = PDBFile(name='nmr.pdb', uploaded_by=self.uploader, personal=True)
        pdb.file.save('nmr.pdb', ContentFile(
            b'MODEL        1\n' + PDB_CONTENT.replace(b'END\n', b'ENDMDL\n') +
            b'MODEL        2\n' + PDB_CONTENT.replace(b'END\n', b'ENDMDL\n') + b'END\n'))
        self.client.force_login(self.uploader)
        url = reverse('pdb_grapher:model_data', args=[pdb.pk])
        self.assertEqual(self.client.get(url, {'model': '2'}).status_code, 200)
        response = self.client.get(url, {'model': '3'})
        self.assertEqual(response.status_code, 400)
        self.assertIn('error', response.json())


class ConfigTests(TestCase):
    """Configuration helpers in ``bondamino/config.py`` used by settings.py."""

    def test_sqlite_urls(self):
        # Three slashes: path relative to the base directory; four slashes: absolute path
        self.assertEqual(Path(parse_database_url('sqlite:///db.sqlite3', '/srv/app')['NAME']),
                         Path('/srv/app/db.sqlite3'))
        config = parse_database_url('sqlite:////data/db.sqlite3', '/srv/app')
        self.assertEqual(config['ENGINE'], 'django.db.backends.sqlite3')
        self.assertEqual(Path(config['NAME']), Path('/data/db.sqlite3'))

    def test_postgres_url_with_encoded_password(self):
        # Percent-encoded password ("p@ss:w/rd"), explicit port and query parameters as OPTIONS
        config = parse_database_url('postgres://bonda:p%40ss%3Aw%2Frd@db.example.org:5433/bondamino?sslmode=require', '/')
        self.assertEqual(config, {
            'ENGINE': 'django.db.backends.postgresql', 'NAME': 'bondamino', 'USER': 'bonda', 'PASSWORD': 'p@ss:w/rd',
            'HOST': 'db.example.org', 'PORT': '5433', 'OPTIONS': {'sslmode': 'require'},
        })

    def test_mariadb_url(self):
        config = parse_database_url('mariadb://u:pw@db/bondamino', '/')
        self.assertEqual((config['ENGINE'], config['PORT'], config['OPTIONS']),
                         ('django.db.backends.mysql', '', {'charset': 'utf8mb4'}))

    def test_unsupported_database(self):
        with self.assertRaises(ValueError):
            parse_database_url('mongodb://localhost/bondamino', '/')

    def test_env_file(self):
        # Comments and blank lines are skipped, quotes removed, "#" and "=" inside quoted values kept;
        # a missing file yields an empty dict
        path = os.path.join(tempfile.mkdtemp(), '.env')
        with open(path, 'w', encoding='utf-8') as f:
            f.write('# Kommentar\nDEBUG=1\nSECRET_KEY="a#b=c"\n\nALLOWED_HOSTS = a.org, b.org\n')
        self.assertEqual(read_env_file(path), {'DEBUG': '1', 'SECRET_KEY': 'a#b=c', 'ALLOWED_HOSTS': 'a.org, b.org'})
        self.assertEqual(read_env_file(path + '.missing'), {})

    def test_empty_values_count_as_unset(self):
        with mock.patch.dict(os.environ, {'X_EMPTY': '', 'X_ON': 'yes', 'X_LIST': 'a.org, b.org,'}):
            self.assertEqual(env('X_EMPTY', 'default'), 'default')
            self.assertTrue(env_bool('X_ON'))
            self.assertFalse(env_bool('X_EMPTY'))
            self.assertEqual(env_list('X_LIST'), ['a.org', 'b.org'])
            self.assertEqual(env_list('X_EMPTY', ['x']), ['x'])

    def test_no_secret_key_in_code(self):
        settings_file = os.path.join(settings.BASE_DIR, 'bondamino', 'settings.py')
        with open(settings_file, encoding='utf-8') as f:
            self.assertNotIn('django-insecure', f.read())


class SetupCommandTests(TestCase):
    """The interactive ``manage.py setup`` command that writes the .env configuration file."""

    def run_setup(self, answers, env_file):
        """Runs the setup command with scripted answers instead of keyboard input.

        ``input()`` and ``getpass.getpass()`` both take the next answer from the list, so the list must
        follow the exact order of the prompts in ``management/commands/setup.py``. An empty string
        accepts the default value. ``skip_finish`` skips the final step (migrate, createcachetable,
        collectstatic, createsuperuser).

        Args:
            answers: Answers to the prompts, in order.
            env_file: Path of the .env file to write.
        """
        answers = iter(answers)
        with mock.patch('builtins.input', lambda prompt='': next(answers)), \
                mock.patch('getpass.getpass', lambda prompt='': next(answers)):
            call_command('setup', env_file=env_file, skip_finish=True, stdout=open(os.devnull, 'w'))

    def test_writes_env_file(self):
        env_file = os.path.join(tempfile.mkdtemp(), '.env')
        self.run_setup([
            '1', 'db.sqlite3',                       # SQLite, database file
            'bondamino.example.org', 'y', 'n', '',   # domain, HTTPS, no DEBUG, time zone
            # SMTP, server, STARTTLS, default port (587), user, password, sender, no test email
            '1', 'smtp.example.org', '1', '', 'mailer', 's3cret', 'Bondamino <noreply@example.org>', '',
            '',                                      # write the file
        ], env_file)
        values = read_env_file(env_file)
        self.assertEqual(values['DATABASE_URL'], 'sqlite:///db.sqlite3')
        self.assertEqual(values['ALLOWED_HOSTS'], 'bondamino.example.org')
        self.assertEqual((values['HTTPS'], values['DEBUG']), ('1', '0'))
        self.assertEqual((values['EMAIL_HOST'], values['EMAIL_PORT'], values['EMAIL_USE_TLS']),
                         ('smtp.example.org', '587', '1'))
        self.assertEqual(values['EMAIL_HOST_PASSWORD'], 's3cret')
        self.assertGreaterEqual(len(values['SECRET_KEY']), 50)

    def test_keeps_secret_key_and_extra_values(self):
        env_file = os.path.join(tempfile.mkdtemp(), '.env')
        with open(env_file, 'w', encoding='utf-8') as f:
            f.write('SECRET_KEY=bestehend\nLOG_LEVEL=INFO\nDEBUG=1\n')
        # Accept all database and web defaults (DEBUG defaults to the existing value), choose
        # "print emails to the console" ('2') and write; the old file is kept as .env.bak
        self.run_setup(['', '', '', '', '', '', '2', ''], env_file)
        values = read_env_file(env_file)
        self.assertEqual((values['SECRET_KEY'], values['LOG_LEVEL'], values['DEBUG']), ('bestehend', 'INFO', '1'))
        self.assertTrue(os.path.exists(env_file + '.bak'))

    def test_localhost_and_time_zone(self):
        """localhost also allows 127.0.0.1; an unknown time zone is asked again."""
        env_file = os.path.join(tempfile.mkdtemp(), '.env')
        self.run_setup(['1', 'db.sqlite3', 'localhost', 'n', 'n', 'Mars/Olympus_Mons', 'Europe/Berlin', '2', ''],
                       env_file)
        values = read_env_file(env_file)
        self.assertEqual((values['ALLOWED_HOSTS'], values['TIME_ZONE']), ('localhost,127.0.0.1', 'Europe/Berlin'))


def atom_line(serial, name, res_seq, x):
    """Formats an ATOM line of an alanine in chain A at (x, 0, 0).

    Args:
        serial: Atom serial number (columns 7-11).
        name: Atom name of up to three characters, e.g. 'CA' (written from column 14 on, as usual for
            names of light elements); its first letter is used as the element symbol.
        res_seq: Residue number (columns 23-26).
        x: x coordinate in Å.

    Returns:
        The PDB line including the trailing newline.
    """
    return (f"ATOM  {serial:5d}  {name:<3} ALA A{res_seq:4d}    "
            f"{x:8.3f}{0:8.3f}{0:8.3f}  1.00  0.00           {name[0]}\n")


def residue_lines(res_seq, start, serial=1):
    """Formats the backbone atoms N, CA, C of an alanine on a straight line along the x axis.

    N is at ``start``, CA at ``start + 1.46`` and C at ``start + 2.98``. A following residue built
    with ``start + 4.31`` therefore has its N 1.33 Å behind this C, i.e. at peptide bond distance.

    Args:
        res_seq: Residue number.
        start: x coordinate of N in Å.
        serial: Atom serial number of N; CA and C get the next two numbers.

    Returns:
        List of three PDB lines.
    """
    return [atom_line(serial + i, name, res_seq, start + dx)
            for i, (name, dx) in enumerate((('N', 0.0), ('CA', 1.46), ('C', 2.98)))]


class ChainGapTests(TestCase):
    """Peptide bonds are only drawn between consecutive residues at bonding distance."""

    def write_pdb(self, lines):
        """Writes PDB lines to a temporary file that is removed after the test.

        Args:
            lines: PDB lines including newlines.

        Returns:
            Path of the temporary file.
        """
        fd, path = tempfile.mkstemp(suffix='.pdb')
        with os.fdopen(fd, 'w') as f:
            f.writelines(lines)
        self.addCleanup(os.remove, path)
        return path

    def peptide_bonds(self, graph):
        """Returns the peptide bonds of a graph as {(source node id, target node id): length}."""
        return {(e['source'], e['target']): e['length'] for e in graph['edges'] if e['type'] == 'peptide'}

    def test_gap_is_no_peptide_bond(self):
        # Residue 2 directly follows 1; residues 3 and 4 are missing before residue 5 (C-N distance ≈ 10.7 Å)
        path = self.write_pdb(residue_lines(1, 0.0) + residue_lines(2, 4.31, 4) + residue_lines(5, 18.0, 7))
        bonds = self.peptide_bonds(build_graph(parse_pdb(path)))
        self.assertEqual(bonds, {('A_1_C', 'A_2_N'): 1.33})

    def test_gap_in_one_model_is_ignored(self):
        # Model 1: residues 1 and 2 are bonded; model 2: residue 2 is far away (C-N 11 Å).
        # The averaged bond length only uses model 1; the ensemble entry is [mean, median, number of models]
        # and the longest bond of all models is CA-C (1.52 Å), not the 11 Å gap.
        lines = (['MODEL        1\n'] + residue_lines(1, 0.0) + residue_lines(2, 4.31, 4) + ['ENDMDL\n'] +
                 ['MODEL        2\n'] + residue_lines(1, 0.0) + residue_lines(2, 14.0, 4) + ['ENDMDL\n'])
        path = self.write_pdb(lines)
        self.assertEqual(self.peptide_bonds(build_average_graph(path, [1, 2])), {('A_1_C', 'A_2_N'): 1.33})
        ensemble = build_bond_ensemble(path, [1, 2])
        self.assertEqual(ensemble['bonds']['A_1_C|A_2_N'], [1.33, 1.33, 1])
        self.assertEqual(ensemble['max']['length'], 1.52)

    def test_models_without_number_are_numbered_by_position(self):
        # MODEL records without a readable serial number (columns 11-14 blank): list_models() and
        # parse_pdb() must agree on the numbering, otherwise selecting model 2 would show no atoms.
        lines = (['MODEL\n'] + residue_lines(1, 0.0) + ['ENDMDL\n'] +
                 ['MODEL\n'] + residue_lines(1, 5.0) + ['ENDMDL\n'])
        path = self.write_pdb(lines)
        self.assertEqual(list_models(path), [1, 2])
        residue = next(iter(parse_pdb(path, model_id=2).values()))
        self.assertEqual(residue['N']['x'], 5.0)


def pdb_line(serial, name, res_name, res_seq, x, alt=' ', i_code=' ', occupancy=1.0, record='ATOM  '):
    """Formats an ATOM/HETATM line in chain A at (x, 0, 0) with all fields relevant for conformations.

    Args:
        serial: Atom serial number (columns 7-11).
        name: Atom name of up to three characters (from column 14 on); its first letter is the element.
        res_name: Residue name, e.g. 'PRO' (columns 18-20).
        res_seq: Residue number (columns 23-26).
        x: x coordinate in Å.
        alt: Alternate location indicator (column 17), ' ' for none.
        i_code: Insertion code (column 27), ' ' for none.
        occupancy: Occupancy (columns 55-60).
        record: Record name padded to six characters, 'ATOM  ' or 'HETATM'.

    Returns:
        The PDB line including the trailing newline.
    """
    return (f"{record}{serial:5d}  {name:<3}{alt}{res_name:>3} A{res_seq:4d}{i_code}   "
            f"{x:8.3f}{0:8.3f}{0:8.3f}{occupancy:6.2f}  0.00           {name[0]}\n")


# Residue 1 with CA in two conformations: A (proline, occupancy 0.3) and B (serine, occupancy 0.7).
# The different residue names model microheterogeneity (two amino acids at the same position).
ALTLOC_LINES = [
    pdb_line(1, 'N', 'PRO', 1, 0.0),
    pdb_line(2, 'CA', 'PRO', 1, 1.46, alt='A', occupancy=0.3),
    pdb_line(3, 'CA', 'SER', 1, 1.50, alt='B', occupancy=0.7),
    pdb_line(4, 'C', 'PRO', 1, 2.98),
]


class ConformationParsingTests(TestCase):
    """Parsing of insertion codes and alternate locations (altlocs) in ``parse_pdb()``."""

    def write_pdb(self, lines):
        """Writes PDB lines to a temporary file that is removed after the test.

        Args:
            lines: PDB lines including newlines.

        Returns:
            Path of the temporary file.
        """
        fd, path = tempfile.mkstemp(suffix='.pdb')
        with os.fdopen(fd, 'w') as f:
            f.writelines(lines)
        self.addCleanup(os.remove, path)
        return path

    def test_insertion_codes_are_separate_residues(self):
        # Residues 52, 52A and 53, each with N, CA, C and spaced for peptide bonds between them
        lines = []
        for i, (seq, i_code) in enumerate(((52, ' '), (52, 'A'), (53, ' '))):
            lines += [pdb_line(3 * i + j + 1, name, 'ALA', seq, 4.31 * i + dx, i_code=i_code)
                      for j, (name, dx) in enumerate((('N', 0.0), ('CA', 1.46), ('C', 2.98)))]
        graph = build_graph(parse_pdb(self.write_pdb(lines)))
        self.assertEqual([r['pos'] for r in graph['residues']], ['52', '52A', '53'])
        self.assertEqual(len(graph['nodes']), 9)
        peptides = [(e['source'], e['target']) for e in graph['edges'] if e['type'] == 'peptide']
        self.assertEqual(peptides, [('A_52_C', 'A_52A_N'), ('A_52A_C', 'A_53_N')])

    def test_conformation_with_highest_occupancy(self):
        path = self.write_pdb(ALTLOC_LINES)
        self.assertEqual(list_altlocs(path), ['A', 'B'])
        residues = parse_pdb(path)
        # Microheterogeneity: one residue, named after the chosen conformation (B, highest occupancy).
        # Residue keys are (chain, number, insertion code, residue name).
        self.assertEqual(list(residues), [('A', 1, '', 'SER')])
        atoms = residues[('A', 1, '', 'SER')]
        self.assertEqual(sorted(atoms), ['C', 'CA', 'N'])
        self.assertEqual((atoms['CA']['alt_loc'], atoms['CA']['x']), ('B', 1.5))
        # For the marker in the sequence bar: amino acid per conformation
        self.assertEqual(build_graph(residues)['residues'][0]['conformations'], {'A': 'PRO', 'B': 'SER'})

    def test_selected_conformation(self):
        path = self.write_pdb(ALTLOC_LINES)
        self.assertEqual(list(parse_pdb(path, altloc='A')), [('A', 1, '', 'PRO')])
        self.assertEqual(parse_pdb(path, altloc='A')[('A', 1, '', 'PRO')]['CA']['x'], 1.46)
        # Residues without the selected conformation use the one with the highest occupancy
        self.assertEqual(list(parse_pdb(path, altloc='C')), [('A', 1, '', 'SER')])

    def test_reused_number_without_altloc_stays_separate(self):
        # A water molecule with the same residue number but without altloc is not merged into the alanine
        lines = residue_lines(1, 0.0) + [pdb_line(4, 'O', 'HOH', 1, 9.0, record='HETATM')]
        self.assertEqual(list(parse_pdb(self.write_pdb(lines))), [('A', 1, '', 'ALA'), ('A', 1, '', 'HOH')])


class ConformationViewTests(BaseTestCase):
    """Conformation selector on the detail page and the ``altloc`` parameter of the model_data view."""

    def setUp(self):
        super().setUp()
        self.pdb = self.make_pdb(personal=True)
        # Replace the default content with the two-conformation fixture
        self.pdb.file.save('altloc.pdb', ContentFile(''.join(ALTLOC_LINES).encode() + b'END\n'))
        self.client.login(username='uploader', password='pw')

    def model_data(self, **params):
        """Requests the JSON data of a model.

        Args:
            **params: Query parameters, e.g. ``altloc`` or ``ensemble``.

        Returns:
            The response of the model_data view.
        """
        return self.client.get(reverse('pdb_grapher:model_data', args=[self.pdb.pk]), params)

    def test_selector_lists_conformations(self):
        response = self.client.get(reverse('pdb_grapher:detail', args=[self.pdb.pk]))
        self.assertContains(response, '<option value="">Konformation: automatisch</option>', html=True)
        self.assertContains(response, '<option value="B">Konformation B</option>', html=True)

    def test_selector_disabled_without_conformations(self):
        pdb = self.make_pdb(personal=True)
        response = self.client.get(reverse('pdb_grapher:detail', args=[pdb.pk]))
        self.assertContains(response, '<option>Eine Konformation</option>', html=True)

    def test_model_data_for_conformation(self):
        # Ensemble statistics are only included when requested with ensemble=1
        data = self.model_data(altloc='A', ensemble='1').json()
        self.assertEqual([r['name'] for r in data['graph']['residues']], ['PRO'])
        self.assertIn('A_1_N|A_1_CA', data['ensemble']['bonds'])
        auto = self.model_data().json()
        self.assertEqual([r['name'] for r in auto['graph']['residues']], ['SER'])
        self.assertNotIn('ensemble', auto)

    def test_unknown_conformation_is_rejected(self):
        self.assertEqual(self.model_data(altloc='Z').status_code, 400)


def xyz_line(serial, name, res_name, res_seq, x, y, z=0.0):
    """Formats an ATOM line in chain A with arbitrary coordinates (for angle and torsion tests).

    Args:
        serial: Atom serial number (columns 7-11).
        name: Atom name of up to three characters (from column 14 on); its first letter is the element.
        res_name: Residue name (columns 18-20).
        res_seq: Residue number (columns 23-26).
        x: x coordinate in Å.
        y: y coordinate in Å.
        z: z coordinate in Å.

    Returns:
        The PDB line including the trailing newline.
    """
    return (f"ATOM  {serial:5d}  {name:<3} {res_name:>3} A{res_seq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {name[0]}\n")


def het_line(serial, name, res_name, x, y=0.0, element=None):
    """Formats a HETATM line of residue 1 in chain A (modified residues, ligands, water, ions).

    Unlike the other helpers, the atom name occupies all four columns 13-16 (left-aligned), and the
    element symbol is right-aligned in columns 77-78, so two-letter elements such as SE or FE fit.

    Args:
        serial: Atom serial number (columns 7-11).
        name: Atom name of up to four characters, e.g. 'SE' or 'CH3'.
        res_name: Residue name, right-aligned in columns 18-20 (e.g. ' MG').
        x: x coordinate in Å.
        y: y coordinate in Å; z is always 0.
        element: Element symbol; defaults to the first letter of ``name``.

    Returns:
        The PDB line including the trailing newline.
    """
    return (f"HETATM{serial:5d} {name:<4} {res_name:>3} A   1    "
            f"{x:8.3f}{y:8.3f}{0:8.3f}  1.00  0.00          {element or name[0]:>2}\n")


# Selenomethionine (MSE, modified methionine) with realistic distances
# (N-CA 1.46, CA-CB 1.53, CB-CG 1.52, CG-SE 1.95, SE-CE 1.95 Å); MSE_SIDECHAIN lists its expected side-chain bonds
MSE_LINES = [
    het_line(1, 'N', 'MSE', 0.0),
    het_line(2, 'CA', 'MSE', 1.46),
    het_line(3, 'C', 'MSE', 1.96, 1.42),
    het_line(4, 'O', 'MSE', 1.36, 2.46),
    het_line(5, 'CB', 'MSE', 1.96, -1.44),
    het_line(6, 'CG', 'MSE', 3.48, -1.44),
    het_line(7, 'SE', 'MSE', 5.43, -1.44, element='SE'),
    het_line(8, 'CE', 'MSE', 7.38, -1.44),
]
MSE_SIDECHAIN = {('CA', 'CB'), ('CB', 'CG'), ('CG', 'SE'), ('SE', 'CE')}


class ModifiedResidueTests(TestCase):
    """Bond detection for modified residues and ligands (HETATM), by distance and from CONECT records."""

    def write_pdb(self, lines):
        """Writes PDB lines to a temporary file that is removed after the test.

        Args:
            lines: PDB lines including newlines.

        Returns:
            Path of the temporary file.
        """
        fd, path = tempfile.mkstemp(suffix='.pdb')
        with os.fdopen(fd, 'w') as f:
            f.writelines(lines)
        self.addCleanup(os.remove, path)
        return path

    def sidechain(self, graph):
        """Returns the side-chain bonds of a graph as a set of (atom name, atom name) pairs."""
        return {(e['source'].split('_')[-1], e['target'].split('_')[-1]) for e in graph['edges'] if e['type'] == 'sidechain'}

    def test_bonds_by_distance_and_parent_from_table(self):
        # Without MODRES the parent amino acid (MET) comes from the built-in table of modified residues
        graph = build_graph(parse_pdb(self.write_pdb(MSE_LINES)))
        self.assertEqual(self.sidechain(graph), MSE_SIDECHAIN)
        self.assertEqual((graph['residues'][0]['parent'], graph['residues'][0]['one_letter']), ('MET', 'M'))

    def test_bonds_from_conect_and_parent_from_modres(self):
        # A deliberately too long CG-SE distance (5.52 Å) is still recognised as a bond via CONECT;
        # the parent amino acid comes from MODRES (here deliberately different from the table: CYS).
        lines = ['MODRES 1ABC MSE A    1  CYS  TEST\n'] + MSE_LINES[:6] + [
            het_line(7, 'SE', 'MSE', 9.0, -1.44, element='SE'), het_line(8, 'CE', 'MSE', 10.95, -1.44)]
        lines += [f'CONECT{a:5d}{b:5d}\n' for a, b in ((2, 5), (5, 6), (6, 7), (7, 8))]
        graph = build_graph(parse_pdb(self.write_pdb(lines)))
        self.assertEqual(self.sidechain(graph), MSE_SIDECHAIN)
        self.assertEqual((graph['residues'][0]['parent'], graph['residues'][0]['one_letter']), ('CYS', 'C'))

    def ligand_edges(self, lines):
        """Builds the graph of the given PDB lines.

        Args:
            lines: PDB lines including newlines.

        Returns:
            Tuple of the graph and the list of its edges as (source, target, type, reference value).
        """
        graph = build_graph(parse_pdb(self.write_pdb(lines)))
        return graph, [(e['source'], e['target'], e['type'], e['ref']) for e in graph['edges']]

    def test_ligand_bonds_by_distance(self):
        # C1-C2 bonded, C3 too far away; ligands have neither reference values nor a parent amino acid
        lines = [het_line(1, 'C1', 'LIG', 0.0), het_line(2, 'C2', 'LIG', 1.5), het_line(3, 'C3', 'LIG', 5.0)]
        graph, edges = self.ligand_edges(lines)
        self.assertEqual(edges, [('A_1_C1', 'A_1_C2', 'ligand', None)])
        self.assertIsNone(graph['residues'][0]['parent'])

    def test_ligand_bonds_from_conect_include_metal(self):
        # Fe-N at 2.0 Å is found by distance thanks to the covalent radius of iron (no CONECT needed);
        # at 2.6 Å, beyond 1.32 + 0.71 + 0.45 Å tolerance, it is only found via CONECT
        lines = [het_line(1, 'FE', 'HEM', 0.0, element='FE'), het_line(2, 'NA', 'HEM', 2.0, element='N')]
        self.assertEqual(self.ligand_edges(lines)[1], [('A_1_FE', 'A_1_NA', 'ligand', None)])
        conect = ['CONECT    1    2\n', 'CONECT    2    1\n']
        far = [het_line(1, 'FE', 'HEM', 0.0, element='FE'), het_line(2, 'NA', 'HEM', 2.6, element='N')]
        self.assertEqual(self.ligand_edges(far + conect)[1], [('A_1_FE', 'A_1_NA', 'ligand', None)])

    def test_ligand_atoms_named_like_backbone_get_no_backbone_bonds(self):
        # The acetyl cap ACE has atoms named C and O but no backbone: C-O only once, as a ligand bond
        lines = [het_line(1, 'C', 'ACE', 0.0), het_line(2, 'O', 'ACE', 1.23, element='O'), het_line(3, 'CH3', 'ACE', -1.5)]
        self.assertEqual(sorted(self.ligand_edges(lines)[1]),
                         [('A_1_C', 'A_1_CH3', 'ligand', None), ('A_1_C', 'A_1_O', 'ligand', None)])

    def test_water_and_ions_have_no_bonds(self):
        # Only 1 Å apart, but water and single ions are never bonded
        lines = [het_line(1, 'O', 'HOH', 0.0, element='O'), het_line(2, 'MG', ' MG', 1.0, element='MG')]
        self.assertEqual(self.ligand_edges(lines)[1], [])


# Header records with ID code, experimental method, resolution and R values (R-work, R-free)
HEADER_LINES = [
    'HEADER    HYDROLASE                               01-JAN-00   9XYZ              \n',
    'EXPDTA    X-RAY DIFFRACTION                                                     \n',
    'REMARK   2 RESOLUTION.    1.80 ANGSTROMS.                                       \n',
    'REMARK   3   R VALUE            (WORKING SET) : 0.193                           \n',
    'REMARK   3   FREE R VALUE                     : 0.234                           \n',
]


class ReferenceValueTests(TestCase):
    """Reference values for bond lengths and angles (edge ``ref``, node ``angle_refs``) and ``parse_header()``.

    A reference value is ``[mean, standard deviation, key]``; the key names the source table entry,
    e.g. ``'backbone:N-CA:GLY'`` (Engh & Huber 2001) or ``'PRO:CA-CB'`` (CCP4 Monomer Library).
    """

    def write_pdb(self, lines):
        """Writes PDB lines to a temporary file that is removed after the test.

        Args:
            lines: PDB lines including newlines.

        Returns:
            Path of the temporary file.
        """
        fd, path = tempfile.mkstemp(suffix='.pdb')
        with os.fdopen(fd, 'w') as f:
            f.writelines(lines)
        self.addCleanup(os.remove, path)
        return path

    def test_reference_classes(self):
        # Glycine before proline: backbone values from Engh & Huber 2001 (Gly class, peptide bond before Pro),
        # proline side chain from the CCP4 Monomer Library
        lines = [pdb_line(i + 1, name, res, seq, x) for i, (name, res, seq, x) in enumerate((
            ('N', 'GLY', 1, 0.0), ('CA', 'GLY', 1, 1.451), ('C', 'GLY', 1, 2.967), ('O', 'GLY', 1, 4.198),
            ('N', 'PRO', 2, 4.308), ('CA', 'PRO', 2, 5.774), ('C', 'PRO', 2, 7.299), ('CB', 'PRO', 2, 9.0),
        ))]
        edges = {(e['source'], e['target']): e['ref'] for e in build_graph(parse_pdb(self.write_pdb(lines)))['edges']}
        self.assertEqual(edges[('A_1_N', 'A_1_CA')], [1.456, 0.015, 'backbone:N-CA:GLY'])
        self.assertEqual(edges[('A_1_CA', 'A_1_C')], [1.514, 0.016, 'backbone:CA-C:GLY'])
        self.assertEqual(edges[('A_1_C', 'A_1_O')], [1.232, 0.016, 'backbone:C-O:GLY'])
        self.assertEqual(edges[('A_1_C', 'A_2_N')], [1.338, 0.019, 'backbone:C-N:X-PRO'])
        self.assertEqual(edges[('A_2_N', 'A_2_CA')], [1.468, 0.017, 'backbone:N-CA:PRO'])
        self.assertEqual(edges[('A_2_CA', 'A_2_CB')], [1.526, 0.016, 'PRO:CA-CB'])

    def test_modified_residue_uses_parent_for_backbone_and_own_sidechain(self):
        # MSE: backbone from the parent class (MET falls into 'other'), side chain from its own MSE entry
        edges = {(e['source'], e['target']): e['ref'] for e in build_graph(parse_pdb(self.write_pdb(MSE_LINES)))['edges']}
        self.assertEqual(edges[('A_1_N', 'A_1_CA')][2], 'backbone:N-CA:other')
        self.assertEqual(edges[('A_1_CA', 'A_1_CB')], [1.532, 0.01, 'MSE:CA-CB'])
        self.assertEqual(edges[('A_1_CG', 'A_1_SE')], [1.963, 0.02, 'MSE:CG-SE'])

    def test_every_standard_sidechain_bond_has_a_reference(self):
        for res_name, bonds in SIDECHAIN_BONDS.items():
            for a, b in bonds:
                with self.subTest(res_name=res_name, bond=(a, b)):
                    self.assertIn((res_name, frozenset((a, b))), SIDECHAIN_BOND_REFERENCES)

    def gly_pro(self, cis):
        """Builds Gly 1 and Pro 2 in the plane z = 0 and returns the angle references per atom.

        Args:
            cis: Whether the peptide bond is cis (ω = CA-C-N-CA = 0°) instead of trans (180°).

        Returns:
            Dict node id -> ``angle_refs`` of that node (key 'A|B' for the angle A-node-B, sorted
            atom names), empty for nodes without references.
        """
        atoms = (
            ('N', 'GLY', 1, -1.2, 1.9), ('CA', 'GLY', 1, 0.0, 1.0), ('C', 'GLY', 1, 0.0, 0.0),
            ('O', 'GLY', 1, -1.0, -0.7), ('N', 'PRO', 2, 1.33, 0.0), ('CA', 'PRO', 2, 1.9, 1.3 if cis else -1.3),
            ('C', 'PRO', 2, 3.4, 1.3 if cis else -1.3), ('CD', 'PRO', 2, 2.0, -1.2 if cis else 1.2),
        )
        lines = [xyz_line(i + 1, name, res, seq, x, y) for i, (name, res, seq, x, y) in enumerate(atoms)]
        return {n['id']: n.get('angle_refs', {}) for n in build_graph(parse_pdb(self.write_pdb(lines)))['nodes']}

    def test_angle_references(self):
        refs = self.gly_pro(cis=False)
        self.assertEqual(refs['A_1_CA']['C|N'], [113.1, 2.5, 'backbone:N-CA-C:GLY'])
        self.assertEqual(refs['A_1_C']['CA|O'], [120.6, 1.8, 'backbone:CA-C-O:GLY'])
        # Across the peptide bond: for CA-C-N the preceding Gly takes precedence, for O-C-N the following Pro
        self.assertEqual(refs['A_1_C']['CA|N'], [116.2, 2.0, 'backbone:CA-C-N:GLY-X'])
        self.assertEqual(refs['A_1_C']['N|O'], [121.1, 1.9, 'backbone:O-C-N:X-PRO'])
        self.assertEqual(refs['A_2_N']['C|CA'], [119.3, 1.5, 'backbone:C-N-CA:X-PRO'])
        # Proline ring from the CCP4 Monomer Library; C-N-CD has no reference
        self.assertEqual(refs['A_2_N']['CA|CD'], [108.0, 1.5, 'PRO:CA-N-CD'])
        self.assertNotIn('C|CD', refs['A_2_N'])

    def test_cis_proline_angle_reference(self):
        self.assertEqual(self.gly_pro(cis=True)['A_2_N']['C|CA'], [127.0, 2.4, 'backbone:C-N-CA:cis-X-PRO'])

    def test_no_angle_references_without_peptide_bond_or_for_ligands(self):
        # Residue 2 starts 20 Å away (no peptide bond), so A_1_C has only one bonded neighbour (CA) and
        # no angle, while A_2_CA keeps N-CA-C; the ligand angle C1-C2-C3 gets no reference
        lines = residue_lines(1, 0.0) + residue_lines(2, 20.0, serial=4) + [
            het_line(7, 'C1', 'LIG', 40.0), het_line(8, 'C2', 'LIG', 41.5), het_line(9, 'C3', 'LIG', 43.0)]
        refs = {n['id']: n.get('angle_refs', {}) for n in build_graph(parse_pdb(self.write_pdb(lines)))['nodes']}
        self.assertEqual(set(refs['A_1_C']), set())
        self.assertEqual(set(refs['A_2_CA']), {'C|N'})

    def test_parse_header(self):
        # All fields are None for a file without header records
        path = self.write_pdb(HEADER_LINES + residue_lines(1, 0.0))
        self.assertEqual(parse_header(path), {
            'id_code': '9XYZ', 'method': 'X-RAY DIFFRACTION', 'resolution': 1.8, 'r_work': 0.193, 'r_free': 0.234,
        })
        self.assertEqual(parse_header(self.write_pdb(residue_lines(1, 0.0))), {
            'id_code': None, 'method': None, 'resolution': None, 'r_work': None, 'r_free': None,
        })


class StructureInfoViewTests(BaseTestCase):
    """Header metadata and reference value sources in the statistics of the detail page."""

    def test_metadata_in_statistics(self):
        # German UI: translated method name and decimal comma
        pdb = self.make_pdb(personal=True)
        pdb.file.save('meta.pdb', ContentFile(''.join(HEADER_LINES).encode() + PDB_CONTENT))
        self.client.login(username='uploader', password='pw')
        response = self.client.get(reverse('pdb_grapher:detail', args=[pdb.pk]))
        self.assertContains(response, '<dd>Röntgenkristallographie</dd>', html=True)
        self.assertContains(response, '<dd>1,80 Å</dd>', html=True)
        self.assertContains(response, '<dd>0,234</dd>', html=True)
        self.assertContains(response, 'Engh &amp; Huber (2001)')
        self.assertContains(response, 'CCP4 Monomer Library')
        self.assertContains(response, 'id="angle-ref-groups"')


def cys_lines(serial, res_seq, x, sg_x):
    """Formats a cysteine with N, CA, C, CB on the x axis starting at ``x`` and SG at (sg_x, 2, 0).

    CB is placed on top of CA; its position is irrelevant because the tests only look at SG-SG
    distances. SG is moved off the axis by overwriting the y field (columns 39-46) of its line.

    Args:
        serial: Atom serial number of N; the other atoms get the next four numbers.
        res_seq: Residue number.
        x: x coordinate of N in Å.
        sg_x: x coordinate of SG in Å.

    Returns:
        List of five PDB lines.
    """
    atoms = [('N', x), ('CA', x + 1.46), ('C', x + 2.98), ('CB', x + 1.46)]
    lines = [pdb_line(serial + i, name, 'CYS', res_seq, ax) for i, (name, ax) in enumerate(atoms)]
    sg = pdb_line(serial + 4, 'SG', 'CYS', res_seq, sg_x)
    return lines + [sg[:38] + f'{2.0:8.3f}' + sg[46:]]


class DisulfideTests(TestCase):
    """Disulfide bridges between cysteines, by SG-SG distance or from SSBOND records."""

    def write_pdb(self, lines):
        """Writes PDB lines to a temporary file that is removed after the test.

        Args:
            lines: PDB lines including newlines.

        Returns:
            Path of the temporary file.
        """
        fd, path = tempfile.mkstemp(suffix='.pdb')
        with os.fdopen(fd, 'w') as f:
            f.writelines(lines)
        self.addCleanup(os.remove, path)
        return path

    def bridges(self, lines):
        """Returns the disulfide edges as (source, target, length, via), via being 'distance' or 'ssbond'."""
        graph = build_graph(parse_pdb(self.write_pdb(lines)))
        return [(e['source'], e['target'], e['length'], e['via']) for e in graph['edges'] if e['type'] == 'disulfide']

    def test_bridge_by_distance(self):
        # SG-SG 2.04 Å between cysteines 1 and 2; cysteine 3 is far away
        lines = cys_lines(1, 1, 0.0, 10.0) + cys_lines(6, 2, 20.0, 12.04) + cys_lines(11, 3, 40.0, 40.0)
        self.assertEqual(self.bridges(lines), [('A_1_SG', 'A_2_SG', 2.04, 'distance')])

    def test_bridge_from_ssbond_regardless_of_distance(self):
        # SG-SG is 3.0 Å, above MAX_DISULFIDE_LENGTH (2.5 Å); symmetry operator 1555 (identity) for both
        ssbond = 'SSBOND   1 CYS A    1    CYS A    3                          1555   1555  3.00  \n'
        lines = [ssbond] + cys_lines(1, 1, 0.0, 10.0) + cys_lines(6, 3, 20.0, 13.0)
        self.assertEqual(self.bridges(lines), [('A_1_SG', 'A_3_SG', 3.0, 'ssbond')])

    def test_ssbond_to_symmetry_copy_is_ignored(self):
        # Symmetry operator 2655 for the second partner: the bridge goes to a symmetry copy not in the file
        ssbond = 'SSBOND   1 CYS A    1    CYS A    3                          1555   2655  2.05  \n'
        lines = [ssbond] + cys_lines(1, 1, 0.0, 10.0) + cys_lines(6, 3, 20.0, 13.0)
        self.assertEqual(self.bridges(lines), [])


class AnglesTabTests(BaseTestCase):
    """Presence of the angles/Ramachandran tab, its controls and export buttons on the detail page."""

    def test_angles_tab_and_export_buttons(self):
        pdb = self.make_pdb(personal=True)
        self.client.login(username='uploader', password='pw')
        response = self.client.get(reverse('pdb_grapher:detail', args=[pdb.pk]))
        self.assertContains(response, 'id="tab-ramachandran"')
        self.assertContains(response, 'id="angle-bond-types"')
        self.assertContains(response, 'id="angle-rama"')
        for fmt in ('png', 'jpeg', 'svg'):
            self.assertContains(response, f'data-ramachandran-format="{fmt}"')
        self.assertContains(response, 'id="rama-type"')
        self.assertContains(response, 'id="export-pdf-bonds"')
        self.assertContains(response, 'rama8000_regions.json')
        self.assertContains(response, 'Williams et al. (2018)')


def _place(a, b, c, length, angle, torsion):
    """Places atom d from three preceding atoms using internal coordinates (NeRF method).

    Args:
        a: Coordinates [x, y, z] of atom a.
        b: Coordinates of atom b.
        c: Coordinates of atom c, to which d is bonded.
        length: Bond length |cd| in Å.
        angle: Bond angle b-c-d in degrees.
        torsion: Torsion (dihedral) angle a-b-c-d in degrees.

    Returns:
        Coordinates [x, y, z] of atom d.
    """
    sub = lambda u, v: [u[k] - v[k] for k in range(3)]  # noqa: E731
    cross = lambda u, v: [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]  # noqa: E731
    unit = lambda u: [x / math.sqrt(sum(y * y for y in u)) for x in u]  # noqa: E731
    bc = unit(sub(c, b))
    n = unit(cross(sub(b, a), bc))
    m = cross(n, bc)
    angle, torsion = math.radians(angle), math.radians(torsion)
    d2 = [-length * math.cos(angle), length * math.sin(angle) * math.cos(torsion), length * math.sin(angle) * math.sin(torsion)]
    return [c[k] + d2[0] * bc[k] + d2[1] * m[k] + d2[2] * n[k] for k in range(3)]


def peptide_lines(residues):
    """Formats the backbone (N, CA, C) of a chain with ideal geometry and the given torsion angles.

    Bond lengths and angles are fixed (N-CA 1.458, CA-C 1.525, C-N 1.329 Å; N-CA-C 111.0°,
    CA-C-N 116.2°, C-N-CA 121.7°), so the computed φ/ψ/ω match the requested values.

    Args:
        residues: List of tuples (residue name, φ, ψ, ω) in degrees. ω is the torsion of the peptide
            bond *before* the residue (180 = trans, 0 = cis). The φ of the first residue, the ψ of the
            last residue and the ω of the first residue are not used.

    Returns:
        List of PDB lines, three per residue, numbered from 1.
    """
    atoms = [[0.0, 0.0, 0.0], [1.458, 0.0, 0.0]]
    atoms.append(_place([0.0, 1.0, 0.0], atoms[0], atoms[1], 1.525, 111.0, -60.0))
    for i, (_, phi, psi, _) in enumerate(residues[:-1]):
        omega = residues[i + 1][3]
        n_, ca, c = atoms[-3:]
        nxt = _place(n_, ca, c, 1.329, 116.2, psi)
        nxt_ca = _place(ca, c, nxt, 1.458, 121.7, omega)
        nxt_c = _place(c, nxt, nxt_ca, 1.525, 111.0, residues[i + 1][1])
        atoms += [nxt, nxt_ca, nxt_c]
    lines = []
    for i, (name, _, _, _) in enumerate(residues):
        for k, atom in enumerate(('N', 'CA', 'C')):
            x, y, z = atoms[3 * i + k]
            lines.append(xyz_line(3 * i + k + 1, atom, name, i + 1, x, y, z))
    return lines


class AveragedModelTests(TestCase):
    """Angles of the averaged model are means over the models, torsions circular means."""

    def write_pdb(self, lines):
        """Writes PDB lines to a temporary file that is removed after the test; returns its path."""
        fd, path = tempfile.mkstemp(suffix='.pdb')
        with os.fdopen(fd, 'w') as f:
            f.writelines(lines)
        self.addCleanup(os.remove, path)
        return path

    def test_angles_are_averaged_over_the_models(self):
        # The two models differ only in ψ of residue 2: 179° and −179°. Their circular mean is
        # 180°; averaging coordinates or plain numbers would give a wrong value (0° for numbers).
        models = [peptide_lines([('ALA', 0, 150, 180), ('ALA', -60, psi, 180), ('ALA', -60, 0, 180)])
                  for psi in (179, -179)]
        path = self.write_pdb(['MODEL        1\n'] + models[0] + ['ENDMDL\n', 'MODEL        2\n'] + models[1]
                              + ['ENDMDL\n'])
        graph = build_average_graph(path, [1, 2])
        rama = graph['residues'][1]['rama']
        self.assertAlmostEqual(abs(rama['psi']), 180, delta=0.1)
        self.assertAlmostEqual(rama['phi'], -60, delta=0.1)
        # Bond angles are identical in both models and keep the ideal value of peptide_lines()
        self.assertAlmostEqual(graph['averaged']['bond_angles']['A_1_C|A_2_N|A_2_CA'], 121.7, delta=0.01)
        self.assertAlmostEqual(abs(graph['averaged']['torsions']['A_1_CA|A_1_C|A_2_N|A_2_CA']), 180, delta=0.1)


class RamachandranTests(TestCase):
    """φ/ψ angles per residue, Ramachandran categories and regions (Top8000 distributions, MolProbity)."""

    def write_pdb(self, lines):
        """Writes PDB lines to a temporary file that is removed after the test.

        Args:
            lines: PDB lines including newlines.

        Returns:
            Path of the temporary file.
        """
        fd, path = tempfile.mkstemp(suffix='.pdb')
        with os.fdopen(fd, 'w') as f:
            f.writelines(lines)
        self.addCleanup(os.remove, path)
        return path

    def rama(self, residues):
        """Builds a chain with ``peptide_lines(residues)`` and returns the ``rama`` entry of each residue.

        Args:
            residues: List of tuples (residue name, φ, ψ, ω), see ``peptide_lines()``.

        Returns:
            List with one dict (keys include 'phi', 'psi', 'type', 'region') or None per residue.
        """
        graph = build_graph(parse_pdb(self.write_pdb(peptide_lines(residues))))
        return [r['rama'] for r in graph['residues']]

    def test_values_and_regions(self):
        # α-helix favoured, left-handed helix favoured for glycine, (100°, -100°) an outlier
        self.assertEqual(ramachandran_region('general', ramachandran_value('general', -63, -43)), 'favored')
        self.assertEqual(ramachandran_region('general', ramachandran_value('general', 100, -100)), 'outlier')
        self.assertEqual(ramachandran_region('glycine', ramachandran_value('glycine', 80, 10)), 'favored')
        # Density thresholds per category as in MolProbity (general 0.0005, cis-proline 0.002, ile-val 0.001)
        self.assertEqual(ramachandran_region('general', 0.0005), 'allowed')
        self.assertEqual(ramachandran_region('cis-proline', 0.0015), 'outlier')
        self.assertEqual(ramachandran_region('ile-val', 0.001), 'allowed')

    def test_wrapping_interpolation_is_continuous(self):
        # φ = -180° and φ = 180° are the same angle; the interpolation wraps around the grid edge
        self.assertAlmostEqual(ramachandran_value('general', -180, 150), ramachandran_value('general', 180, 150))

    def test_phi_psi_and_types(self):
        helix = (-63.0, -43.0)
        rama = self.rama([('ALA', *helix, 180), ('GLY', *helix, 180), ('VAL', *helix, 180), ('LEU', *helix, 180),
                          ('PRO', *helix, 180), ('SER', *helix, 180), ('PRO', *helix, 0), ('ALA', *helix, 180),
                          ('ALA', *helix, 180)])
        self.assertIsNone(rama[0])  # no φ
        self.assertIsNone(rama[-1])  # no ψ
        self.assertAlmostEqual(rama[1]['phi'], -63.0, places=1)
        self.assertAlmostEqual(rama[1]['psi'], -43.0, places=1)
        # GLY -> glycine, VAL -> ile-val, LEU before PRO -> pre-proline, PRO with trans ω -> trans-proline,
        # SER before PRO -> pre-proline, PRO with ω = 0 -> cis-proline, ALA -> general
        self.assertEqual([r['type'] for r in rama[1:-1]],
                         ['glycine', 'ile-val', 'pre-proline', 'trans-proline', 'pre-proline', 'cis-proline', 'general'])

    def test_outlier_and_chain_gap(self):
        rama = self.rama([('ALA', -63, -43, 180), ('ALA', 100, -100, 180), ('ALA', -63, -43, 180)])
        self.assertEqual(rama[1]['region'], 'outlier')
        lines = peptide_lines([('ALA', -63, -43, 180)] * 3)
        # Shift residue 3 (residue number in columns 23-26) by 30 Å in x (columns 31-38) to create a gap
        moved = [line[:30] + f'{float(line[30:38]) + 30:8.3f}' + line[38:] if int(line[22:26]) == 3 else line for line in lines]
        graph = build_graph(parse_pdb(self.write_pdb(moved)))
        self.assertIsNone(graph['residues'][1]['rama'])  # chain gap after residue 2, so no ψ
