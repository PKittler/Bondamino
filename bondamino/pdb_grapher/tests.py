import os
import shutil
import tempfile
import urllib.error
from unittest import mock

from django.conf import settings
from django.contrib.auth import get_user_model
from django.core.files.base import ContentFile
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase, override_settings
from django.urls import reverse

from .models import Language, Note, PDBFile, WorkGroup
from .utils import (
    build_average_graph, build_bond_ensemble, build_graph, list_altlocs, parse_header, parse_pdb,
)

PDB_CONTENT = (
    b"ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N\n"
    b"ATOM      2  CA  ALA A   1      11.639   6.071  -5.147  1.00  0.00           C\n"
    b"ATOM      3  C   ALA A   1      13.161   6.195  -5.161  1.00  0.00           C\n"
    b"ATOM      4  O   ALA A   1      13.773   5.880  -6.191  1.00  0.00           O\n"
    b"END\n"
)

MEDIA_ROOT = tempfile.mkdtemp()
User = get_user_model()


@override_settings(MEDIA_ROOT=MEDIA_ROOT)
class BaseTestCase(TestCase):
    """Users, Arbeitsgruppen und Hilfsmethoden für alle Tests."""

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()
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
        pdb = PDBFile(name='test.pdb', uploaded_by=self.uploader, personal=personal)
        pdb.file.save('test.pdb', ContentFile(PDB_CONTENT), save=False)
        pdb.save()
        pdb.groups.set(groups)
        return pdb


class PermissionTests(BaseTestCase):
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
        personal = self.make_pdb(personal=True)
        shared = self.make_pdb(personal=True, groups=[self.group_a])
        self.uploader.delete()
        self.assertFalse(PDBFile.objects.filter(pk=personal.pk).exists())
        shared.refresh_from_db()
        self.assertIsNone(shared.uploaded_by)
        self.assertFalse(shared.personal)

    def test_index_lists_files_with_actions(self):
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
    def setUp(self):
        super().setUp()
        self.pdb = self.make_pdb(groups=[self.group_a, self.group_b])

    def make_note(self, author, groups=()):
        note = Note.objects.create(pdb=self.pdb, author=author, text=f'Notiz von {author}')
        note.groups.set(groups)
        return note

    def post_note(self, **data):
        return self.client.post(reverse('pdb_grapher:note_create', args=[self.pdb.pk]),
                                {f'note-new-{key}': value for key, value in data.items()})

    def notes_on_detail_page(self, username):
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
        self.client.login(username='member', password='pw')
        self.assertEqual(self.post_note(text='x', visibility='groups').status_code, 400)
        self.assertEqual(self.post_note(text='x', visibility='groups', groups=[self.group_b.pk]).status_code, 400)
        self.assertEqual(self.post_note(text='', visibility='personal').status_code, 400)
        self.assertFalse(Note.objects.exists())

    def test_visibility(self):
        personal = self.make_note(self.member)
        shared_a = self.make_note(self.member, groups=[self.group_a])
        shared_b = self.make_note(self.outsider, groups=[self.group_b])
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
        for username in ('member', 'outsider'):  # Verfasser ohne Leitung bzw. Leiter einer anderen Gruppe
            self.client.login(username=username, password='pw')
            self.assertEqual(self.client.post(url).status_code, 403)
        self.client.login(username='leader', password='pw')
        self.assertEqual(self.client.post(url).status_code, 200)
        self.assertEqual(list(note.groups.all()), [self.group_b])

    def test_removal_from_last_group_keeps_note_for_author(self):
        note = self.make_note(self.member, groups=[self.group_a])
        orphan = self.make_note(self.member, groups=[self.group_a])
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
        note = self.make_note(self.member, groups=[self.group_a])
        Note.objects.filter(pk=note.pk).update(author=None)
        kept = self.make_note(self.member, groups=[self.group_a])
        self.group_a.delete()
        self.assertFalse(Note.objects.filter(pk=note.pk).exists())
        self.assertTrue(Note.objects.filter(pk=kept.pk).exists())


class IndexTests(BaseTestCase):
    def index_files(self, **params):
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
        hidden = Note.objects.create(pdb=pdb, author=self.member, text='b')  # nur für den Verfasser
        Note.objects.filter(pk=hidden.pk).update(created_at=visible.created_at.replace(year=2099))
        [row] = self.index_files()
        self.assertEqual(row.latest_note_at, visible.created_at)

    def test_failed_upload_reopens_dialog(self):
        self.client.login(username='uploader', password='pw')
        response = self.client.post(reverse('pdb_grapher:index'), {'pdb_file': SimpleUploadedFile('x.pdb', PDB_CONTENT)})
        self.assertTrue(response.context['show_upload'])
        self.assertFalse(self.client.get(reverse('pdb_grapher:index')).context['show_upload'])


class RCSBImportTests(BaseTestCase):
    def post_upload(self, **data):
        self.client.login(username='uploader', password='pw')
        return self.client.post(reverse('pdb_grapher:index'), {'personal': 'on', **data})

    def mock_rcsb(self, content=PDB_CONTENT, error=None):
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
    def setUp(self):
        super().setUp()
        self.client.login(username='uploader', password='pw')

    def index(self, **extra):
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
        self.assertContains(response, 'Log out')
        self.assertContains(response, 'All available')

    def test_only_languages_from_admin(self):
        response = self.client.post(reverse('pdb_grapher:set_language'), {'language': 'fr'})
        self.assertNotIn(settings.LANGUAGE_COOKIE_NAME, response.cookies)
        self.assertContains(self.index(HTTP_ACCEPT_LANGUAGE='fr'), '<html lang="de"')
        self.client.cookies[settings.LANGUAGE_COOKIE_NAME] = 'en-gb'
        Language.objects.filter(code='en-gb').delete()
        self.assertContains(self.index(), '<html lang="de"')

    def test_browser_language_matches_variant(self):
        self.assertContains(self.index(HTTP_ACCEPT_LANGUAGE='en-US,en;q=0.8'), '<html lang="en-gb"')

    def test_unsafe_next_url_is_ignored(self):
        response = self.client.post(reverse('pdb_grapher:set_language'),
                                    {'language': 'en-gb', 'next': 'https://example.com/'})
        self.assertRedirects(response, reverse('pdb_grapher:index'), fetch_redirect_response=False)

    def test_selector_lists_admin_languages(self):
        self.assertContains(self.index(), '<option value="en-gb" lang="en-gb" >English</option>', html=False)
        Language.objects.filter(code='en-gb').delete()
        self.assertNotContains(self.index(), 'name="language"')

    def test_javascript_catalog(self):
        self.client.cookies[settings.LANGUAGE_COOKIE_NAME] = 'en-gb'
        response = self.client.get(reverse('javascript-catalog'))
        self.assertContains(response, 'Aspartic acid')

    def test_server_messages_are_translated(self):
        self.client.cookies[settings.LANGUAGE_COOKIE_NAME] = 'en-gb'
        pdb = self.make_pdb(personal=True)
        response = self.client.post(reverse('pdb_grapher:delete', args=[pdb.pk]), follow=True)
        self.assertContains(response, '“test.pdb” has been deleted.')


def atom_line(serial, name, res_seq, x):
    return (f"ATOM  {serial:5d}  {name:<3} ALA A{res_seq:4d}    "
            f"{x:8.3f}{0:8.3f}{0:8.3f}  1.00  0.00           {name[0]}\n")


def residue_lines(res_seq, start, serial=1):
    """Backbone-Atome N, CA, C auf einer Linie; das N des nächsten Residuums liegt bei start + 4.31,
    also 1,33 Å hinter C."""
    return [atom_line(serial + i, name, res_seq, start + dx)
            for i, (name, dx) in enumerate((('N', 0.0), ('CA', 1.46), ('C', 2.98)))]


class ChainGapTests(TestCase):
    def write_pdb(self, lines):
        fd, path = tempfile.mkstemp(suffix='.pdb')
        with os.fdopen(fd, 'w') as f:
            f.writelines(lines)
        self.addCleanup(os.remove, path)
        return path

    def peptide_bonds(self, graph):
        return {(e['source'], e['target']): e['length'] for e in graph['edges'] if e['type'] == 'peptide'}

    def test_gap_is_no_peptide_bond(self):
        # Residuum 2 folgt direkt auf 1, vor Residuum 5 fehlen zwei Residuen (C–N-Abstand ≈ 10 Å)
        path = self.write_pdb(residue_lines(1, 0.0) + residue_lines(2, 4.31, 4) + residue_lines(5, 18.0, 7))
        bonds = self.peptide_bonds(build_graph(parse_pdb(path)))
        self.assertEqual(bonds, {('A_1_C', 'A_2_N'): 1.33})

    def test_gap_in_one_model_is_ignored(self):
        # Model 1: Residuen 1 und 2 verbunden; Model 2: Residuum 2 liegt weit entfernt
        lines = (['MODEL        1\n'] + residue_lines(1, 0.0) + residue_lines(2, 4.31, 4) + ['ENDMDL\n'] +
                 ['MODEL        2\n'] + residue_lines(1, 0.0) + residue_lines(2, 14.0, 4) + ['ENDMDL\n'])
        path = self.write_pdb(lines)
        self.assertEqual(self.peptide_bonds(build_average_graph(path, [1, 2])), {('A_1_C', 'A_2_N'): 1.33})
        ensemble = build_bond_ensemble(path, [1, 2])
        self.assertEqual(ensemble['bonds']['A_1_C|A_2_N'], [1.33, 1.33, 1])
        self.assertEqual(ensemble['max']['length'], 1.52)


def pdb_line(serial, name, res_name, res_seq, x, alt=' ', i_code=' ', occupancy=1.0, record='ATOM  '):
    return (f"{record}{serial:5d}  {name:<3}{alt}{res_name:>3} A{res_seq:4d}{i_code}   "
            f"{x:8.3f}{0:8.3f}{0:8.3f}{occupancy:6.2f}  0.00           {name[0]}\n")


# Residuum 1: CA in zwei Konformationen, A (Prolin, Belegung 0,3) und B (Serin, Belegung 0,7)
ALTLOC_LINES = [
    pdb_line(1, 'N', 'PRO', 1, 0.0),
    pdb_line(2, 'CA', 'PRO', 1, 1.46, alt='A', occupancy=0.3),
    pdb_line(3, 'CA', 'SER', 1, 1.50, alt='B', occupancy=0.7),
    pdb_line(4, 'C', 'PRO', 1, 2.98),
]


class ConformationParsingTests(TestCase):
    def write_pdb(self, lines):
        fd, path = tempfile.mkstemp(suffix='.pdb')
        with os.fdopen(fd, 'w') as f:
            f.writelines(lines)
        self.addCleanup(os.remove, path)
        return path

    def test_insertion_codes_are_separate_residues(self):
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
        # Mikroheterogenität: ein Residuum, Name aus der gewählten Konformation
        self.assertEqual(list(residues), [('A', 1, '', 'SER')])
        atoms = residues[('A', 1, '', 'SER')]
        self.assertEqual(sorted(atoms), ['C', 'CA', 'N'])
        self.assertEqual((atoms['CA']['alt_loc'], atoms['CA']['x']), ('B', 1.5))
        # Für die Markierung in der Sequenzleiste: Aminosäure je Konformation
        self.assertEqual(build_graph(residues)['residues'][0]['conformations'], {'A': 'PRO', 'B': 'SER'})

    def test_selected_conformation(self):
        path = self.write_pdb(ALTLOC_LINES)
        self.assertEqual(list(parse_pdb(path, altloc='A')), [('A', 1, '', 'PRO')])
        self.assertEqual(parse_pdb(path, altloc='A')[('A', 1, '', 'PRO')]['CA']['x'], 1.46)
        # Residuen ohne die gewählte Konformation nehmen die mit der höchsten Belegung
        self.assertEqual(list(parse_pdb(path, altloc='C')), [('A', 1, '', 'SER')])

    def test_reused_number_without_altloc_stays_separate(self):
        lines = residue_lines(1, 0.0) + [pdb_line(4, 'O', 'HOH', 1, 9.0, record='HETATM')]
        self.assertEqual(list(parse_pdb(self.write_pdb(lines))), [('A', 1, '', 'ALA'), ('A', 1, '', 'HOH')])


class ConformationViewTests(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.pdb = self.make_pdb(personal=True)
        self.pdb.file.save('altloc.pdb', ContentFile(''.join(ALTLOC_LINES).encode() + b'END\n'))
        self.client.login(username='uploader', password='pw')

    def model_data(self, **params):
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
        data = self.model_data(altloc='A', ensemble='1').json()
        self.assertEqual([r['name'] for r in data['graph']['residues']], ['PRO'])
        self.assertIn('A_1_N|A_1_CA', data['ensemble']['bonds'])
        auto = self.model_data().json()
        self.assertEqual([r['name'] for r in auto['graph']['residues']], ['SER'])
        self.assertNotIn('ensemble', auto)

    def test_unknown_conformation_is_rejected(self):
        self.assertEqual(self.model_data(altloc='Z').status_code, 400)


def het_line(serial, name, res_name, x, y=0.0, element=None):
    return (f"HETATM{serial:5d} {name:<4} {res_name:>3} A   1    "
            f"{x:8.3f}{y:8.3f}{0:8.3f}  1.00  0.00          {element or name[0]:>2}\n")


# Selenomethionin mit realistischen Abständen (N–CA 1,46, CA–CB 1,53, CB–CG 1,52, CG–SE 1,95, SE–CE 1,95 Å)
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
    def write_pdb(self, lines):
        fd, path = tempfile.mkstemp(suffix='.pdb')
        with os.fdopen(fd, 'w') as f:
            f.writelines(lines)
        self.addCleanup(os.remove, path)
        return path

    def sidechain(self, graph):
        return {(e['source'].split('_')[-1], e['target'].split('_')[-1]) for e in graph['edges'] if e['type'] == 'sidechain'}

    def test_bonds_by_distance_and_parent_from_table(self):
        graph = build_graph(parse_pdb(self.write_pdb(MSE_LINES)))
        self.assertEqual(self.sidechain(graph), MSE_SIDECHAIN)
        self.assertEqual((graph['residues'][0]['parent'], graph['residues'][0]['one_letter']), ('MET', 'M'))

    def test_bonds_from_conect_and_parent_from_modres(self):
        # Ein absichtlich zu langer CG–SE-Abstand wird über CONECT trotzdem als Bindung erkannt,
        # die Ausgangs-Aminosäure kommt aus MODRES (hier abweichend von der Tabelle).
        lines = ['MODRES 1ABC MSE A    1  CYS  TEST\n'] + MSE_LINES[:6] + [
            het_line(7, 'SE', 'MSE', 9.0, -1.44, element='SE'), het_line(8, 'CE', 'MSE', 10.95, -1.44)]
        lines += [f'CONECT{a:5d}{b:5d}\n' for a, b in ((2, 5), (5, 6), (6, 7), (7, 8))]
        graph = build_graph(parse_pdb(self.write_pdb(lines)))
        self.assertEqual(self.sidechain(graph), MSE_SIDECHAIN)
        self.assertEqual((graph['residues'][0]['parent'], graph['residues'][0]['one_letter']), ('CYS', 'C'))

    def test_ligands_stay_without_bonds(self):
        lines = [het_line(1, 'C1', 'LIG', 0.0), het_line(2, 'C2', 'LIG', 1.5)]
        graph = build_graph(parse_pdb(self.write_pdb(lines)))
        self.assertEqual(graph['edges'], [])
        self.assertIsNone(graph['residues'][0]['parent'])


HEADER_LINES = [
    'HEADER    HYDROLASE                               01-JAN-00   9XYZ              \n',
    'EXPDTA    X-RAY DIFFRACTION                                                     \n',
    'REMARK   2 RESOLUTION.    1.80 ANGSTROMS.                                       \n',
    'REMARK   3   R VALUE            (WORKING SET) : 0.193                           \n',
    'REMARK   3   FREE R VALUE                     : 0.234                           \n',
]


class ReferenceValueTests(TestCase):
    def write_pdb(self, lines):
        fd, path = tempfile.mkstemp(suffix='.pdb')
        with os.fdopen(fd, 'w') as f:
            f.writelines(lines)
        self.addCleanup(os.remove, path)
        return path

    def test_reference_classes(self):
        # Glycin vor Prolin: N–CA und CA–C der Gly-Klasse, Peptidbindung mit Prolin-Wert
        lines = [pdb_line(i + 1, name, res, seq, x) for i, (name, res, seq, x) in enumerate((
            ('N', 'GLY', 1, 0.0), ('CA', 'GLY', 1, 1.451), ('C', 'GLY', 1, 2.967), ('O', 'GLY', 1, 4.198),
            ('N', 'PRO', 2, 4.308), ('CA', 'PRO', 2, 5.774), ('C', 'PRO', 2, 7.299), ('CB', 'PRO', 2, 9.0),
        ))]
        edges = {(e['source'], e['target']): e['ref'] for e in build_graph(parse_pdb(self.write_pdb(lines)))['edges']}
        self.assertEqual(edges[('A_1_N', 'A_1_CA')], [1.451, 0.016, 'N:CA:GLY'])
        self.assertEqual(edges[('A_1_CA', 'A_1_C')], [1.516, 0.018, 'CA:C:GLY'])
        self.assertEqual(edges[('A_1_C', 'A_1_O')], [1.231, 0.020, 'C:O:all'])
        self.assertEqual(edges[('A_1_C', 'A_2_N')], [1.341, 0.016, 'C:N:PRO'])
        self.assertEqual(edges[('A_2_N', 'A_2_CA')], [1.466, 0.015, 'N:CA:PRO'])
        self.assertEqual(edges[('A_2_CA', 'A_2_CB')], [1.530, 0.020, 'CA:CB:other'])

    def test_modified_residue_uses_parent_and_sidechain_has_no_reference(self):
        edges = {(e['source'], e['target']): e['ref'] for e in build_graph(parse_pdb(self.write_pdb(MSE_LINES)))['edges']}
        self.assertEqual(edges[('A_1_N', 'A_1_CA')][2], 'N:CA:other')
        self.assertEqual(edges[('A_1_CA', 'A_1_CB')][2], 'CA:CB:other')
        self.assertIsNone(edges[('A_1_CG', 'A_1_SE')])

    def test_parse_header(self):
        path = self.write_pdb(HEADER_LINES + residue_lines(1, 0.0))
        self.assertEqual(parse_header(path), {
            'id_code': '9XYZ', 'method': 'X-RAY DIFFRACTION', 'resolution': 1.8, 'r_work': 0.193, 'r_free': 0.234,
        })
        self.assertEqual(parse_header(self.write_pdb(residue_lines(1, 0.0))), {
            'id_code': None, 'method': None, 'resolution': None, 'r_work': None, 'r_free': None,
        })


class StructureInfoViewTests(BaseTestCase):
    def test_metadata_in_statistics(self):
        pdb = self.make_pdb(personal=True)
        pdb.file.save('meta.pdb', ContentFile(''.join(HEADER_LINES).encode() + PDB_CONTENT))
        self.client.login(username='uploader', password='pw')
        response = self.client.get(reverse('pdb_grapher:detail', args=[pdb.pk]))
        self.assertContains(response, '<dd>Röntgenkristallographie</dd>', html=True)
        self.assertContains(response, '<dd>1,80 Å</dd>', html=True)
        self.assertContains(response, '<dd>0,234</dd>', html=True)
        self.assertContains(response, 'Engh &amp; Huber (1991)')
