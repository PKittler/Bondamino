"""Views of Bondamino (URL routes in ``urls.py``).

Request flow:

* ``index`` lists the visible files and handles uploads (local file or PDB ID fetched from RCSB, see
  ``forms.PDBUploadForm``). A new upload is test-parsed once and discarded if parsing fails.
* ``detail`` parses the file with the helpers from ``utils.py`` and renders ``graph.html``, which draws the graph
  in the browser. Switching model or alternate conformation (altloc) later fetches JSON from ``model_data``.
* The notes tab posts to ``note_create``/``note_edit``/``note_delete``/``note_remove_from_group``; these return the
  re-rendered note list (``_note_list.html``) as an HTML fragment, or JSON errors with status 400.
* ``edit``/``delete``/``remove_from_group`` manage a file's availability; ``account``, ``set_language`` and
  ``set_theme`` manage user settings; the password views handle invitations and forgotten passwords.

Permission model:

* Every view except ``set_language`` and the password views requires login.
* Files and notes are always looked up through ``visible_to`` (see ``models.py``), so a user who may not see an
  object gets a 404 rather than a 403 and cannot probe for existing IDs.
* Changing an object additionally checks ownership (``PDBFile.can_manage``, ``Note.editable_by``) or group
  leadership (``WorkGroup.leader``) and answers 403 otherwise.
"""
import logging

from django.conf import settings
from django.contrib import messages
from django.contrib.auth import update_session_auth_hash
from django.contrib.auth import views as auth_views
from django.contrib.auth.decorators import login_required
from django.contrib.auth.forms import PasswordChangeForm
from django.db import transaction
from django.db.models import OuterRef, Q, Subquery
from django.http import FileResponse, HttpResponseForbidden, JsonResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse, reverse_lazy
from django.utils import timezone
from django.utils.http import url_has_allowed_host_and_scheme
from django.utils.translation import gettext as _
from django.views.decorators.http import require_POST

from .forms import AvailabilityForm, NoteForm, PasswordResetForm, PDBUploadForm, UserSettingsForm
from .models import Language, Note, PDBFile, UserSettings, WorkGroup
from .utils import (
    RAMACHANDRAN_SOURCE, REFERENCE_SOURCES, build_average_graph, build_bond_ensemble, build_graph, list_altlocs, list_models,
    parse_header, parse_pdb,
)

logger = logging.getLogger(__name__)

AVERAGE_MODEL = 'avg'  # Pseudo model: bond lengths averaged over all models


def _method_label(method):
    """Translates the structure determination method (EXPDTA record) into a readable label.

    Args:
        method: Raw EXPDTA value; several methods are separated by ";". Unknown methods are kept as they are.

    Returns:
        The translated labels joined with "; ", or None if ``method`` is empty.
    """
    if not method:
        return None
    labels = {
        'X-RAY DIFFRACTION': _('Röntgenkristallographie'),
        'NEUTRON DIFFRACTION': _('Neutronenbeugung'),
        'ELECTRON CRYSTALLOGRAPHY': _('Elektronenkristallographie'),
        'ELECTRON MICROSCOPY': _('Elektronenmikroskopie'),
        'SOLUTION NMR': _('NMR (Lösung)'),
        'SOLID-STATE NMR': _('Festkörper-NMR'),
    }
    return '; '.join(labels.get(m.strip(), m.strip()) for m in method.split(';'))


def _visible_pdb_or_404(user, pk):
    """Returns the PDB file ``pk`` if ``user`` may see it.

    Raises:
        Http404: If the file does not exist or is not visible to the user.
    """
    return get_object_or_404(PDBFile.objects.visible_to(user), pk=pk)


@login_required
def index(request):
    """Start page: list of all files visible to the user, plus the upload dialog.

    GET parameters ``group`` (``'personal'`` for the user's own uploads or a group pk) and ``q`` (substring of the
    file name) filter the list. A POST uploads a new file (``PDBUploadForm``); the file is test-parsed and, on
    success, the user is redirected to its detail page. On failure the page is re-rendered with the upload
    dialog open and the errors shown.
    """
    form = PDBUploadForm(user=request.user)
    if request.method == 'POST':
        form = PDBUploadForm(request.POST, request.FILES, user=request.user)
        if form.is_valid():
            pdb = form.save(commit=False)
            pdb.file = form.cleaned_data['pdb_file']
            pdb.name = pdb.file.name
            pdb.uploaded_by = request.user
            pdb.save()
            form.save_m2m()
            try:
                models = list_models(pdb.file.path)
                parse_pdb(pdb.file.path, model_id=models[0] if models else None)
            except Exception:
                # Details go to the server log only, not to the user
                logger.exception('Could not process uploaded file %s', pdb.name)
                pdb.delete()
                form.add_error('pdb_file', _('Die Datei konnte nicht verarbeitet werden.'))
            else:
                return redirect('pdb_grapher:detail', pk=pdb.pk)

    user_groups = WorkGroup.objects.filter(members=request.user)
    # Date of the newest note the user may see (same condition as Note.objects.visible_to)
    latest_note = (Note.objects.filter(pdb=OuterRef('pk'))
                   .filter(Q(author=request.user) | Q(groups__members=request.user))
                   .order_by('-created_at').values('created_at')[:1])
    pdb_files = (PDBFile.objects.visible_to(request.user)
                 .annotate(latest_note_at=Subquery(latest_note))
                 .prefetch_related('groups'))
    selected_group = request.GET.get('group')
    if selected_group == 'personal':
        pdb_files = pdb_files.filter(uploaded_by=request.user)
    elif selected_group and selected_group.isdigit():
        pdb_files = pdb_files.filter(groups__pk=selected_group)
    query = request.GET.get('q', '').strip()
    if query:
        pdb_files = pdb_files.filter(name__icontains=query)
    return render(request, 'pdb_grapher/index.html', {
        'form': form,
        'show_upload': form.is_bound,  # reopen the upload dialog showing the errors
        'pdb_files': pdb_files,
        'user_groups': user_groups,
        'selected_group': selected_group,
        'query': query,
        'led_group_ids': set(request.user.led_workgroups.values_list('pk', flat=True)),
    })


def _safe_next(request):
    """Returns the redirect target from the POST field ``next``.

    Prevents open redirects: the target is used only if it points to this host.

    Returns:
        The ``next`` URL if it is safe, otherwise the URL of the start page.
    """
    next_url = request.POST.get('next', '')
    if not url_has_allowed_host_and_scheme(next_url, allowed_hosts={request.get_host()},
                                           require_https=request.is_secure()):
        next_url = reverse('pdb_grapher:index')
    return next_url


def _set_language_cookie(response, code):
    """Sets Django's language cookie on ``response`` (read by ``middleware.select_language``)."""
    response.set_cookie(
        settings.LANGUAGE_COOKIE_NAME, code,
        max_age=settings.LANGUAGE_COOKIE_AGE, path=settings.LANGUAGE_COOKIE_PATH,
        domain=settings.LANGUAGE_COOKIE_DOMAIN, secure=settings.LANGUAGE_COOKIE_SECURE,
        httponly=settings.LANGUAGE_COOKIE_HTTPONLY, samesite=settings.LANGUAGE_COOKIE_SAMESITE,
    )


@require_POST
def set_language(request):
    """Stores the selected language (POST only, login not required).

    POST parameters: ``language`` (code of a language created in the admin; other values are ignored) and
    ``next`` (where to go back to). The language is stored in the cookie and, if the user is logged in, also in
    the user's settings.

    Returns:
        A redirect to ``next`` (or the start page).
    """
    response = redirect(_safe_next(request))
    code = request.POST.get('language')
    if Language.objects.filter(code=code).exists():
        _set_language_cookie(response, code)
        if request.user.is_authenticated:
            UserSettings.objects.update_or_create(user=request.user, defaults={'language': code})
    return response


@login_required
@require_POST
def set_theme(request):
    """Stores the theme (light/dark) in the user's settings (POST only).

    POST parameters: ``theme`` and ``next``. Requests sent by the theme switch via fetch (header
    ``X-Requested-With: fetch``) get JSON instead of a redirect.

    Returns:
        For fetch: ``{"theme": ...}``, or ``{"error": ...}`` with status 400 for an unknown theme.
        Otherwise a redirect to ``next`` (an unknown theme is silently ignored).
    """
    theme = request.POST.get('theme')
    if theme in dict(UserSettings.THEME_CHOICES):
        UserSettings.objects.update_or_create(user=request.user, defaults={'theme': theme})
    elif request.headers.get('X-Requested-With') == 'fetch':
        return JsonResponse({'error': 'invalid theme'}, status=400)
    if request.headers.get('X-Requested-With') == 'fetch':
        return JsonResponse({'theme': UserSettings.of(request.user).theme})
    return redirect(_safe_next(request))


class PasswordResetView(auth_views.PasswordResetView):
    """Forgotten password: sends a link by email.

    Also works for users who never accepted their invitation; rate limited (see ``forms.PasswordResetForm``).
    """
    form_class = PasswordResetForm
    template_name = 'pdb_grapher/password_reset_form.html'
    subject_template_name = 'pdb_grapher/email/password_reset_subject.txt'
    email_template_name = 'pdb_grapher/email/password_reset_email.txt'
    success_url = reverse_lazy('pdb_grapher:password_reset_done')
    extra_email_context = {'valid_hours': settings.PASSWORD_RESET_TIMEOUT // 3600}


class PasswordResetDoneView(auth_views.PasswordResetDoneView):
    """Confirmation page after requesting a password reset, including how long the link is valid."""
    template_name = 'pdb_grapher/password_reset_done.html'
    extra_context = {'valid_hours': settings.PASSWORD_RESET_TIMEOUT // 3600}


class SetPasswordView(auth_views.PasswordResetConfirmView):
    """Sets a password via the link from an invitation or from "forgot password".

    The same view serves both URLs; ``invitation=True`` (set in ``urls.py``) only changes the texts of the
    template. Afterwards the user is logged in. The link works only once because the token depends on the
    password hash.

    Attributes:
        invitation: True for the invitation link, False for the password reset link.
    """
    template_name = 'pdb_grapher/password_set.html'
    post_reset_login = True
    success_url = reverse_lazy('pdb_grapher:index')
    invitation = False

    def get_form(self, form_class=None):
        form = super().get_form(form_class)
        for field in form.fields.values():
            field.widget.attrs['class'] = 'form-control'
        return form

    def form_valid(self, form):
        response = super().form_valid(form)
        messages.success(self.request, _('Dein Passwort wurde gespeichert.'))
        return response

    def get_context_data(self, **kwargs):
        return super().get_context_data(invitation=self.invitation, **kwargs)


@login_required
def account(request):
    """Account page: user data (maintained by the admin), password change and settings.

    The page contains two forms; the POST field ``action`` (``'password'`` or ``'settings'``) tells which one was
    submitted. After a successful save the user is redirected back to this page (post/redirect/get); invalid
    forms are re-rendered with their errors.
    """
    user_settings = UserSettings.of(request.user)
    action = request.POST.get('action') if request.method == 'POST' else None
    password_form = PasswordChangeForm(request.user, request.POST if action == 'password' else None)
    settings_form = UserSettingsForm(request.POST if action == 'settings' else None, instance=user_settings,
                                     current_language=request.LANGUAGE_CODE)
    for field in password_form.fields.values():
        field.widget.attrs['class'] = 'form-control'

    if action == 'password' and password_form.is_valid():
        password_form.save()
        # Keep the current session valid after the change (Django would otherwise log the user out)
        update_session_auth_hash(request, password_form.user)
        messages.success(request, _('Das Passwort wurde geändert.'))
        return redirect('pdb_grapher:account')
    if action == 'settings' and settings_form.is_valid():
        saved = settings_form.save()
        response = redirect('pdb_grapher:account')
        _set_language_cookie(response, saved.language)
        messages.success(request, _('Die Einstellungen wurden gespeichert.'))
        return response

    return render(request, 'pdb_grapher/account.html', {
        'password_form': password_form,
        'settings_form': settings_form,
        'groups': WorkGroup.objects.filter(members=request.user),
    })


@login_required
def detail(request, pk):
    """Graph page of a PDB file (requires the file to be visible to the user, otherwise 404).

    Parses header, models and altlocs and builds the graph of the first model plus the bond statistics over all
    models. If the file cannot be processed, the user is redirected to the start page with an error message.
    """
    pdb =_visible_pdb_or_404(request.user, pk)
    file_path = pdb.file.path
    try:
        models = list_models(file_path)
        altlocs = list_altlocs(file_path)
        structure_info = parse_header(file_path)
        structure_info.update(method_label=_method_label(structure_info['method']), reference_sources=REFERENCE_SOURCES,
                              ramachandran_source=RAMACHANDRAN_SOURCE)
        selected_model = models[0] if models else None
        graph_data = build_graph(parse_pdb(file_path, model_id=selected_model))
        bond_ensemble = build_bond_ensemble(file_path, models)
    except Exception:
        logger.exception('Could not process file %s (pk %s)', pdb.name, pdb.pk)
        messages.error(request, _('Die Datei „%(name)s“ konnte nicht verarbeitet werden.') % {'name': pdb.name})
        return redirect('pdb_grapher:index')
    return render(request, 'pdb_grapher/graph.html', {
        'pdb': pdb,
        'graph_data': graph_data,
        'filename': pdb.name,
        'models': models,
        'selected_model': selected_model,
        'average_model': AVERAGE_MODEL,
        'altlocs': altlocs,
        'structure_info': structure_info,
        'bond_ensemble': bond_ensemble,
        'note_form': NoteForm(user=request.user, pdb=pdb, prefix='note-new'),
        'user_groups': WorkGroup.objects.filter(members=request.user),
        **_notes_context(request, pdb),
    })


@login_required
def model_data(request, pk):
    """Returns the graph data for a given model and conformation of an already uploaded PDB file.

    Called via AJAX by the model/altloc selectors on the graph page. The file must be visible to the user,
    otherwise 404.

    GET parameters: ``model`` (model number or ``AVERAGE_MODEL``; default: first model), ``altloc`` (empty = for
    each residue the conformation with the highest occupancy) and ``ensemble=1`` to also return the bond
    statistics over all models (they change with the conformation).

    Returns:
        JSON ``{"graph": ..., "ensemble": ...}`` (``ensemble`` only if requested), or ``{"error": ...}`` with
        status 400 for invalid parameters or a file that cannot be processed.
    """
    file_path = _visible_pdb_or_404(request.user, pk).file.path

    altloc = request.GET.get('altloc') or None
    if altloc is not None and altloc not in list_altlocs(file_path):
        return JsonResponse({'error': _('Ungültige Konformation.')}, status=400)

    models = list_models(file_path)
    model_param = request.GET.get('model')
    average = model_param == AVERAGE_MODEL
    if average:
        model_id = None
    elif model_param is not None:
        try:
            model_id = int(model_param)
        except ValueError:
            return JsonResponse({'error': _('Ungültiges Model.')}, status=400)
        if model_id not in models:
            # Would otherwise silently return an empty graph.
            return JsonResponse({'error': _('Ungültiges Model.')}, status=400)
    else:
        model_id = models[0] if models else None

    try:
        if average:
            graph_data = build_average_graph(file_path, models, altloc=altloc)
        else:
            graph_data = build_graph(parse_pdb(file_path, model_id=model_id, altloc=altloc))
        data = {'graph': graph_data}
        if request.GET.get('ensemble'):
            data['ensemble'] = build_bond_ensemble(file_path, models, altloc=altloc)
    except Exception:
        logger.exception('Could not process file pk %s (model %s, altloc %s)', pk, model_param, altloc)
        return JsonResponse({'error': _('Die Datei konnte nicht verarbeitet werden.')}, status=400)

    return JsonResponse(data)


@login_required
def download(request, pk):
    """Sends the PDB file as a download, only to users who may see it (otherwise 404)."""
    pdb = _visible_pdb_or_404(request.user, pk)
    return FileResponse(pdb.file.open('rb'), as_attachment=True, filename=pdb.name)


def _notes_context(request, pdb):
    """Builds the template context for the note list of a file.

    Returns:
        A dict with ``notes`` (notes visible to the user; each with ``user_can_edit`` and, for the author, an
        ``edit_form``) and ``led_group_ids`` (pks of the groups the user leads, for the "remove from group"
        buttons).
    """
    notes = list(pdb.notes.visible_to(request.user)
                 .select_related('author').prefetch_related('groups'))
    for note in notes:
        note.user_can_edit = note.editable_by(request.user)
        if note.user_can_edit:
            note.edit_form = NoteForm(instance=note, user=request.user, pdb=pdb, prefix=f'note-{note.pk}')
    return {
        'notes': notes,
        'led_group_ids': set(request.user.led_workgroups.values_list('pk', flat=True)),
    }


def _render_notes(request, pdb):
    """Renders the updated note list as an HTML fragment (used via AJAX by the notes tab)."""
    return render(request, 'pdb_grapher/_note_list.html', _notes_context(request, pdb))


def _form_errors(form):
    """Returns all error messages of ``form`` as JSON ``{"errors": [...]}`` with status 400."""
    return JsonResponse({'errors': [e for errors in form.errors.values() for e in errors]}, status=400)


@login_required
@require_POST
def note_create(request, pk):
    """Creates a note on a visible file (POST, fields of ``NoteForm`` with prefix ``note-new``).

    Returns:
        The updated note list, or JSON errors with status 400.
    """
    pdb =_visible_pdb_or_404(request.user, pk)
    form = NoteForm(request.POST, user=request.user, pdb=pdb, prefix='note-new')
    if not form.is_valid():
        return _form_errors(form)
    note = form.save(commit=False)
    note.pdb = pdb
    note.author = request.user
    note.save()
    form.save_m2m()
    return _render_notes(request, pdb)


@login_required
@require_POST
def note_edit(request, pk, note_pk):
    """Edits a note; author only (POST, fields of ``NoteForm`` with prefix ``note-<pk>``).

    Returns:
        The updated note list, JSON errors with status 400, or 403 if the user is not the author.
    """
    pdb = _visible_pdb_or_404(request.user, pk)
    note = get_object_or_404(pdb.notes.visible_to(request.user), pk=note_pk)
    if not note.editable_by(request.user):
        return HttpResponseForbidden(_('Nur der Verfasser darf die Anmerkung bearbeiten.'))
    form = NoteForm(request.POST, instance=note, user=request.user, pdb=pdb, prefix=f'note-{note.pk}')
    if not form.is_valid():
        return _form_errors(form)
    note = form.save(commit=False)
    note.edited_at = timezone.now()
    note.save()
    form.save_m2m()
    return _render_notes(request, pdb)


@login_required
@require_POST
def note_delete(request, pk, note_pk):
    """Deletes a note; author only (POST).

    Returns:
        The updated note list, or 403 if the user is not the author.
    """
    pdb = _visible_pdb_or_404(request.user, pk)
    note = get_object_or_404(pdb.notes.visible_to(request.user), pk=note_pk)
    if not note.editable_by(request.user):
        return HttpResponseForbidden(_('Nur der Verfasser darf die Anmerkung löschen.'))
    note.delete()
    return _render_notes(request, pdb)


@login_required
@require_POST
def note_remove_from_group(request, pk, note_pk, group_pk):
    """Removes a note from a group; group leader only (POST).

    The note stays available to its author. If it has no author and is no longer shared with any group, it is
    deleted.

    Returns:
        The updated note list, or 403 if the user is not the leader of the group.
    """
    pdb = _visible_pdb_or_404(request.user, pk)
    group = get_object_or_404(WorkGroup, pk=group_pk)
    if group.leader_id != request.user.pk:
        return HttpResponseForbidden(_('Nur der Gruppenleiter darf Anmerkungen aus der Gruppe entfernen.'))
    note = get_object_or_404(pdb.notes.filter(groups=group), pk=note_pk)
    with transaction.atomic():
        note.groups.remove(group)
        note.delete_if_orphaned()
    return _render_notes(request, pdb)


@login_required
def edit(request, pk):
    """Edits a file's availability (personal and/or groups); uploader only.

    GET shows ``AvailabilityForm``; a valid POST saves it and redirects to the start page. Answers 403 if the user
    is not the uploader.
    """
    pdb = _visible_pdb_or_404(request.user, pk)
    if not pdb.can_manage(request.user):
        return HttpResponseForbidden(_('Nur der Uploader darf die Verfügbarkeit bearbeiten.'))
    form = AvailabilityForm(request.POST or None, instance=pdb, user=request.user)
    if request.method == 'POST' and form.is_valid():
        form.save()
        messages.success(request, _('Verfügbarkeit von „%(name)s“ gespeichert.') % {'name': pdb.name})
        return redirect('pdb_grapher:index')
    return render(request, 'pdb_grapher/edit.html', {'form': form, 'pdb': pdb})


@login_required
@require_POST
def delete(request, pk):
    """Deletes a file permanently, for all groups; uploader only (POST, 403 otherwise).

    Redirects to the start page. The stored file is removed by ``signals.delete_stored_file``.
    """
    pdb = _visible_pdb_or_404(request.user, pk)
    if not pdb.can_manage(request.user):
        return HttpResponseForbidden(_('Nur der Uploader darf die Datei löschen.'))
    pdb.delete()
    messages.success(request, _('„%(name)s“ wurde gelöscht.') % {'name': pdb.name})
    return redirect('pdb_grapher:index')


@login_required
@require_POST
def remove_from_group(request, pk, group_pk):
    """Removes a file from a group; group leader only (POST, 403 otherwise).

    The leader need not be the uploader. If the file is then available nowhere (not personal and in no other
    group), it is deleted. Redirects to the start page.
    """
    group = get_object_or_404(WorkGroup, pk=group_pk)
    if group.leader_id != request.user.pk:
        return HttpResponseForbidden(_('Nur der Gruppenleiter darf Dateien aus der Gruppe entfernen.'))
    pdb = get_object_or_404(group.pdb_files, pk=pk)
    with transaction.atomic():
        pdb.groups.remove(group)
        deleted = pdb.delete_if_orphaned()
    if deleted:
        messages.success(request, _('„%(name)s“ wurde aus „%(group)s“ entfernt und gelöscht, '
                                    'da sie nirgends mehr verfügbar war.') % {'name': pdb.name, 'group': group})
    else:
        messages.success(request, _('„%(name)s“ wurde aus „%(group)s“ entfernt.') % {'name': pdb.name, 'group': group})
    return redirect('pdb_grapher:index')
