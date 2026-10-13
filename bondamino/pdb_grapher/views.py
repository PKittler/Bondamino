import json

from django.conf import settings
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.db import transaction
from django.db.models import OuterRef, Q, Subquery
from django.http import FileResponse, HttpResponseForbidden, JsonResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse
from django.utils import timezone
from django.utils.http import url_has_allowed_host_and_scheme
from django.utils.translation import gettext as _
from django.views.decorators.http import require_POST

from .forms import AvailabilityForm, NoteForm, PDBUploadForm
from .models import Language, Note, PDBFile, WorkGroup
from .utils import parse_pdb, build_graph, build_average_graph, build_bond_ensemble, list_altlocs, list_models

AVERAGE_MODEL = 'avg'  # Pseudo-Model: Mittelwert der Bindungslängen über alle Models


def _visible_pdb_or_404(user, pk):
    return get_object_or_404(PDBFile.objects.visible_to(user), pk=pk)


@login_required
def index(request):
    """Upload neuer PDB-Dateien und Liste aller für den User verfügbaren Dateien."""
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
            except Exception as e:
                pdb.delete()
                form.add_error('pdb_file', _('Fehler beim Verarbeiten der Datei: %(error)s') % {'error': e})
            else:
                return redirect('pdb_grapher:detail', pk=pdb.pk)

    user_groups = WorkGroup.objects.filter(members=request.user)
    # Datum der neuesten Anmerkung, die der User sehen darf
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
        'show_upload': form.is_bound,  # Upload-Dialog mit Fehlermeldungen wieder öffnen
        'pdb_files': pdb_files,
        'user_groups': user_groups,
        'selected_group': selected_group,
        'query': query,
        'led_group_ids': set(request.user.led_workgroups.values_list('pk', flat=True)),
    })


@require_POST
def set_language(request):
    """Merkt sich die gewählte Sprache im Cookie. Nur im Admin angelegte Sprachen sind gültig."""
    next_url = request.POST.get('next', '')
    if not url_has_allowed_host_and_scheme(next_url, allowed_hosts={request.get_host()},
                                           require_https=request.is_secure()):
        next_url = reverse('pdb_grapher:index')
    response = redirect(next_url)
    code = request.POST.get('language')
    if Language.objects.filter(code=code).exists():
        response.set_cookie(
            settings.LANGUAGE_COOKIE_NAME, code,
            max_age=settings.LANGUAGE_COOKIE_AGE, path=settings.LANGUAGE_COOKIE_PATH,
            domain=settings.LANGUAGE_COOKIE_DOMAIN, secure=settings.LANGUAGE_COOKIE_SECURE,
            httponly=settings.LANGUAGE_COOKIE_HTTPONLY, samesite=settings.LANGUAGE_COOKIE_SAMESITE,
        )
    return response


@login_required
def detail(request, pk):
    pdb = _visible_pdb_or_404(request.user, pk)
    file_path = pdb.file.path
    try:
        models = list_models(file_path)
        altlocs = list_altlocs(file_path)
        selected_model = models[0] if models else None
        graph_data = build_graph(parse_pdb(file_path, model_id=selected_model))
        bond_ensemble = build_bond_ensemble(file_path, models)
    except Exception as e:
        messages.error(request, _('Fehler beim Verarbeiten der Datei „%(name)s“: %(error)s')
                       % {'name': pdb.name, 'error': e})
        return redirect('pdb_grapher:index')
    return render(request, 'pdb_grapher/graph.html', {
        'pdb': pdb,
        'graph_data': json.dumps(graph_data),
        'filename': pdb.name,
        'models': models,
        'selected_model': selected_model,
        'average_model': AVERAGE_MODEL,
        'altlocs': altlocs,
        'bond_ensemble': bond_ensemble,
        'note_form': NoteForm(user=request.user, prefix='note-new'),
        'user_groups': WorkGroup.objects.filter(members=request.user),
        **_notes_context(request, pdb),
    })


@login_required
def model_data(request, pk):
    """Liefert die Graph-Daten für ein bestimmtes Model und eine Konformation derselben,
    bereits hochgeladenen PDB-Datei (per AJAX von den Auswahlfeldern verwendet).

    Parameter: `model` (Nummer oder AVERAGE_MODEL), `altloc` (leer = je Residuum die
    Konformation mit der höchsten Belegung) und `ensemble=1`, um zusätzlich die
    Bindungsstatistik über alle Models mitzuliefern (ändert sich mit der Konformation).
    Antwort: {"graph": ..., "ensemble": ...}.
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
    except Exception as e:
        return JsonResponse({'error': _('Fehler beim Verarbeiten der Datei: %(error)s') % {'error': e}}, status=400)

    return JsonResponse(data)


@login_required
def download(request, pk):
    """Liefert die PDB-Datei als Download aus – nur für User, die sie sehen dürfen."""
    pdb = _visible_pdb_or_404(request.user, pk)
    return FileResponse(pdb.file.open('rb'), as_attachment=True, filename=pdb.name)


def _notes_context(request, pdb):
    """Für den User sichtbare Anmerkungen samt Berechtigungen und Bearbeitungsformularen."""
    notes = list(pdb.notes.visible_to(request.user)
                 .select_related('author').prefetch_related('groups'))
    for note in notes:
        note.user_can_edit = note.editable_by(request.user)
        if note.user_can_edit:
            note.edit_form = NoteForm(instance=note, user=request.user, prefix=f'note-{note.pk}')
    return {
        'notes': notes,
        'led_group_ids': set(request.user.led_workgroups.values_list('pk', flat=True)),
    }


def _render_notes(request, pdb):
    """Aktualisierte Anmerkungsliste (per AJAX vom Anmerkungen-Tab verwendet)."""
    return render(request, 'pdb_grapher/_note_list.html', _notes_context(request, pdb))


def _form_errors(form):
    return JsonResponse({'errors': [e for errors in form.errors.values() for e in errors]}, status=400)


@login_required
@require_POST
def note_create(request, pk):
    pdb = _visible_pdb_or_404(request.user, pk)
    form = NoteForm(request.POST, user=request.user, prefix='note-new')
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
    """Anmerkung bearbeiten – nur für den Verfasser."""
    pdb = _visible_pdb_or_404(request.user, pk)
    note = get_object_or_404(pdb.notes.visible_to(request.user), pk=note_pk)
    if not note.editable_by(request.user):
        return HttpResponseForbidden(_('Nur der Verfasser darf die Anmerkung bearbeiten.'))
    form = NoteForm(request.POST, instance=note, user=request.user, prefix=f'note-{note.pk}')
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
    """Anmerkung löschen – nur für den Verfasser."""
    pdb = _visible_pdb_or_404(request.user, pk)
    note = get_object_or_404(pdb.notes.visible_to(request.user), pk=note_pk)
    if not note.editable_by(request.user):
        return HttpResponseForbidden(_('Nur der Verfasser darf die Anmerkung löschen.'))
    note.delete()
    return _render_notes(request, pdb)


@login_required
@require_POST
def note_remove_from_group(request, pk, note_pk, group_pk):
    """Gruppenleiter entfernt eine Anmerkung aus seiner Gruppe. Für den Verfasser bleibt
    sie erhalten; hat sie keinen Verfasser und ist für keine Gruppe mehr lesbar, wird sie
    gelöscht."""
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
    """Verfügbarkeit bearbeiten – nur für den Uploader."""
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
    """Datei endgültig löschen – nur für den Uploader."""
    pdb = _visible_pdb_or_404(request.user, pk)
    if not pdb.can_manage(request.user):
        return HttpResponseForbidden(_('Nur der Uploader darf die Datei löschen.'))
    pdb.delete()
    messages.success(request, _('„%(name)s“ wurde gelöscht.') % {'name': pdb.name})
    return redirect('pdb_grapher:index')


@login_required
@require_POST
def remove_from_group(request, pk, group_pk):
    """Gruppenleiter entfernt eine Datei aus seiner Gruppe. Ist sie danach nirgends
    mehr verfügbar, wird sie gelöscht."""
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
