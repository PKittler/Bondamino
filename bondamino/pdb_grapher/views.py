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

AVERAGE_MODEL = 'avg'  # Pseudo-Model: Mittelwert der Bindungslängen über alle Models


def _method_label(method):
    """Lesbare Bezeichnung der Strukturbestimmungsmethode (EXPDTA), mehrere mit "; " getrennt."""
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
            except Exception:
                # Details nur ins Server-Log, nicht an den User
                logger.exception('Could not process uploaded file %s', pdb.name)
                pdb.delete()
                form.add_error('pdb_file', _('Die Datei konnte nicht verarbeitet werden.'))
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


def _safe_next(request):
    """Ziel aus dem Feld `next`, falls es auf diese Seite zeigt, sonst die Startseite."""
    next_url = request.POST.get('next', '')
    if not url_has_allowed_host_and_scheme(next_url, allowed_hosts={request.get_host()},
                                           require_https=request.is_secure()):
        next_url = reverse('pdb_grapher:index')
    return next_url


def _set_language_cookie(response, code):
    response.set_cookie(
        settings.LANGUAGE_COOKIE_NAME, code,
        max_age=settings.LANGUAGE_COOKIE_AGE, path=settings.LANGUAGE_COOKIE_PATH,
        domain=settings.LANGUAGE_COOKIE_DOMAIN, secure=settings.LANGUAGE_COOKIE_SECURE,
        httponly=settings.LANGUAGE_COOKIE_HTTPONLY, samesite=settings.LANGUAGE_COOKIE_SAMESITE,
    )


@require_POST
def set_language(request):
    """Merkt sich die gewählte Sprache im Cookie und, wenn angemeldet, in den Einstellungen des
    Users. Nur im Admin angelegte Sprachen sind gültig."""
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
    """Darstellung (hell/dunkel) in den Einstellungen speichern. Per fetch ohne Weiterleitung."""
    theme = request.POST.get('theme')
    if theme in dict(UserSettings.THEME_CHOICES):
        UserSettings.objects.update_or_create(user=request.user, defaults={'theme': theme})
    elif request.headers.get('X-Requested-With') == 'fetch':
        return JsonResponse({'error': 'invalid theme'}, status=400)
    if request.headers.get('X-Requested-With') == 'fetch':
        return JsonResponse({'theme': UserSettings.of(request.user).theme})
    return redirect(_safe_next(request))


class PasswordResetView(auth_views.PasswordResetView):
    """Passwort vergessen: Link per E-Mail (auch für User, die ihre Einladung nicht angenommen haben)."""
    form_class = PasswordResetForm
    template_name = 'pdb_grapher/password_reset_form.html'
    subject_template_name = 'pdb_grapher/email/password_reset_subject.txt'
    email_template_name = 'pdb_grapher/email/password_reset_email.txt'
    success_url = reverse_lazy('pdb_grapher:password_reset_done')
    extra_email_context = {'valid_hours': settings.PASSWORD_RESET_TIMEOUT // 3600}


class PasswordResetDoneView(auth_views.PasswordResetDoneView):
    template_name = 'pdb_grapher/password_reset_done.html'
    extra_context = {'valid_hours': settings.PASSWORD_RESET_TIMEOUT // 3600}


class SetPasswordView(auth_views.PasswordResetConfirmView):
    """Passwort über den Link aus der Einladung (`invitation`) oder aus "Passwort vergessen"
    festlegen; danach ist der User angemeldet. Der Link gilt nur einmal."""
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
    """Konto des Users: Angaben (vom Admin gepflegt), Passwort ändern, Einstellungen."""
    user_settings = UserSettings.of(request.user)
    action = request.POST.get('action') if request.method == 'POST' else None
    password_form = PasswordChangeForm(request.user, request.POST if action == 'password' else None)
    settings_form = UserSettingsForm(request.POST if action == 'settings' else None, instance=user_settings,
                                     current_language=request.LANGUAGE_CODE)
    for field in password_form.fields.values():
        field.widget.attrs['class'] = 'form-control'

    if action == 'password' and password_form.is_valid():
        password_form.save()
        # Sitzung bleibt nach der Änderung gültig
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
    pdb = _visible_pdb_or_404(request.user, pk)
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
    except Exception:
        logger.exception('Could not process file pk %s (model %s, altloc %s)', pk, model_param, altloc)
        return JsonResponse({'error': _('Die Datei konnte nicht verarbeitet werden.')}, status=400)

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
