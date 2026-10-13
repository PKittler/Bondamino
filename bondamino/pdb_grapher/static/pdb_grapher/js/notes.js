'use strict';

// ===== Anmerkungen-Tab =====
// Anlegen, Bearbeiten und Löschen laufen per AJAX; der Server liefert jeweils die
// aktualisierte Liste zurück, damit die Graph-Seite nicht neu geladen werden muss.
(function () {
    const root = document.getElementById('notes');
    if (!root) return;
    const list = document.getElementById('note-list');
    const filter = document.getElementById('notes-filter');
    const empty = document.getElementById('notes-empty');

    // Gruppenauswahl nur zeigen, wenn die Anmerkung für Arbeitsgruppen sichtbar sein soll.
    function syncVisibility(form) {
        const groups = form.querySelector('[data-note-groups]');
        if (!groups) return;
        const checked = form.querySelector('input[name$="visibility"]:checked');
        groups.hidden = !checked || checked.value !== 'groups';
    }

    // "personal": eigene Anmerkungen ohne Arbeitsgruppen; Gruppen-ID: für diese Gruppe lesbar.
    function applyFilter() {
        const value = filter.value;
        const notes = list.querySelectorAll('.note');
        let shown = 0;
        notes.forEach(note => {
            const groups = note.dataset.groups ? note.dataset.groups.split(' ') : [];
            const visible = value === 'all'
                || (value === 'personal' ? groups.length === 0 : groups.includes(value));
            note.hidden = !visible;
            if (visible) shown++;
        });
        empty.hidden = shown > 0;
        empty.textContent = notes.length
            ? gettext('Keine Anmerkungen für diese Auswahl.')
            : gettext('Noch keine Anmerkungen vorhanden.');
    }

    function refresh() {
        root.querySelectorAll('form[data-note-action]').forEach(syncVisibility);
        applyFilter();
    }

    function showError(form, message) {
        const errors = form.querySelector('.note-errors');
        if (errors) errors.textContent = message;
        else alert(message);
    }

    root.addEventListener('submit', async event => {
        const form = event.target.closest('form[data-note-action]');
        if (!form) return;
        event.preventDefault();
        const action = form.dataset.noteAction;

        const button = form.querySelector('[type="submit"]');
        button.disabled = true;
        try {
            const response = await fetch(form.action, { method: 'POST', body: new FormData(form) });
            if (response.ok) {
                list.innerHTML = await response.text();
                if (action === 'create') {
                    form.reset();
                    showError(form, '');
                }
                refresh();
            } else {
                const data = await response.json().catch(() => ({}));
                showError(form, (data.errors || [gettext('Die Aktion konnte nicht ausgeführt werden.')]).join(' '));
            }
        } catch {
            showError(form, gettext('Verbindung zum Server fehlgeschlagen.'));
        } finally {
            button.disabled = false;
        }
    });

    // Bearbeiten: Text wird durch das Formular ersetzt; Abbrechen setzt es zurück.
    root.addEventListener('click', event => {
        const toggle = event.target.closest('[data-note-edit], [data-note-cancel]');
        if (!toggle) return;
        const note = toggle.closest('.note');
        const form = note.querySelector('form[data-note-action="edit"]');
        const editing = toggle.matches('[data-note-edit]');
        note.querySelector('.note-text').hidden = editing;
        note.querySelector('[data-note-edit]').hidden = editing;
        form.hidden = !editing;
        if (editing) {
            form.querySelector('textarea').focus();
        } else {
            form.reset();
            showError(form, '');
            syncVisibility(form);
        }
    });

    root.addEventListener('change', event => {
        if (event.target.matches('input[name$="visibility"]')) syncVisibility(event.target.form);
    });
    filter.addEventListener('change', applyFilter);

    refresh();
})();
