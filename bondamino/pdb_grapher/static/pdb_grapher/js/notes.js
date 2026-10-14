/**
 * @fileoverview Notes tab of the graph page (template
 * `pdb_grapher/templates/pdb_grapher/graph.html`, container `#notes`).
 *
 * Creating, editing, deleting and removing notes from a group run via AJAX:
 * every form with `data-note-action` is posted with fetch to its `action` URL
 * (views `note_create`, `note_edit`, `note_delete`, `note_remove_from_group`).
 * On success the server returns the re-rendered `_note_list.html`, which
 * replaces `#note-list`, so the graph page does not need to be reloaded; on
 * failure it returns JSON `{"errors": [...]}`.
 *
 * Also handles the edit/cancel toggles of each note, shows the group
 * checkboxes only for the visibility "groups", and filters the list by the
 * `#notes-filter` select (using each note's `data-groups`).
 *
 * Globals: `gettext` (Django JavaScript catalog, loaded in base.html).
 */
'use strict';

// ===== Notes tab =====
// Create, edit and delete run via AJAX; the server returns the updated list each
// time so that the graph page does not have to be reloaded.
(function () {
    const root = document.getElementById('notes');
    if (!root) return;
    const list = document.getElementById('note-list');
    const filter = document.getElementById('notes-filter');
    const empty = document.getElementById('notes-empty');

    /**
     * Shows the group selection only if the note is to be visible to work groups.
     * @param {!HTMLFormElement} form Note form containing the visibility fields.
     */
    function syncVisibility(form) {
        const groups = form.querySelector('[data-note-groups]');
        if (!groups) return;
        const checked = form.querySelector('input[name$="visibility"]:checked');
        groups.hidden = !checked || checked.value !== 'groups';
    }

    /**
     * Filters the notes by the value of `#notes-filter` and updates the empty
     * message. "all": every visible note; "personal": the user's own notes that
     * are not shared with any work group; a group id: notes readable by that group.
     */
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

    /** Re-applies visibility toggles and the filter after the list has changed. */
    function refresh() {
        root.querySelectorAll('form[data-note-action]').forEach(syncVisibility);
        applyFilter();
    }

    /**
     * Shows an error message in the form's `.note-errors` element, or as an
     * alert if the form has none (e.g. the delete and remove buttons).
     * @param {!HTMLFormElement} form Form the message belongs to.
     * @param {string} message Message to show; an empty string clears it.
     */
    function showError(form, message) {
        const errors = form.querySelector('.note-errors');
        if (errors) errors.textContent = message;
        else alert(message);
    }

    // All note forms (also those inside the replaced list) via event delegation.
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

    // Edit: the text is replaced by the form; Cancel resets the form and shows the text again.
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
