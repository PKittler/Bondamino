/**
 * @fileoverview Behaviour of the start page (template
 * `pdb_grapher/templates/pdb_grapher/index.html`, view `views.index`).
 *
 * Two independent parts, each in its own IIFE:
 *   - Upload dialog (`#upload-modal`): file selection by button or drag & drop,
 *     mutual exclusion of file and PDB ID, busy state while submitting, and
 *     reopening the dialog after a failed upload (`data-show`).
 *   - File list: row selection (mouse and keyboard) enabling the edit/delete
 *     toolbar for files the user uploaded (`data-edit-url`, `data-delete-url`),
 *     and instant client-side filtering by the search field.
 *
 * Globals: `bootstrap` (Bootstrap bundle) and `gettext`/`interpolate` (Django
 * JavaScript catalog, loaded in base.html).
 */
'use strict';

// ===== Upload dialog =====
(function () {
    const modal = document.getElementById('upload-modal');
    const input = document.getElementById(modal.dataset.fileInput);
    const pdbIdInput = document.getElementById(modal.dataset.pdbIdInput);
    const display = document.getElementById('file-name-display');
    const dropArea = document.getElementById('file-drop-area');
    const helpText = display.textContent;

    /**
     * Shows the chosen file name in the drop area and clears the PDB ID field.
     * @param {string} name File name to display.
     */
    function showFile(name) {
        display.textContent = name;
        dropArea.classList.add('has-file');
        pdbIdInput.value = '';
    }

    // File and PDB ID are mutually exclusive: the most recently chosen source wins.
    input.addEventListener('change', () => {
        if (input.files.length > 0) showFile(input.files[0].name);
    });
    pdbIdInput.addEventListener('input', () => {
        if (pdbIdInput.value.trim() && input.files.length > 0) {
            input.value = '';
            display.textContent = helpText;
            dropArea.classList.remove('has-file');
        }
    });

    ['dragover', 'dragenter'].forEach(evt => {
        dropArea.addEventListener(evt, e => { e.preventDefault(); dropArea.classList.add('drag-over'); });
    });
    ['dragleave', 'drop'].forEach(evt => {
        dropArea.addEventListener(evt, () => dropArea.classList.remove('drag-over'));
    });
    dropArea.addEventListener('drop', e => {
        e.preventDefault();
        if (e.dataTransfer.files.length > 0) {
            input.files = e.dataTransfer.files;
            showFile(e.dataTransfer.files[0].name);
        }
    });

    // Analysis or download from RCSB can take a few seconds: disable the button, show a spinner.
    input.form.addEventListener('submit', () => {
        const submit = document.getElementById('upload-submit');
        submit.disabled = true;
        submit.querySelector('.spinner-border').hidden = false;
    });

    // Failed upload: reopen the dialog showing the error messages.
    if ('show' in modal.dataset) bootstrap.Modal.getOrCreateInstance(modal).show();
})();

// ===== File list: selection and search =====
(function () {
    const rows = Array.from(document.querySelectorAll('.file-row'));
    const editLink = document.getElementById('edit-selected');
    const deleteForm = document.getElementById('delete-selected-form');
    const deleteButton = document.getElementById('delete-selected');
    const search = document.getElementById('file-search');
    const empty = document.getElementById('file-empty');
    let selected = null;

    /**
     * Selects a row and updates the toolbar. Edit and delete act only on the
     * selected row and only if the user uploaded the file (the row then carries
     * `data-edit-url` and `data-delete-url`).
     * @param {?HTMLTableRowElement} row Row to select, or null to clear the selection.
     */
    function select(row) {
        if (selected) {
            selected.classList.remove('table-active');
            selected.setAttribute('aria-selected', 'false');
        }
        selected = row;
        if (row) {
            row.classList.add('table-active');
            row.setAttribute('aria-selected', 'true');
        }
        const editUrl = row && row.dataset.editUrl;
        editLink.classList.toggle('disabled', !editUrl);
        editLink.setAttribute('aria-disabled', String(!editUrl));
        if (editUrl) editLink.href = editUrl; else editLink.removeAttribute('href');
        const deleteUrl = row && row.dataset.deleteUrl;
        deleteButton.disabled = !deleteUrl;
        if (deleteUrl) {
            deleteForm.action = deleteUrl;
            // The confirmation dialog is handled by the data-confirm handler in base.html.
            deleteForm.dataset.confirm = interpolate(gettext('„%s“ endgültig löschen?'), [row.dataset.fileName]);
        } else {
            deleteForm.removeAttribute('action');
            delete deleteForm.dataset.confirm;
        }
    }

    rows.forEach(row => {
        row.addEventListener('click', () => select(row));
        row.addEventListener('keydown', e => {
            if (e.target === row && (e.key === 'Enter' || e.key === ' ')) {
                e.preventDefault();
                select(row);
            }
        });
    });

    deleteForm.addEventListener('submit', e => {
        if (!selected || !selected.dataset.deleteUrl) e.preventDefault();
    });

    // The search filters the list instantly; on submit (Enter) it is kept in the URL (GET q).
    search.addEventListener('input', () => {
        const term = search.value.trim().toLowerCase();
        let shown = 0;
        rows.forEach(row => {
            row.hidden = !row.dataset.name.includes(term);
            if (!row.hidden) shown++;
        });
        if (selected && selected.hidden) select(null);
        empty.hidden = shown > 0;
        if (rows.length) empty.cells[0].textContent = gettext('Keine PDB-Dateien für diese Auswahl.');
    });
})();
