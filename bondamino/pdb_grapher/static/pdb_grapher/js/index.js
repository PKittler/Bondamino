'use strict';

// ===== Upload-Dialog =====
(function () {
    const modal = document.getElementById('upload-modal');
    const input = document.getElementById(modal.dataset.fileInput);
    const pdbIdInput = document.getElementById(modal.dataset.pdbIdInput);
    const display = document.getElementById('file-name-display');
    const dropArea = document.getElementById('file-drop-area');
    const helpText = display.textContent;

    function showFile(name) {
        display.textContent = name;
        dropArea.classList.add('has-file');
        pdbIdInput.value = '';
    }

    // Datei und PDB-ID schließen sich aus: Die zuletzt gewählte Quelle gilt.
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

    // Analyse bzw. Abruf bei RCSB kann einige Sekunden dauern.
    input.form.addEventListener('submit', () => {
        const submit = document.getElementById('upload-submit');
        submit.disabled = true;
        submit.querySelector('.spinner-border').hidden = false;
    });

    // Fehlgeschlagener Upload: Dialog mit den Fehlermeldungen wieder öffnen.
    if ('show' in modal.dataset) bootstrap.Modal.getOrCreateInstance(modal).show();
})();

// ===== Dateiliste: Auswahl und Suche =====
(function () {
    const rows = Array.from(document.querySelectorAll('.file-row'));
    const editLink = document.getElementById('edit-selected');
    const deleteForm = document.getElementById('delete-selected-form');
    const deleteButton = document.getElementById('delete-selected');
    const search = document.getElementById('file-search');
    const empty = document.getElementById('file-empty');
    let selected = null;

    // Bearbeiten und Löschen wirken nur auf die ausgewählte Zeile und nur,
    // wenn der User die Datei hochgeladen hat.
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
            // Rückfrage übernimmt der data-confirm-Handler in base.html.
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

    // Suche filtert die Liste sofort; beim Absenden (Enter) bleibt sie in der URL erhalten.
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
