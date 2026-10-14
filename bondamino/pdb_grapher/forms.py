from django import forms
from django.core.files.uploadedfile import SimpleUploadedFile
from django.utils.translation import gettext_lazy as _

from . import rcsb
from .models import Language, Note, PDBFile, UserSettings, WorkGroup


class UserSettingsForm(forms.ModelForm):
    """Einstellungen auf der Account-Seite; zur Wahl stehen die im Admin angelegten Sprachen."""
    language = forms.ChoiceField(label=_('Sprache'), widget=forms.Select(attrs={'class': 'form-select'}))

    class Meta:
        model = UserSettings
        fields = ['language', 'theme']
        widgets = {'theme': forms.Select(attrs={'class': 'form-select'})}

    def __init__(self, *args, current_language, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['language'].choices = [(lang.code, lang.name) for lang in Language.objects.all()]
        if not self.instance.language:
            self.initial['language'] = current_language


class AvailabilityForm(forms.ModelForm):
    """Verfügbarkeit einer PDB-Datei: persönlich und/oder für Arbeitsgruppen.

    Auswählbar sind die Gruppen des Users sowie Gruppen, in denen die Datei bereits
    freigegeben ist (damit der Uploader sie auch nach einem Gruppenwechsel entfernen kann).
    """
    personal = forms.BooleanField(
        label=_('Für mich persönlich'), required=False, initial=True,
        widget=forms.CheckboxInput(attrs={'class': 'form-check-input'}),
    )
    groups = forms.ModelMultipleChoiceField(
        label=_('Arbeitsgruppen'), queryset=WorkGroup.objects.none(), required=False,
        widget=forms.CheckboxSelectMultiple(attrs={'class': 'form-check-input'}),
    )

    class Meta:
        model = PDBFile
        fields = ['personal', 'groups']

    def __init__(self, *args, user, **kwargs):
        super().__init__(*args, **kwargs)
        groups = WorkGroup.objects.filter(members=user)
        if self.instance.pk:
            groups = groups | self.instance.groups.all()
        self.fields['groups'].queryset = groups.distinct()

    def clean(self):
        cleaned = super().clean()
        if not cleaned.get('personal') and not cleaned.get('groups'):
            raise forms.ValidationError(_(
                'Die Datei muss für dich persönlich oder für mindestens eine Arbeitsgruppe '
                'verfügbar sein. Um sie zu entfernen, lösche sie.'
            ))
        return cleaned


class PDBUploadForm(AvailabilityForm):
    """Neue PDB-Datei: entweder hochgeladen oder per PDB-ID von RCSB geladen."""
    pdb_file = forms.FileField(
        label=_('PDB-Datei'), required=False,
        help_text=_('Unterstützt .pdb Dateien'),
    )
    pdb_id = forms.CharField(
        label=_('PDB-ID'), required=False, max_length=10,
        widget=forms.TextInput(attrs={
            'class': 'form-control font-monospace', 'placeholder': _('z. B. 1LYZ'),
            'autocomplete': 'off', 'spellcheck': 'false',
        }),
    )

    class Meta(AvailabilityForm.Meta):
        fields = ['pdb_file', 'pdb_id', 'personal', 'groups']

    def clean_pdb_file(self):
        f = self.cleaned_data['pdb_file']
        if not f:
            return f
        if not f.name.lower().endswith('.pdb'):
            raise forms.ValidationError(_('Nur .pdb Dateien sind erlaubt.'))
        if f.size > rcsb.MAX_SIZE:
            raise forms.ValidationError(_('Datei darf maximal 20 MB groß sein.'))
        return f

    def clean_pdb_id(self):
        pdb_id = self.cleaned_data['pdb_id']
        if not pdb_id:
            return pdb_id
        try:
            return rcsb.normalize_pdb_id(pdb_id)
        except rcsb.RCSBError as e:
            raise forms.ValidationError(str(e))

    def clean(self):
        cleaned = super().clean()
        if self.errors:
            return cleaned  # Kein Abruf bei RCSB, solange das Formular ohnehin ungültig ist
        pdb_file, pdb_id = cleaned.get('pdb_file'), cleaned.get('pdb_id')
        if pdb_file and pdb_id:
            self.add_error('pdb_file', _('Bitte entweder eine Datei auswählen oder eine PDB-ID eingeben, nicht beides.'))
        elif not pdb_file and not pdb_id:
            self.add_error('pdb_file', _('Bitte eine Datei auswählen oder eine PDB-ID eingeben.'))
        elif pdb_id:
            try:
                content = rcsb.fetch_pdb(pdb_id)
            except rcsb.RCSBError as e:
                self.add_error('pdb_id', str(e))
            else:
                cleaned['pdb_file'] = SimpleUploadedFile(f'{pdb_id}.pdb', content)
        return cleaned


class NoteForm(forms.ModelForm):
    """Anmerkung mit Sichtbarkeit: nur für den Verfasser oder zusätzlich für ausgewählte
    Arbeitsgruppen.

    Auswählbar sind die Gruppen des Users sowie Gruppen, für die die Anmerkung bereits
    lesbar ist.
    """
    PERSONAL, GROUPS = 'personal', 'groups'

    visibility = forms.ChoiceField(
        label=_('Sichtbar für'), initial=PERSONAL,
        choices=[(PERSONAL, _('Nur für mich')), (GROUPS, _('Für Arbeitsgruppen'))],
        widget=forms.RadioSelect(attrs={'class': 'form-check-input'}),
    )
    groups = forms.ModelMultipleChoiceField(
        label=_('Arbeitsgruppen'), queryset=WorkGroup.objects.none(), required=False,
        widget=forms.CheckboxSelectMultiple(attrs={'class': 'form-check-input'}),
    )

    class Meta:
        model = Note
        fields = ['text', 'groups']
        widgets = {
            'text': forms.Textarea(attrs={
                'class': 'form-control', 'rows': 3,
                'placeholder': _('Anmerkung schreiben …'),
            }),
        }

    def __init__(self, *args, user, **kwargs):
        super().__init__(*args, **kwargs)
        groups = WorkGroup.objects.filter(members=user)
        if self.instance.pk:
            current = self.instance.groups.all()
            groups = groups | current
            self.initial['visibility'] = self.GROUPS if current else self.PERSONAL
        self.fields['groups'].queryset = groups.distinct()

    def clean(self):
        cleaned = super().clean()
        if cleaned.get('visibility') == self.GROUPS:
            if not cleaned.get('groups') and 'groups' not in self.errors:
                raise forms.ValidationError(
                    _('Wähle mindestens eine Arbeitsgruppe aus, die die Anmerkung lesen darf.')
                )
        else:
            cleaned['groups'] = WorkGroup.objects.none()
        return cleaned
