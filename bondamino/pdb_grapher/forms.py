from django import forms

class PDBUploadForm(forms.Form):
    pdb_file = forms.FileField(
        label='Select PDB file (.pdb, .ent)',
        widget=forms.ClearableFileInput(attrs={'accept': '.pdb,.ent'}),
        required=True
    )

    def clean_pdb_file(self):
        file = self.cleaned_data.get('pdb_file')
        if file:
            filename = file.name
            if not (filename.lower().endswith('.pdb') or filename.lower().endswith('.ent')):
                raise forms.ValidationError("Only .pdb or .ent files are allowed.")
        return file


class GraphOptionsForm(forms.Form):
    model_id = forms.ChoiceField(
        label='Choose model',
        widget=forms.Select,
        required=True
    )
    label_type = forms.ChoiceField(
        label='Label Nodes',
        choices=[('atom_name', 'atom name (e.g. CA, CB)'), ('element', 'element symbol (e.g. C, N, O)')],
        widget=forms.RadioSelect,
        initial='atom_name',
        required=True
    )

    def __init__(self, *args, **kwargs):
        num_models = kwargs.pop('num_models', 1)
        super().__init__(*args, **kwargs)
        self.fields['model_id'].choices = [(i, f"Model {i}") for i in range(num_models)]