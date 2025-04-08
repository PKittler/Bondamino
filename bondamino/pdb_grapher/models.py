import uuid
from django.db import models
from django.utils import timezone

def pdb_upload_path(instance, filename):
    # pdb file will be uploaded to MEDIA_ROOT/uploads/<uuid>.<ext>

    ext = filename.split('.')[-1]
    return f'uploads/{instance.id}.{ext}'

def graph_upload_path(instance, filename):
    # graph file will be uploaded to MEDIA_ROOT/generated_graphs/<uuid>.html
    return f'generated_graphs/{instance.id}.html'

class PDBFile(models.Model):

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    original_filename = models.CharField(max_length=255)
    pdb_file = models.FileField(upload_to=pdb_upload_path)
    uploaded_at = models.DateTimeField(default=timezone.now)
    num_models = models.IntegerField(null=True, blank=True)

    def __str__(self):
        return self.original_filename

    def get_id(self):
        return self.id

    def get_original_filename(self):
        return self.original_filename

class GeneratedGraph(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    pdb_file = models.ForeignKey(PDBFile, on_delete=models.CASCADE, related_name='graphs')
    model_id_used = models.IntegerField()
    label_type_used = models.CharField(max_length=20) # 'atom_name' or 'element'
    graph_html_file = models.FileField(upload_to=graph_upload_path)
    created_at = models.DateTimeField(default=timezone.now)

    def __str__(self):
        return f"Graph for {self.pdb_file.original_filename} (Model {self.model_id_used})"