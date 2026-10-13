from django.apps import AppConfig


class PdbGrapherConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'pdb_grapher'
    verbose_name = 'PDB Grapher'

    def ready(self):
        from . import signals  # noqa: F401
