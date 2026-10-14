"""App configuration of pdb_grapher, the single Django app that makes up Bondamino."""
from django.apps import AppConfig


class PdbGrapherConfig(AppConfig):
    """Configuration of the pdb_grapher app."""
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'pdb_grapher'
    verbose_name = 'PDB Grapher'

    def ready(self):
        """Connects the signal receivers in ``signals.py``."""
        from . import signals  # noqa: F401  (imported for the side effect of registering the receivers)
