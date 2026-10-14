#!/usr/bin/env python
"""Django's command-line utility for administrative tasks.

Run it from this directory. Commands commonly used with Bondamino:

    python manage.py setup              Guided setup: writes .env and prepares the database
                                        (pdb_grapher/management/commands/setup.py).
    python manage.py runserver          Development server (with DEBUG=1).
    python manage.py migrate            Apply database migrations.
    python manage.py createcachetable   Create the cache table used for sign-in throttling.
    python manage.py test pdb_grapher   Run the test suite.
"""
import os
import sys


def main():
    """Runs the administrative task given on the command line.

    Raises:
        ImportError: If Django is not installed or the virtual environment is not active.
    """
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'bondamino.settings')
    try:
        from django.core.management import execute_from_command_line
    except ImportError as exc:
        raise ImportError(
            "Couldn't import Django. Are you sure it's installed and "
            "available on your PYTHONPATH environment variable? Did you "
            "forget to activate a virtual environment?"
        ) from exc
    execute_from_command_line(sys.argv)


if __name__ == '__main__':
    main()
