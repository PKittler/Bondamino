"""WSGI entry point of the Bondamino project.

Exposes the WSGI callable as the module-level variable `application`. This is what the
production servers load: gunicorn in the Docker image (`bondamino.wsgi:application`, see
docker/entrypoint.sh) and waitress on Windows (see README).

See https://docs.djangoproject.com/en/5.2/howto/deployment/wsgi/ for details.
"""

import os

from django.core.wsgi import get_wsgi_application

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'bondamino.settings')

application = get_wsgi_application()
