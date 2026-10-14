"""ASGI entry point of the Bondamino project.

Exposes the ASGI callable as the module-level variable `application`. Bondamino is deployed
via WSGI (see wsgi.py); this module is kept for ASGI servers such as uvicorn or daphne.

See https://docs.djangoproject.com/en/5.2/howto/deployment/asgi/ for details.
"""

import os

from django.core.asgi import get_asgi_application

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'bondamino.settings')

application = get_asgi_application()
