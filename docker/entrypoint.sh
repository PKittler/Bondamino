#!/bin/sh
# Start im Container: Secret Key bereitstellen, Datenbank migrieren, ersten Administrator
# anlegen (optional), dann den Webserver starten.
set -e
cd /app/bondamino

mkdir -p "$(dirname "${SECRET_KEY_FILE:-/data/secret_key}")" "${MEDIA_ROOT:-/data/media}"

# Ohne SECRET_KEY: Schlüssel einmalig erzeugen und im Datenvolume speichern
if [ -z "$SECRET_KEY" ]; then
    export SECRET_KEY_FILE="${SECRET_KEY_FILE:-/data/secret_key}"
    if [ ! -s "$SECRET_KEY_FILE" ]; then
        python -c "import secrets; print(secrets.token_urlsafe(50))" > "$SECRET_KEY_FILE"
        chmod 600 "$SECRET_KEY_FILE"
        echo "Generated a new secret key in $SECRET_KEY_FILE"
    fi
fi

# Die Datenbank (eigener Container) braucht beim ersten Start etwas länger
attempt=1
until python manage.py migrate --noinput; do
    if [ "$attempt" -ge 30 ]; then
        echo "Database not reachable, giving up." >&2
        exit 1
    fi
    echo "Waiting for the database ($attempt/30) ..."
    attempt=$((attempt + 1))
    sleep 2
done
python manage.py createcachetable

# Erster Administrator aus DJANGO_SUPERUSER_USERNAME/_EMAIL/_PASSWORD, falls noch keiner existiert
if [ -n "$DJANGO_SUPERUSER_USERNAME" ] && [ -n "$DJANGO_SUPERUSER_PASSWORD" ]; then
    python manage.py shell -c "
from django.contrib.auth import get_user_model
import sys
sys.exit(0 if get_user_model().objects.filter(is_superuser=True).exists() else 1)
" || python manage.py createsuperuser --noinput
fi

if [ "$#" -gt 0 ]; then
    exec "$@"
fi
exec gunicorn bondamino.wsgi:application \
    --bind 0.0.0.0:8000 \
    --workers "${GUNICORN_WORKERS:-3}" \
    --timeout "${GUNICORN_TIMEOUT:-120}" \
    --access-logfile - \
    --no-control-socket
