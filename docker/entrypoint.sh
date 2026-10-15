#!/bin/sh
# Container start: provide the secret key, migrate the database, create the cache table,
# create the first administrator (optional), then start the web server (gunicorn).
# Arguments passed to the container (e.g. `docker run ... python manage.py shell`) are run
# instead of gunicorn after the preparation steps.
set -e
cd /app/bondamino

mkdir -p "$(dirname "${SECRET_KEY_FILE:-/data/secret_key}")" "${MEDIA_ROOT:-/data/media}"

# Without SECRET_KEY: generate a key once and store it in the data volume (read by settings.py
# via SECRET_KEY_FILE)
if [ -z "$SECRET_KEY" ]; then
    export SECRET_KEY_FILE="${SECRET_KEY_FILE:-/data/secret_key}"
    if [ ! -s "$SECRET_KEY_FILE" ]; then
        python -c "import secrets; print(secrets.token_urlsafe(50))" > "$SECRET_KEY_FILE"
        chmod 600 "$SECRET_KEY_FILE"
        echo "Generated a new secret key in $SECRET_KEY_FILE"
    fi
fi

# The database (separate container) takes a while on the first start: retry for up to ~60 s
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
# Cache table for the sign-in throttling counters (no-op if it already exists)
python manage.py createcachetable

# First administrator from DJANGO_SUPERUSER_USERNAME/_EMAIL/_PASSWORD, if no superuser exists yet.
# Problems only produce a warning: a failing start would make the container restart endlessly.
if [ -n "$DJANGO_SUPERUSER_USERNAME" ] && [ -n "$DJANGO_SUPERUSER_PASSWORD" ]; then
    python manage.py shell <<'PY' || echo "Warning: the first administrator could not be created." >&2
import os
from django.contrib.auth import get_user_model
from django.contrib.auth.password_validation import validate_password
from django.core.exceptions import ValidationError
from django.utils import translation

translation.activate('en')
User = get_user_model()
name = os.environ['DJANGO_SUPERUSER_USERNAME']
email = os.environ.get('DJANGO_SUPERUSER_EMAIL', '')
password = os.environ['DJANGO_SUPERUSER_PASSWORD']
if User.objects.filter(is_superuser=True).exists():
    pass
elif User.objects.filter(username=name).exists():
    print(f'Warning: the user "{name}" already exists but is no administrator; ADMIN_USERNAME is ignored. '
          'Give the user administrator rights or choose another ADMIN_USERNAME.')
else:
    try:
        validate_password(password, User(username=name, email=email))
    except ValidationError as error:
        print('Warning: ADMIN_PASSWORD is weak: ' + ' '.join(error.messages))
    User.objects.create_superuser(name, email, password)
    print('Superuser created successfully.')
PY
fi

# Run a command given as container arguments instead of the web server
if [ "$#" -gt 0 ]; then
    exec "$@"
fi
exec gunicorn bondamino.wsgi:application \
    --bind 0.0.0.0:8000 \
    --workers "${GUNICORN_WORKERS:-3}" \
    --timeout "${GUNICORN_TIMEOUT:-120}" \
    --access-logfile - \
    --no-control-socket
