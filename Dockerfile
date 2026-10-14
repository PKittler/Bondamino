# Bondamino als Container-Image (PostgreSQL, MariaDB/MySQL und SQLite).
# Bauen: docker build -t bondamino .   Konfiguration über Umgebungsvariablen, siehe README.

FROM python:3.13-slim AS build
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential pkg-config default-libmysqlclient-dev \
    && rm -rf /var/lib/apt/lists/*
COPY requirements*.txt /tmp/
RUN pip wheel --no-cache-dir --wheel-dir /wheels \
    -r /tmp/requirements-postgresql.txt -r /tmp/requirements-mariadb.txt

FROM python:3.13-slim
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    DATABASE_URL=sqlite:////data/db.sqlite3 \
    MEDIA_ROOT=/data/media \
    STATIC_ROOT=/app/staticfiles
RUN apt-get update \
    && apt-get install -y --no-install-recommends libmariadb3 \
    && rm -rf /var/lib/apt/lists/* \
    && useradd --system --uid 1000 --home-dir /app bondamino \
    && mkdir -p /data /app \
    && chown bondamino:bondamino /data
COPY --from=build /wheels /wheels
RUN pip install --no-cache-dir /wheels/* && rm -rf /wheels

COPY --chown=bondamino:bondamino bondamino /app/bondamino
COPY --chmod=755 docker/entrypoint.sh /usr/local/bin/bondamino-entrypoint
# Neben manage.py, damit `docker compose exec bondamino python manage.py ...` funktioniert
WORKDIR /app/bondamino
# Statische Dateien beim Bauen sammeln; der Schlüssel wird dafür nicht gebraucht
RUN SECRET_KEY=collectstatic python manage.py collectstatic --noinput --verbosity 0 \
    && chown -R bondamino:bondamino /app/staticfiles

USER bondamino
VOLUME /data
EXPOSE 8000
ENTRYPOINT ["bondamino-entrypoint"]
