# Bondamino as a container image (PostgreSQL, MariaDB/MySQL and SQLite).
# Build: docker build -t bondamino .   Configuration via environment variables, see README.
#
# Two stages: "build" compiles the Python wheels (mysqlclient needs a compiler and the MariaDB
# headers), the final stage only installs the wheels and the MariaDB runtime library.
# On start, docker/entrypoint.sh runs migrations and then gunicorn on port 8000.

# --- Build stage: wheels for all dependencies, including both database drivers ---
FROM python:3.13-slim AS build
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential pkg-config default-libmysqlclient-dev \
    && rm -rf /var/lib/apt/lists/*
COPY requirements*.txt /tmp/
RUN pip wheel --no-cache-dir --wheel-dir /wheels \
    -r /tmp/requirements-postgresql.txt -r /tmp/requirements-mariadb.txt

# --- Runtime stage ---
FROM python:3.13-slim
# Defaults for the container: database, uploads (and the generated secret key) live in /data
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    DATABASE_URL=sqlite:////data/db.sqlite3 \
    MEDIA_ROOT=/data/media \
    STATIC_ROOT=/app/staticfiles
# Runtime library for mysqlclient; unprivileged user that owns the data volume
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
# Next to manage.py so that `docker compose exec bondamino python manage.py ...` works
WORKDIR /app/bondamino
# Collect static files at build time; a real secret key is not needed for that (dummy value)
RUN SECRET_KEY=collectstatic python manage.py collectstatic --noinput --verbosity 0 \
    && chown -R bondamino:bondamino /app/staticfiles

USER bondamino
VOLUME /data
EXPOSE 8000
ENTRYPOINT ["bondamino-entrypoint"]
