# Bondamino

Web tool to visualize bond lengths in amino acids as graph.

Bondamino is a Django application. It runs standalone on a server, in Docker, in Portainer
or on any platform that can run a container image. All variants are configured the same way:
through environment variables (standalone: the file `bondamino/.env`).

- [Requirements](#requirements)
- [Choosing a database](#choosing-a-database)
- [Installation](#installation)
  - [Standalone](#standalone)
  - [Docker Compose](#docker-compose)
  - [Portainer](#portainer)
  - [Other platforms](#other-platforms)
- [Configuration reference](#configuration-reference)
- [Updating](#updating)
- [Backup](#backup)
- [Development](#development)

## Requirements

- Standalone: Python 3.12 or newer, and for PostgreSQL or MariaDB a database server
- Docker / Portainer: Docker Engine 24 or newer with Compose v2
- A mail server (SMTP) for invitations and forgotten passwords. Without one, emails are only
  printed to the server console.

## Choosing a database

| Database | Use |
|---|---|
| PostgreSQL | Recommended for production |
| MariaDB / MySQL | Supported |
| SQLite | Testing, single users and small groups; no separate server needed |

MongoDB and DuckDB are not supported: Bondamino's data (users, groups, files, notes) is
relational, and DuckDB is designed for analytics, not for concurrent writes of a web application.

## Installation

### Standalone

1. Get the code and create a virtual environment:

   ```sh
   git clone https://github.com/PKittler/Bondamino.git
   cd Bondamino
   python3 -m venv .venv
   . .venv/bin/activate            # Windows: .venv\Scripts\activate
   ```

2. Install the packages for your database:

   ```sh
   pip install -r requirements.txt               # SQLite
   pip install -r requirements-postgresql.txt    # PostgreSQL
   pip install -r requirements-mariadb.txt       # MariaDB / MySQL
   ```

   `mysqlclient` needs the MariaDB/MySQL client library and a compiler, e.g. on Debian/Ubuntu
   `sudo apt install build-essential pkg-config default-libmysqlclient-dev`.

   For PostgreSQL or MariaDB, create an empty database and a user first, e.g. in PostgreSQL:

   ```sql
   CREATE USER bondamino WITH PASSWORD 'choose-a-password';
   CREATE DATABASE bondamino OWNER bondamino;
   ```

   or in MariaDB:

   ```sql
   CREATE DATABASE bondamino CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;
   CREATE USER 'bondamino'@'localhost' IDENTIFIED BY 'choose-a-password';
   GRANT ALL PRIVILEGES ON bondamino.* TO 'bondamino'@'localhost';
   ```

3. Run the guided setup:

   ```sh
   cd bondamino
   python manage.py setup
   ```

   It asks step by step for the database (and tests the connection), the domain, HTTPS and the
   mail server (and can send a test email). It then writes `bondamino/.env`, creates the
   database tables, collects the static files and creates the first administrator. Run it again
   at any time to change the configuration; the previous file is kept as `.env.bak`.

4. Start the web server:

   ```sh
   gunicorn bondamino.wsgi:application --bind 127.0.0.1:8000 --workers 3 --timeout 120
   ```

   On Windows, use waitress instead (`pip install waitress`):

   ```sh
   waitress-serve --listen=127.0.0.1:8000 bondamino.wsgi:application
   ```

5. Put a reverse proxy with HTTPS in front (nginx, Apache, Caddy, Traefik, …), answer "HTTPS"
   with yes in the setup and forward the original host and protocol. Example for nginx:

   ```nginx
   server {
       listen 443 ssl;
       server_name bondamino.example.org;
       # ssl_certificate ... / ssl_certificate_key ...
       client_max_body_size 50m;

       location / {
           proxy_pass http://127.0.0.1:8000;
           proxy_set_header Host $host;
           proxy_set_header X-Forwarded-Proto $scheme;
           proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
           proxy_read_timeout 120s;
       }
   }
   ```

6. To start Bondamino automatically, e.g. with systemd (`/etc/systemd/system/bondamino.service`):

   ```ini
   [Unit]
   Description=Bondamino
   After=network.target

   [Service]
   User=bondamino
   WorkingDirectory=/opt/Bondamino/bondamino
   ExecStart=/opt/Bondamino/.venv/bin/gunicorn bondamino.wsgi:application --bind 127.0.0.1:8000 --workers 3 --timeout 120
   Restart=on-failure

   [Install]
   WantedBy=multi-user.target
   ```

   then `sudo systemctl enable --now bondamino`.

### Docker Compose

There is one Compose file per database:

| File | Database |
|---|---|
| `docker-compose.yml` | PostgreSQL (recommended) |
| `docker-compose.mariadb.yml` | MariaDB |
| `docker-compose.sqlite.yml` | SQLite, single container |

1. Copy the Compose file of your choice and `.env.example` into a folder, or clone the
   repository. Rename `.env.example` to `.env` and fill in at least:

   - `DB_PASSWORD` (not needed for SQLite): letters and digits only, e.g. from
     `python -c "import secrets; print(secrets.token_hex(24))"`
   - `ALLOWED_HOSTS`: the domain under which Bondamino is reached
   - `ADMIN_USERNAME`, `ADMIN_EMAIL`, `ADMIN_PASSWORD`: the first administrator
   - the `EMAIL_*` settings of your mail server

   Leave `SECRET_KEY` empty: it is generated on the first start and stored in the data volume.
   `DATABASE_URL` is set by the Compose file; the entry in `.env` is ignored.

2. Start:

   ```sh
   docker compose up -d                                   # PostgreSQL
   docker compose -f docker-compose.mariadb.yml up -d     # MariaDB
   docker compose -f docker-compose.sqlite.yml up -d      # SQLite
   ```

   Bondamino is then reachable on port 8000 (change it with `BONDAMINO_PORT`). The image is
   pulled from `ghcr.io/pkittler/bondamino`; to build it yourself from the cloned repository,
   add `--build`.

3. For HTTPS, put a reverse proxy in front (see the nginx example above, or Traefik/Caddy in
   Docker) and set `HTTPS=1`.

Logs: `docker compose logs -f bondamino`. Management commands run inside the container, e.g.
`docker compose exec bondamino python manage.py createsuperuser`.

### Portainer

1. **Stacks → Add stack**, give it a name (e.g. `bondamino`).
2. Choose one of:
   - **Repository**: URL `https://github.com/PKittler/Bondamino`, reference `refs/heads/main`,
     Compose path `docker-compose.yml` (or `docker-compose.mariadb.yml`,
     `docker-compose.sqlite.yml`).
   - **Web editor**: paste the content of one of the Compose files and delete the line
     `build: .` (the Web editor cannot build images; the image is pulled instead).
3. Under **Environment variables**, add the variables from step 1 of
   [Docker Compose](#docker-compose) (`DB_PASSWORD`, `ALLOWED_HOSTS`, `ADMIN_*`, `EMAIL_*`, …).
   Portainer can also load them from a file ("Load variables from .env file").
4. **Deploy the stack.** Bondamino listens on port 8000 of the host.

### Other platforms

The image `ghcr.io/pkittler/bondamino` runs on any container platform (Kubernetes, Nomad,
Synology, Unraid, cloud services, …):

- Port: `8000` (HTTP)
- Volume: `/data` (uploaded PDB files, the generated secret key, and the SQLite database if used)
- Configuration: the environment variables in the [reference](#configuration-reference); for
  PostgreSQL or MariaDB set `DATABASE_URL`, without it SQLite in `/data` is used
- First administrator: `DJANGO_SUPERUSER_USERNAME`, `DJANGO_SUPERUSER_EMAIL`,
  `DJANGO_SUPERUSER_PASSWORD`
- On every start the container applies pending database migrations.

When mounting a host folder to `/data` instead of a named volume, it must be writable for
user ID 1000.

## Configuration reference

| Variable | Default | Meaning |
|---|---|---|
| `SECRET_KEY` | – (Docker: generated) | Long random value; signs sessions and password links. Keep it secret. |
| `SECRET_KEY_FILE` | Docker: `/data/secret_key` | File with the secret key, if `SECRET_KEY` is empty |
| `DEBUG` | `0` | `1` = development mode with detailed error pages. Never on a public server. |
| `ALLOWED_HOSTS` | `localhost` with `DEBUG=1`, otherwise none | Domains under which Bondamino is reached, comma-separated |
| `HTTPS` | `0` | `1` = served via HTTPS through a reverse proxy: secure cookies, `X-Forwarded-Proto` is trusted |
| `CSRF_TRUSTED_ORIGINS` | derived from `ALLOWED_HOSTS` | Only needed if the proxy changes the host name, e.g. `https://bondamino.example.org` |
| `DATABASE_URL` | `sqlite:///db.sqlite3` | `sqlite:///file` (relative to `manage.py`), `sqlite:////absolute/path`, `postgres://user:password@host:5432/name`, `mariadb://user:password@host:3306/name`. URL-encode special characters in user and password (`@` → `%40`). Options after `?`, e.g. `?sslmode=require`. |
| `CONN_MAX_AGE` | `60` | Seconds a database connection is reused |
| `MEDIA_ROOT` | `bondamino/media` (Docker: `/data/media`) | Uploaded PDB files |
| `STATIC_ROOT` | `bondamino/staticfiles` | Target of `collectstatic` |
| `TIME_ZONE` | `UTC` | e.g. `Europe/Berlin` |
| `EMAIL_BACKEND` | console | `django.core.mail.backends.smtp.EmailBackend` to send emails |
| `EMAIL_HOST`, `EMAIL_PORT` | `localhost`, `25` | SMTP server |
| `EMAIL_HOST_USER`, `EMAIL_HOST_PASSWORD` | – | SMTP login |
| `EMAIL_USE_TLS`, `EMAIL_USE_SSL` | `0` | `1` = STARTTLS (port 587) or SSL/TLS (port 465) |
| `DEFAULT_FROM_EMAIL` | `Bondamino <noreply@localhost>` | Sender of all emails |
| `LOG_LEVEL` | `WARNING` | Log level of the console output |
| `GUNICORN_WORKERS`, `GUNICORN_TIMEOUT` | `3`, `120` | Docker only: web server processes and request timeout |

Standalone, variables set in the environment take precedence over `bondamino/.env`. Another
file can be used with `BONDAMINO_ENV_FILE=/path/to/file`.

Password links (invitation, forgotten password) are valid for 24 hours and only once. Links in
emails use the domain under which the page was requested, so `ALLOWED_HOSTS` and the proxy's
`Host` header must contain the public domain.

## Updating

Standalone:

```sh
git pull
pip install -r requirements-postgresql.txt     # the file of your database
cd bondamino
python manage.py migrate
python manage.py collectstatic --noinput
sudo systemctl restart bondamino
```

Docker Compose: `docker compose pull && docker compose up -d` (with `-f` for MariaDB/SQLite).
Portainer: open the stack and choose **Pull and redeploy** (Repository: **Update the stack**
with "Re-pull image").

## Backup

Back up two things:

- the database: `pg_dump` (PostgreSQL), `mariadb-dump` (MariaDB), or the SQLite file
- the uploaded PDB files: `MEDIA_ROOT` (Docker: the volume `/data`, which also holds the
  secret key)

Docker example for PostgreSQL:

```sh
docker compose exec db pg_dump -U bondamino bondamino > bondamino.sql
docker run --rm -v bondamino_bondamino-data:/data -v "$PWD":/backup alpine tar czf /backup/bondamino-data.tgz -C /data .
```

## Development

```sh
pip install -r requirements.txt
cd bondamino
python manage.py setup        # SQLite, answer "Development mode" with yes
python manage.py runserver
python manage.py test pdb_grapher
```

With `DEBUG=1`, `ALLOWED_HOSTS` defaults to `localhost`.
