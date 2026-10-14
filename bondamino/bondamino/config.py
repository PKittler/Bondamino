"""Konfiguration aus Umgebungsvariablen und der Datei .env (für settings.py und `manage.py setup`).

Werte aus der Umgebung haben Vorrang vor der .env-Datei, damit Docker, Portainer und ähnliche
Umgebungen einzelne Werte überschreiben können.
"""
import os
from pathlib import Path
from urllib.parse import parse_qsl, unquote, urlsplit


def read_env_file(path):
    """Zeilen `NAME=wert` als dict; leere Zeilen und Kommentare (#) werden übersprungen.
    Werte dürfen in einfache oder doppelte Anführungszeichen gesetzt sein."""
    values = {}
    path = Path(path)
    if not path.is_file():
        return values
    for line in path.read_text(encoding='utf-8').splitlines():
        line = line.strip()
        if not line or line.startswith('#') or '=' not in line:
            continue
        name, value = line.split('=', 1)
        value = value.strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in '"\'':
            value = value[1:-1]
        values[name.strip()] = value
    return values


def load_env_file(path):
    """Werte aus der .env-Datei in os.environ übernehmen, ohne gesetzte Variablen zu überschreiben."""
    for name, value in read_env_file(path).items():
        os.environ.setdefault(name, value)


def env(name, default=''):
    """Wert einer Variable; leer gilt als nicht gesetzt (Docker Compose setzt leere Werte)."""
    value = os.environ.get(name)
    return default if value is None or value == '' else value


def env_bool(name, default=False):
    value = os.environ.get(name)
    if value is None or value == '':
        return default
    return value.strip().lower() in ('1', 'true', 'yes', 'on')


def env_list(name, default=()):
    """Kommagetrennte Liste, z. B. ALLOWED_HOSTS=bondamino.example.org,localhost."""
    value = os.environ.get(name)
    if value is None or value.strip() == '':
        return list(default)
    return [item.strip() for item in value.split(',') if item.strip()]


DATABASE_ENGINES = {
    'sqlite': 'django.db.backends.sqlite3',
    'postgres': 'django.db.backends.postgresql',
    'postgresql': 'django.db.backends.postgresql',
    'mysql': 'django.db.backends.mysql',
    'mariadb': 'django.db.backends.mysql',
}


def parse_database_url(url, base_dir):
    """DATABASE_URL als Eintrag für DATABASES.

    sqlite:///db.sqlite3 (relativ zu base_dir), sqlite:////data/db.sqlite3 (absolut),
    postgres://user:passwort@host:5432/name?sslmode=require, mariadb://user:passwort@host:3306/name.
    Sonderzeichen in Benutzer und Passwort müssen URL-kodiert sein (z. B. @ als %40).
    Weitere Parameter nach ? landen in OPTIONS.
    """
    parts = urlsplit(url)
    scheme = parts.scheme.lower()
    if scheme not in DATABASE_ENGINES:
        raise ValueError(f'Unsupported database type "{scheme}" in DATABASE_URL '
                         f'(supported: sqlite, postgres, mariadb, mysql)')
    engine = DATABASE_ENGINES[scheme]
    if scheme == 'sqlite':
        # sqlite:///relativ → Pfad "/relativ"; sqlite:////absolut → Pfad "//absolut"
        path = unquote(parts.path)
        path = path[1:] if path.startswith('/') else path
        if not path:
            raise ValueError('DATABASE_URL for SQLite needs a file name, e.g. sqlite:///db.sqlite3')
        name = Path(path) if Path(path).is_absolute() or path.startswith('/') else Path(base_dir) / path
        return {'ENGINE': engine, 'NAME': str(name)}

    options = dict(parse_qsl(parts.query))
    if engine == 'django.db.backends.mysql':
        options.setdefault('charset', 'utf8mb4')
    return {
        'ENGINE': engine,
        'NAME': unquote(parts.path.lstrip('/')),
        'USER': unquote(parts.username or ''),
        'PASSWORD': unquote(parts.password or ''),
        'HOST': parts.hostname or '',
        'PORT': str(parts.port or ''),
        'OPTIONS': options,
    }


def read_secret_key(path):
    """Secret Key aus einer Datei (Docker: im Datenvolume, siehe docker-entrypoint.sh)."""
    try:
        return Path(path).read_text(encoding='utf-8').strip()
    except OSError:
        return ''
