"""Configuration from environment variables and the .env file.

Small, dependency-free helpers used by settings.py and by the management command
`manage.py setup` (pdb_grapher/management/commands/setup.py), which writes the .env file.

Values from the environment take precedence over the .env file so that Docker, Portainer and
similar environments can override individual values.

Typical usage in settings.py:

    load_env_file(BASE_DIR / '.env')
    DEBUG = env_bool('DEBUG', False)
    ALLOWED_HOSTS = env_list('ALLOWED_HOSTS', ['localhost'])
    DATABASES = {'default': parse_database_url(env('DATABASE_URL', 'sqlite:///db.sqlite3'), BASE_DIR)}
"""
import os
from pathlib import Path
from urllib.parse import parse_qsl, unquote, urlsplit


def read_env_file(path):
    """Reads a .env file into a dict.

    Lines of the form `NAME=value` are parsed; empty lines, comments (#) and lines without `=`
    are skipped. Only the first `=` separates name and value, so values may contain `=` and `#`.
    Values may be enclosed in single or double quotes, which are removed.

    Args:
        path: Path of the .env file (str or Path).

    Returns:
        A dict mapping variable names to string values; empty if the file does not exist.
    """
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
    """Copies the values of a .env file into os.environ.

    Variables that are already set in the environment are not overwritten.

    Args:
        path: Path of the .env file (str or Path); a missing file is ignored.
    """
    for name, value in read_env_file(path).items():
        os.environ.setdefault(name, value)


def env(name, default=''):
    """Returns the value of an environment variable.

    An empty value counts as not set, because Docker Compose passes unset variables as empty
    strings.

    Args:
        name: Name of the environment variable.
        default: Value returned if the variable is unset or empty.

    Returns:
        The value of the variable as a string, or `default`.
    """
    value = os.environ.get(name)
    return default if value is None or value == '' else value


def env_bool(name, default=False):
    """Returns an environment variable interpreted as a boolean.

    Args:
        name: Name of the environment variable.
        default: Value returned if the variable is unset or empty.

    Returns:
        True for `1`, `true`, `yes` or `on` (case-insensitive), False for any other non-empty
        value, `default` if the variable is unset or empty.
    """
    value = os.environ.get(name)
    if value is None or value == '':
        return default
    return value.strip().lower() in ('1', 'true', 'yes', 'on')


def env_list(name, default=()):
    """Returns a comma-separated environment variable as a list.

    Example: `ALLOWED_HOSTS=bondamino.example.org,localhost`. Items are stripped and empty
    items are dropped.

    Args:
        name: Name of the environment variable.
        default: Iterable returned (as a list) if the variable is unset or blank.

    Returns:
        A list of strings.
    """
    value = os.environ.get(name)
    if value is None or value.strip() == '':
        return list(default)
    return [item.strip() for item in value.split(',') if item.strip()]


# URL scheme of DATABASE_URL -> Django database backend. To support another database, add its
# scheme here, its driver to a requirements file and, for the guided setup, an entry to
# DATABASES in pdb_grapher/management/commands/setup.py.
DATABASE_ENGINES = {
    'sqlite': 'django.db.backends.sqlite3',
    'postgres': 'django.db.backends.postgresql',
    'postgresql': 'django.db.backends.postgresql',
    'mysql': 'django.db.backends.mysql',
    'mariadb': 'django.db.backends.mysql',
}


def parse_database_url(url, base_dir):
    """Converts a DATABASE_URL into an entry for Django's DATABASES setting.

    Supported forms:
        sqlite:///db.sqlite3 (relative to base_dir), sqlite:////data/db.sqlite3 (absolute),
        postgres://user:password@host:5432/name?sslmode=require,
        mariadb://user:password@host:3306/name (also mysql://).

    Special characters in user and password must be URL-encoded (e.g. @ as %40). Further
    parameters after ? end up in OPTIONS; for MariaDB/MySQL the charset defaults to utf8mb4.

    Args:
        url: The database URL.
        base_dir: Directory that relative SQLite paths are resolved against.

    Returns:
        A dict with ENGINE and NAME (SQLite) or ENGINE, NAME, USER, PASSWORD, HOST, PORT and
        OPTIONS (server databases).

    Raises:
        ValueError: If the scheme is not supported or a SQLite URL has no file name.
    """
    parts = urlsplit(url)
    scheme = parts.scheme.lower()
    if scheme not in DATABASE_ENGINES:
        raise ValueError(f'Unsupported database type "{scheme}" in DATABASE_URL '
                         f'(supported: sqlite, postgres, mariadb, mysql)')
    engine = DATABASE_ENGINES[scheme]
    if scheme == 'sqlite':
        # sqlite:///relative -> path "/relative"; sqlite:////absolute -> path "//absolute"
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
    """Reads the secret key from a file.

    In Docker the file is created in the data volume by docker/entrypoint.sh.

    Args:
        path: Path of the file containing the key.

    Returns:
        The key without surrounding whitespace, or an empty string if the file cannot be read.
    """
    try:
        return Path(path).read_text(encoding='utf-8').strip()
    except OSError:
        return ''
