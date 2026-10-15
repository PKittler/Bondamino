"""Django settings for the Bondamino project.

All deployment-specific values come from environment variables, so the same code runs
standalone, in Docker/Portainer or on other platforms without editing this file. When
running standalone, a `.env` file next to `manage.py` (written by `python manage.py setup`,
template: `.env.example` in the repository root) is loaded first; variables that are already
set in the environment take precedence. The helpers that read and parse the variables live in
`bondamino/config.py`.

Environment variables (defaults and meaning are documented in INSTALL.md, section
"Configuration reference"):

    BONDAMINO_ENV_FILE: Path of the .env file to load instead of `bondamino/.env`.
    SECRET_KEY, SECRET_KEY_FILE: Secret key, or a file containing it (Docker).
    DEBUG: Development mode with detailed error pages.
    ALLOWED_HOSTS: Comma-separated host names under which Bondamino is reached.
    HTTPS: Running behind an HTTPS reverse proxy (secure cookies, X-Forwarded-*).
    CSRF_TRUSTED_ORIGINS: Allowed form origins; derived from ALLOWED_HOSTS if empty.
    DATABASE_URL, CONN_MAX_AGE: Database connection (see config.parse_database_url).
    MEDIA_ROOT, STATIC_ROOT: Directories for uploaded PDB files and collected static files.
    TIME_ZONE: Default time zone for display (users can choose their own in their settings).
    EMAIL_BACKEND, EMAIL_HOST, EMAIL_PORT, EMAIL_HOST_USER, EMAIL_HOST_PASSWORD,
        EMAIL_USE_TLS, EMAIL_USE_SSL, DEFAULT_FROM_EMAIL: Sending of invitation and
        password reset emails.
    LOG_LEVEL: Level of the console log output.

To add a new setting, read it with `env()`, `env_bool()` or `env_list()` from config.py,
give it a safe default, and document it in INSTALL.md ("Configuration reference"), in
`.env.example` and, if it is needed in containers, in the docker-compose files.

See https://docs.djangoproject.com/en/5.2/topics/settings/ and
https://docs.djangoproject.com/en/5.2/ref/settings/ for Django's own settings.
"""

from pathlib import Path
import os

from .config import env, env_bool, env_list, load_env_file, parse_database_url, read_secret_key

# Build paths inside the project like this: BASE_DIR / 'subdir'.
BASE_DIR = Path(__file__).resolve().parent.parent

# Configuration comes from environment variables, supplemented locally by the .env file next to
# manage.py (create it with `python manage.py setup`, template: .env.example). Variables that are
# already set take precedence.
load_env_file(env('BONDAMINO_ENV_FILE', str(BASE_DIR / '.env')))

# Never in the code: taken from SECRET_KEY or from the file SECRET_KEY_FILE (the Docker
# entrypoint creates it in the data volume). Django refuses to start without a key.
SECRET_KEY = env('SECRET_KEY') or (read_secret_key(env('SECRET_KEY_FILE')) if env('SECRET_KEY_FILE') else '')

# Enable for development only (DEBUG=1)
DEBUG = env_bool('DEBUG', False)

# Domains under which Bondamino is reached, e.g. ALLOWED_HOSTS=bondamino.example.org;
# with DEBUG the local addresses are allowed by default.
ALLOWED_HOSTS = env_list('ALLOWED_HOSTS', ['localhost', '127.0.0.1', '[::1]'] if DEBUG else [])

# Bondamino runs behind a reverse proxy with HTTPS (HTTPS=1): secure cookies and trust the
# X-Forwarded-Proto header set by the proxy. pdb_grapher/throttle.py additionally takes the
# client address from X-Forwarded-For when HTTPS is set.
HTTPS = env_bool('HTTPS', False)
SESSION_COOKIE_SECURE = HTTPS
CSRF_COOKIE_SECURE = HTTPS
if HTTPS:
    SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTO', 'https')
# Forms may be submitted from these origins; default: the ALLOWED_HOSTS (without wildcards)
CSRF_TRUSTED_ORIGINS = env_list('CSRF_TRUSTED_ORIGINS', [
    f"{'https' if HTTPS else 'http'}://{host.lstrip('.')}" for host in ALLOWED_HOSTS
    if host not in ('*', '[::1]') and not host.startswith('.')
])


# Application definition

INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'pdb_grapher',
]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    # Serve static files without a separate web server
    'whitenoise.middleware.WhiteNoiseMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    # Instead of django.middleware.locale.LocaleMiddleware: only the languages created in the
    # admin; placed after AuthenticationMiddleware because it reads the user's language setting
    'pdb_grapher.middleware.LanguageMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = 'bondamino.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
                'django.template.context_processors.i18n',
                'pdb_grapher.context_processors.languages',
            ],
        },
    },
]

WSGI_APPLICATION = 'bondamino.wsgi.application'


# Database
# https://docs.djangoproject.com/en/5.2/ref/settings/#databases

# DATABASE_URL, e.g. postgres://user:password@host:5432/bondamino or
# mariadb://user:password@host:3306/bondamino; if unset, SQLite next to manage.py.
# Parsed by config.parse_database_url().
DATABASES = {
    'default': {
        **parse_database_url(env('DATABASE_URL', 'sqlite:///db.sqlite3'), BASE_DIR),
        'CONN_MAX_AGE': int(env('CONN_MAX_AGE', '60')),
    }
}

# Cache in the database so that all server processes see the same counters for limiting sign-in
# and password reset attempts (pdb_grapher/throttle.py). Create the table with
# `manage.py createcachetable` (the Docker entrypoint does this on every start).
CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.db.DatabaseCache',
        'LOCATION': 'bondamino_cache',
    }
}


# Password validation
# https://docs.djangoproject.com/en/5.2/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]


# Internationalization
# https://docs.djangoproject.com/en/5.2/topics/i18n/

# Default language and source language of the texts (translation ids are German). Selectable
# are the languages created in the admin (pdb_grapher.models.Language); LANGUAGES stays Django's
# complete list so that any language with a translation can be created there.
LANGUAGE_CODE = 'de'

TIME_ZONE = env('TIME_ZONE', 'UTC')

USE_I18N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/5.2/howto/static-files/

STATIC_URL = '/static/'
# Target of `manage.py collectstatic`; WhiteNoise serves the files from there
STATIC_ROOT = env('STATIC_ROOT', str(BASE_DIR / 'staticfiles'))
STORAGES = {
    'default': {'BACKEND': 'django.core.files.storage.FileSystemStorage'},
    'staticfiles': {'BACKEND': 'whitenoise.storage.CompressedStaticFilesStorage'},
}

# Default primary key field type
# https://docs.djangoproject.com/en/5.2/ref/settings/#default-auto-field

DEFAULT_AUTO_FIELD = 'django.db.models.BigAutoField'

# Uploaded PDB files are stored here (not served publicly, see bondamino/urls.py)
MEDIA_URL = '/media/'
MEDIA_ROOT = env('MEDIA_ROOT', os.path.join(BASE_DIR, 'media'))

# Login: users are created by the admin only, there is no self-registration
LOGIN_URL = 'login'
LOGIN_REDIRECT_URL = 'pdb_grapher:index'
LOGOUT_REDIRECT_URL = 'login'

# Links for setting the password (invitation, forgotten password) are valid for 24 hours
PASSWORD_RESET_TIMEOUT = 24 * 60 * 60

# Sending emails (invitation, forgotten password). Without configuration the emails are printed
# to the server console; set the environment variables to send them via SMTP.
EMAIL_BACKEND = env('EMAIL_BACKEND', 'django.core.mail.backends.console.EmailBackend')
EMAIL_HOST = env('EMAIL_HOST', 'localhost')
EMAIL_PORT = int(env('EMAIL_PORT', '25'))
EMAIL_HOST_USER = env('EMAIL_HOST_USER')
EMAIL_HOST_PASSWORD = env('EMAIL_HOST_PASSWORD')
EMAIL_USE_TLS = env_bool('EMAIL_USE_TLS')
EMAIL_USE_SSL = env_bool('EMAIL_USE_SSL')
DEFAULT_FROM_EMAIL = env('DEFAULT_FROM_EMAIL', 'Bondamino <noreply@localhost>')
SERVER_EMAIL = DEFAULT_FROM_EMAIL

# Errors and warnings go to the console (Docker: docker logs), also without DEBUG
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'handlers': {'console': {'class': 'logging.StreamHandler'}, 'null': {'class': 'logging.NullHandler'}},
    'root': {'handlers': ['console'], 'level': env('LOG_LEVEL', 'WARNING')},
    # Django answers requests with an unknown host (mostly bots) with 400; do not log them
    'loggers': {'django.security.DisallowedHost': {'handlers': ['null'], 'propagate': False}},
}
