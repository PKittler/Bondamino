"""Geführte Einrichtung: fragt Datenbank, Domain und E-Mail ab, schreibt die Datei .env und richtet
danach die Datenbank ein (Migrationen, statische Dateien, erster Administrator).

Für Docker, Portainer u. ä. nicht nötig: Dort kommt die Konfiguration aus Umgebungsvariablen
(siehe docker-compose.yml und README).
"""
import getpass
import os
import secrets
import shutil
import subprocess
import sys
from pathlib import Path
from urllib.parse import quote, urlsplit

from django.conf import settings
from django.contrib.auth import get_user_model
from django.core.mail import EmailMessage, get_connection
from django.core.management import BaseCommand, call_command

from bondamino.config import parse_database_url, read_env_file

DATABASES = [
    # (Schlüssel, Bezeichnung, URL-Schema, Standardport, Paketdatei für den Treiber)
    ('sqlite', 'SQLite (single file, for testing and small installations)', 'sqlite', '', None),
    ('postgres', 'PostgreSQL (recommended)', 'postgres', '5432', 'requirements-postgresql.txt'),
    ('mariadb', 'MariaDB / MySQL', 'mariadb', '3306', 'requirements-mariadb.txt'),
]
SMTP_BACKEND = 'django.core.mail.backends.smtp.EmailBackend'
CONSOLE_BACKEND = 'django.core.mail.backends.console.EmailBackend'


class Command(BaseCommand):
    help = 'Interactive setup: database, domain and email; writes .env and prepares the database.'

    def add_arguments(self, parser):
        parser.add_argument('--env-file', default=os.environ.get('BONDAMINO_ENV_FILE', str(settings.BASE_DIR / '.env')),
                            help='Configuration file to write (default: .env next to manage.py)')
        parser.add_argument('--skip-finish', action='store_true',
                            help='Only write the configuration; do not migrate the database')
        parser.add_argument('--finish', action='store_true',
                            help='Only prepare the database with the current configuration')

    # ----- Eingaben -----
    def ask(self, prompt, default=''):
        suffix = f' [{default}]' if default else ''
        answer = input(f'{prompt}{suffix}: ').strip()
        return answer or default

    def ask_secret(self, prompt, current):
        hint = ' [keep current]' if current else ''
        # Ohne Terminal (Eingabe per Pipe) liest getpass unter Windows nicht von stdin
        read = getpass.getpass if sys.stdin.isatty() else input
        answer = read(f'{prompt}{hint}: ')
        return answer or current

    def ask_yes_no(self, prompt, default):
        options = 'Y/n' if default else 'y/N'
        while True:
            answer = input(f'{prompt} [{options}]: ').strip().lower()
            if not answer:
                return default
            if answer in ('y', 'yes', 'j', 'ja'):
                return True
            if answer in ('n', 'no', 'nein'):
                return False

    def ask_choice(self, prompt, choices, default):
        """choices: [(Schlüssel, Bezeichnung)]; Eingabe der Nummer."""
        for number, (_, label) in enumerate(choices, 1):
            self.stdout.write(f'  {number}) {label}')
        keys = [key for key, _ in choices]
        while True:
            answer = self.ask(prompt, str(keys.index(default) + 1))
            if answer.isdigit() and 1 <= int(answer) <= len(choices):
                return keys[int(answer) - 1]

    def section(self, title):
        self.stdout.write('')
        self.stdout.write(self.style.MIGRATE_HEADING(title))

    # ----- Schritte -----
    def configure_database(self, current):
        self.section('1. Database')
        url = current.get('DATABASE_URL', 'sqlite:///db.sqlite3')
        parts = urlsplit(url)
        scheme = {'postgresql': 'postgres', 'mysql': 'mariadb'}.get(parts.scheme, parts.scheme)
        keys = [key for key, *_ in DATABASES]
        kind = self.ask_choice('Database', [(key, label) for key, label, *_ in DATABASES],
                               scheme if scheme in keys else 'sqlite')
        _, _, url_scheme, default_port, requirements = next(d for d in DATABASES if d[0] == kind)
        same = scheme == kind
        while True:
            if kind == 'sqlite':
                default_file = parts.path[1:] if same and parts.path else 'db.sqlite3'
                url = f"sqlite:///{self.ask('Database file (relative to manage.py or absolute)', default_file)}"
            else:
                host = self.ask('Server', (parts.hostname if same else '') or 'localhost')
                port = self.ask('Port', str(parts.port or '') if same and parts.port else default_port)
                name = self.ask('Database name', (parts.path.lstrip('/') if same else '') or 'bondamino')
                user = self.ask('Database user', (parts.username if same else '') or 'bondamino')
                password = self.ask_secret('Database password', parts.password if same and parts.password else '')
                url = f'{url_scheme}://{quote(user, safe="")}:{quote(password, safe="")}@{host}:{port}/{quote(name)}'
            error = self.check_database(url, requirements)
            if not error:
                self.stdout.write(self.style.SUCCESS('Connection OK.'))
                return url
            self.stdout.write(self.style.ERROR(f'Connection failed: {error}'))
            if not self.ask_yes_no('Enter the database settings again?', True):
                return url
            parts, same = urlsplit(url), True

    def check_database(self, url, requirements):
        """Fehlermeldung oder None, wenn die Verbindung klappt."""
        from django.db.utils import ConnectionHandler
        try:
            config = parse_database_url(url, settings.BASE_DIR)
        except ValueError as error:
            return str(error)
        if config['ENGINE'].endswith('sqlite3'):
            folder = Path(config['NAME']).parent
            return None if folder.is_dir() else f'Folder {folder} does not exist'
        try:
            connection = ConnectionHandler({'default': config})['default']
        except Exception as error:  # fehlender Treiber
            return f'{error}\nInstall the driver with: pip install -r {requirements}'
        try:
            connection.ensure_connection()
        except Exception as error:
            return str(error).strip()
        finally:
            connection.close()
        return None

    def configure_web(self, current):
        self.section('2. Web address')
        hosts = self.ask('Domain(s) under which Bondamino is reached, comma-separated',
                         current.get('ALLOWED_HOSTS', 'localhost'))
        https = self.ask_yes_no('Is Bondamino served via HTTPS (through a reverse proxy)?',
                                current.get('HTTPS', '0') == '1')
        debug = self.ask_yes_no('Development mode (DEBUG)? Never use it on a public server',
                                current.get('DEBUG', '0') == '1')
        time_zone = self.ask('Time zone', current.get('TIME_ZONE', 'Europe/Berlin'))
        return {'ALLOWED_HOSTS': hosts.replace(' ', ''), 'HTTPS': '1' if https else '0',
                'DEBUG': '1' if debug else '0', 'TIME_ZONE': time_zone}

    def configure_email(self, current):
        self.section('3. Email (invitations, forgotten passwords)')
        smtp = current.get('EMAIL_BACKEND') == SMTP_BACKEND
        kind = self.ask_choice('Sending', [
            ('smtp', 'Send via SMTP server'),
            ('console', 'Do not send; print emails in the server console (testing only)'),
        ], 'smtp' if smtp else 'console')
        if kind == 'console':
            return {'EMAIL_BACKEND': CONSOLE_BACKEND}
        values = {
            'EMAIL_BACKEND': SMTP_BACKEND,
            'EMAIL_HOST': self.ask('SMTP server', current.get('EMAIL_HOST', '')),
        }
        security = self.ask_choice('Encryption', [('tls', 'STARTTLS (usually port 587)'),
                                                  ('ssl', 'SSL/TLS (usually port 465)'),
                                                  ('none', 'None (usually port 25)')],
                                   'ssl' if current.get('EMAIL_USE_SSL') == '1'
                                   else 'none' if smtp and current.get('EMAIL_USE_TLS') != '1' else 'tls')
        values['EMAIL_USE_TLS'] = '1' if security == 'tls' else '0'
        values['EMAIL_USE_SSL'] = '1' if security == 'ssl' else '0'
        values['EMAIL_PORT'] = self.ask('Port', current.get('EMAIL_PORT') if smtp and current.get('EMAIL_PORT')
                                        else {'tls': '587', 'ssl': '465', 'none': '25'}[security])
        values['EMAIL_HOST_USER'] = self.ask('SMTP user (empty if none)', current.get('EMAIL_HOST_USER', ''))
        values['EMAIL_HOST_PASSWORD'] = (self.ask_secret('SMTP password', current.get('EMAIL_HOST_PASSWORD', ''))
                                         if values['EMAIL_HOST_USER'] else '')
        values['DEFAULT_FROM_EMAIL'] = self.ask('Sender address', current.get('DEFAULT_FROM_EMAIL', ''))
        recipient = self.ask('Send a test email to (empty to skip)')
        if recipient:
            self.send_test_email(values, recipient)
        return values

    def send_test_email(self, values, recipient):
        try:
            connection = get_connection(
                SMTP_BACKEND, host=values['EMAIL_HOST'], port=int(values['EMAIL_PORT']),
                username=values['EMAIL_HOST_USER'], password=values['EMAIL_HOST_PASSWORD'],
                use_tls=values['EMAIL_USE_TLS'] == '1', use_ssl=values['EMAIL_USE_SSL'] == '1', timeout=20,
            )
            EmailMessage('Bondamino test email', 'The email settings of Bondamino work.',
                         values['DEFAULT_FROM_EMAIL'] or None, [recipient], connection=connection).send()
            self.stdout.write(self.style.SUCCESS(f'Test email sent to {recipient}.'))
        except Exception as error:
            self.stdout.write(self.style.ERROR(f'Test email failed: {error}'))

    # ----- Datei -----
    def write_env(self, path, values):
        path = Path(path)
        if path.exists():
            backup = path.with_name(path.name + '.bak')
            shutil.copyfile(path, backup)
            self.stdout.write(f'Previous configuration saved as {backup}')
        lines = ['# Bondamino configuration, written by "python manage.py setup".',
                 '# Contains passwords and the secret key: do not share or commit this file.']
        lines += [f'{name}={value}' for name, value in values.items()]
        path.write_text('\n'.join(lines) + '\n', encoding='utf-8')
        try:
            path.chmod(0o600)
        except OSError:
            pass

    def summary(self, values):
        self.section('Summary')
        for name, value in values.items():
            if name in ('SECRET_KEY', 'EMAIL_HOST_PASSWORD') and value:
                value = '********'
            elif name == 'DATABASE_URL':
                password = urlsplit(value).password
                if password:
                    value = value.replace(f':{password}@', ':********@')
            self.stdout.write(f'  {name}={value}')

    # ----- Datenbank einrichten -----
    def finish(self):
        self.section('Preparing the database')
        call_command('migrate', interactive=False)
        call_command('createcachetable')
        call_command('collectstatic', interactive=False, verbosity=0)
        self.stdout.write(self.style.SUCCESS('Static files collected.'))
        if not get_user_model().objects.filter(is_superuser=True).exists():
            self.section('First administrator')
            call_command('createsuperuser')
        self.stdout.write('')
        self.stdout.write(self.style.SUCCESS('Bondamino is set up.'))

    def handle(self, *args, env_file, skip_finish, finish, **options):
        if finish:
            self.finish()
            return
        current = read_env_file(env_file)
        self.stdout.write(self.style.MIGRATE_HEADING('Bondamino setup'))
        self.stdout.write(f'Configuration file: {env_file}')
        self.stdout.write('Press Enter to accept the value in brackets.')

        values = {'DATABASE_URL': self.configure_database(current)}
        values.update(self.configure_web(current))
        values.update(self.configure_email(current))
        values['SECRET_KEY'] = current.get('SECRET_KEY') or secrets.token_urlsafe(50)
        # Weitere, von Hand eingetragene Werte bleiben erhalten
        values.update({name: value for name, value in current.items() if name not in values})

        self.summary(values)
        if not self.ask_yes_no(f'Write {env_file}?', True):
            self.stdout.write('Nothing written.')
            return
        self.write_env(env_file, values)
        self.stdout.write(self.style.SUCCESS(f'{env_file} written.'))
        if skip_finish:
            return
        # In einem neuen Prozess, damit die neue Konfiguration gilt. Werte, die aus der alten
        # .env-Datei in die Umgebung geladen wurden, dürfen die neuen nicht überdecken.
        child_env = {name: value for name, value in os.environ.items() if current.get(name) != value}
        child_env['BONDAMINO_ENV_FILE'] = str(env_file)
        subprocess.run([sys.executable, str(settings.BASE_DIR / 'manage.py'), 'setup', '--finish'],
                       env=child_env, check=False)
