"""Begrenzung von Anmeldeversuchen und Anfragen nach einem neuen Passwort.

Die Zähler liegen im Cache (CACHES in settings.py: Datenbank-Cache, damit alle Server-Prozesse
dieselben Zähler sehen). Ist der Cache nicht erreichbar, z. B. weil die Tabelle nach einem
Update noch fehlt (`manage.py createcachetable`), wird nicht begrenzt und ein Fehler protokolliert.
"""
import hashlib
import logging

from django.conf import settings
from django.core.cache import cache
from django.db import DatabaseError

logger = logging.getLogger(__name__)

# Anmeldung: Fehlversuche je IP-Adresse und Benutzername sowie je IP-Adresse insgesamt
LOGIN_FAILURES_PER_USER = 5
LOGIN_FAILURES_PER_IP = 20
LOGIN_LOCK_SECONDS = 15 * 60

# Passwort vergessen: Anfragen je E-Mail-Adresse und je IP-Adresse
RESET_REQUESTS_PER_EMAIL = 3
RESET_REQUESTS_PER_IP = 10
RESET_WINDOW_SECONDS = 60 * 60


def client_ip(request):
    """IP-Adresse des Clients. Hinter einem Reverse Proxy (HTTPS=1) der letzte Eintrag von
    X-Forwarded-For, den der Proxy selbst anhängt; frühere Einträge kann der Client fälschen."""
    if settings.HTTPS:
        forwarded = request.META.get('HTTP_X_FORWARDED_FOR', '')
        if forwarded.strip():
            return forwarded.split(',')[-1].strip()
    return request.META.get('REMOTE_ADDR', '')


def _key(kind, *parts):
    digest = hashlib.sha256('\0'.join(str(p).casefold() for p in parts).encode()).hexdigest()
    return f'throttle:{kind}:{digest}'


def _count(key):
    try:
        return cache.get(key, 0)
    except DatabaseError:
        logger.exception('Throttle cache not available; run "manage.py createcachetable"')
        return 0


def _increment(key, timeout):
    """Zähler erhöhen; das Zeitfenster beginnt mit dem ersten Eintrag."""
    try:
        if not cache.add(key, 1, timeout):
            cache.incr(key)
    except ValueError:  # zwischen add und incr abgelaufen
        cache.add(key, 1, timeout)
    except DatabaseError:
        logger.exception('Throttle cache not available; run "manage.py createcachetable"')


def _delete(key):
    try:
        cache.delete(key)
    except DatabaseError:
        logger.exception('Throttle cache not available; run "manage.py createcachetable"')


def login_locked(request, username):
    """True, solange für diese IP-Adresse (und diesen Benutzernamen) zu viele Fehlversuche vorliegen."""
    ip = client_ip(request)
    return (_count(_key('login', ip, username)) >= LOGIN_FAILURES_PER_USER
            or _count(_key('login-ip', ip)) >= LOGIN_FAILURES_PER_IP)


def login_failed(request, username):
    ip = client_ip(request)
    _increment(_key('login', ip, username), LOGIN_LOCK_SECONDS)
    _increment(_key('login-ip', ip), LOGIN_LOCK_SECONDS)


def login_succeeded(request, username):
    _delete(_key('login', client_ip(request), username))


def reset_allowed(request, email):
    """Zählt die Anfrage und gibt zurück, ob noch eine E-Mail gesendet werden darf."""
    keys = ((_key('reset', email), RESET_REQUESTS_PER_EMAIL), (_key('reset-ip', client_ip(request)), RESET_REQUESTS_PER_IP))
    allowed = all(_count(key) < limit for key, limit in keys)
    for key, _limit in keys:
        _increment(key, RESET_WINDOW_SECONDS)
    return allowed
