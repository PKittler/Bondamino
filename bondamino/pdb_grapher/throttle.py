"""Rate limiting of login attempts and password reset requests.

The counters live in the cache (CACHES in settings.py: database cache, so that all server processes see the same
counters). Keys contain only a SHA-256 hash of IP address, user name or email address. If the cache is not
available, e.g. because the table is still missing after an update (``manage.py createcachetable``), nothing is
limited and an error is logged.

Used by ``forms.ThrottledLoginMixin`` (login) and ``forms.PasswordResetForm`` (forgotten password).
"""
import hashlib
import logging

from django.conf import settings
from django.core.cache import cache
from django.db import DatabaseError

logger = logging.getLogger(__name__)

# Login: failed attempts per IP address and user name, and per IP address in total
LOGIN_FAILURES_PER_USER = 5
LOGIN_FAILURES_PER_IP = 20
LOGIN_LOCK_SECONDS = 15 * 60

# Forgotten password: requests per email address and per IP address
RESET_REQUESTS_PER_EMAIL = 3
RESET_REQUESTS_PER_IP = 10
RESET_WINDOW_SECONDS = 60 * 60


def client_ip(request):
    """Returns the client's IP address.

    Behind a reverse proxy (``HTTPS=1``) this is the last entry of X-Forwarded-For, which the proxy appends
    itself; earlier entries can be forged by the client. Otherwise REMOTE_ADDR.
    """
    if settings.HTTPS:
        forwarded = request.META.get('HTTP_X_FORWARDED_FOR', '')
        if forwarded.strip():
            return forwarded.split(',')[-1].strip()
    return request.META.get('REMOTE_ADDR', '')


def _key(kind, *parts):
    """Builds a cache key from a hash of the case-folded parts (no personal data in the cache)."""
    digest = hashlib.sha256('\0'.join(str(p).casefold() for p in parts).encode()).hexdigest()
    return f'throttle:{kind}:{digest}'


def _count(key):
    """Returns the counter value, 0 if unset or the cache is unavailable."""
    try:
        return cache.get(key, 0)
    except DatabaseError:
        logger.exception('Throttle cache not available; run "manage.py createcachetable"')
        return 0


def _increment(key, timeout):
    """Increments a counter; the time window starts with the first entry and is not extended by later ones."""
    try:
        if not cache.add(key, 1, timeout):
            cache.incr(key)
    except ValueError:  # expired between add and incr
        cache.add(key, 1, timeout)
    except DatabaseError:
        logger.exception('Throttle cache not available; run "manage.py createcachetable"')


def _delete(key):
    """Resets a counter."""
    try:
        cache.delete(key)
    except DatabaseError:
        logger.exception('Throttle cache not available; run "manage.py createcachetable"')


def login_locked(request, username):
    """Checks whether logins are currently locked for this client.

    Args:
        request: The login request (for the client IP).
        username: The entered user name.

    Returns:
        True while there are too many failed attempts for this IP address and user name, or for this IP address
        in total.
    """
    ip = client_ip(request)
    return (_count(_key('login', ip, username)) >= LOGIN_FAILURES_PER_USER
            or _count(_key('login-ip', ip)) >= LOGIN_FAILURES_PER_IP)


def login_failed(request, username):
    """Counts a failed login for the IP address and user name and for the IP address in total."""
    ip = client_ip(request)
    _increment(_key('login', ip, username), LOGIN_LOCK_SECONDS)
    _increment(_key('login-ip', ip), LOGIN_LOCK_SECONDS)


def login_succeeded(request, username):
    """Resets the failure counter for IP address and user name.

    The per-IP total is not reset, so a successful login with one account does not lift the limit for guessing
    the passwords of other accounts from the same IP address.
    """
    _delete(_key('login', client_ip(request), username))


def reset_allowed(request, email):
    """Counts a password reset request and checks whether an email may still be sent.

    Args:
        request: The reset request (for the client IP).
        email: The entered email address.

    Returns:
        True if neither the limit per email address nor the limit per IP address was reached before this
        request.
    """
    keys = ((_key('reset', email), RESET_REQUESTS_PER_EMAIL), (_key('reset-ip', client_ip(request)), RESET_REQUESTS_PER_IP))
    allowed = all(_count(key) < limit for key, limit in keys)
    for key, _limit in keys:
        _increment(key, RESET_WINDOW_SECONDS)
    return allowed
