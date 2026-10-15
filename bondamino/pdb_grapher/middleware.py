"""Language and time zone selection for each request.

Bondamino replaces Django's ``LocaleMiddleware`` with ``LanguageMiddleware`` so that only languages created in the
admin (``models.Language``) can be activated, and so that the language stored in the user's settings takes
precedence over cookie and browser. The same middleware activates the user's time zone, so that dates and times in
templates are shown in it.
"""
import zoneinfo

from django.conf import settings
from django.utils import timezone, translation
from django.utils.cache import patch_vary_headers
from django.utils.translation.trans_real import parse_accept_lang_header

from .models import Language, UserSettings


def activate_time_zone(name):
    """Activates the time zone ``name`` for the current request, or the server default if empty or unknown.

    Args:
        name: IANA time zone name from the user's settings, or ''.
    """
    if not name:
        timezone.deactivate()
        return
    try:
        timezone.activate(zoneinfo.ZoneInfo(name))
    except (zoneinfo.ZoneInfoNotFoundError, ValueError):
        # E.g. a name entered by hand in the admin, or removed from the time zone database
        timezone.deactivate()


def select_language(request, codes, preferred=''):
    """Chooses the language for a request.

    Order: the user's settings (``preferred``), then the selected language (cookie), then the best match from the
    browser's Accept-Language header, then ``settings.LANGUAGE_CODE``, then the first available language.

    Args:
        request: The current request.
        codes: Codes of the languages created in the admin; only these are considered.
        preferred: Language code from the user's settings, empty if none.

    Returns:
        A language code from ``codes``, or ``settings.LANGUAGE_CODE`` if ``codes`` is empty.
    """
    if not codes:
        return settings.LANGUAGE_CODE
    for chosen in (preferred, request.COOKIES.get(settings.LANGUAGE_COOKIE_NAME)):
        if chosen in codes:
            return chosen
    for accepted, _ in parse_accept_lang_header(request.META.get('HTTP_ACCEPT_LANGUAGE', '')):
        accepted = accepted.lower()
        if accepted in codes:
            return accepted
        # "en-us" or "en" also matches "en-gb"
        generic = accepted.split('-')[0]
        for code in codes:
            if code.split('-')[0] == generic:
                return code
    return settings.LANGUAGE_CODE if settings.LANGUAGE_CODE in codes else codes[0]


class LanguageMiddleware:
    """Replaces Django's LocaleMiddleware: only languages created in the admin are activated.

    Also activates the time zone from the user's settings (``activate_time_zone``).

    Stores the available languages and the user's settings on the request as ``available_languages`` and
    ``user_settings`` (exposed to templates by ``context_processors.languages``, e.g. for the navbar). Must come
    after the AuthenticationMiddleware in ``MIDDLEWARE`` because the settings belong to the user.
    """

    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        request.available_languages = list(Language.objects.all())
        request.user_settings = UserSettings.of(request.user) if request.user.is_authenticated else None
        preferred = request.user_settings.language if request.user_settings else ''
        language = select_language(request, [lang.code for lang in request.available_languages], preferred)
        translation.activate(language)
        request.LANGUAGE_CODE = translation.get_language()
        activate_time_zone(request.user_settings.time_zone if request.user_settings else '')

        response = self.get_response(request)
        # The response depends on the browser language and the language cookie, so caches must not mix them up
        patch_vary_headers(response, ('Accept-Language', 'Cookie'))
        response.headers.setdefault('Content-Language', request.LANGUAGE_CODE)
        return response
