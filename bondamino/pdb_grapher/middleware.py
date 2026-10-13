from django.conf import settings
from django.utils import translation
from django.utils.cache import patch_vary_headers
from django.utils.translation.trans_real import parse_accept_lang_header

from .models import Language


def select_language(request, codes):
    """Gewählte Sprache (Cookie), sonst die passendste laut Browser, sonst die Standardsprache.
    Berücksichtigt werden nur die im Admin angelegten Sprachen `codes`."""
    if not codes:
        return settings.LANGUAGE_CODE
    chosen = request.COOKIES.get(settings.LANGUAGE_COOKIE_NAME)
    if chosen in codes:
        return chosen
    for accepted, _ in parse_accept_lang_header(request.META.get('HTTP_ACCEPT_LANGUAGE', '')):
        accepted = accepted.lower()
        if accepted in codes:
            return accepted
        # „en-us“ oder „en“ passt auch zu „en-gb“
        generic = accepted.split('-')[0]
        for code in codes:
            if code.split('-')[0] == generic:
                return code
    return settings.LANGUAGE_CODE if settings.LANGUAGE_CODE in codes else codes[0]


class LanguageMiddleware:
    """Ersetzt Djangos LocaleMiddleware: Aktiviert werden nur Sprachen, die im Admin angelegt
    sind. Die Sprachen merkt sich der Request für die Auswahl in der Navbar."""

    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        request.available_languages = list(Language.objects.all())
        language = select_language(request, [lang.code for lang in request.available_languages])
        translation.activate(language)
        request.LANGUAGE_CODE = translation.get_language()

        response = self.get_response(request)
        patch_vary_headers(response, ('Accept-Language', 'Cookie'))
        response.headers.setdefault('Content-Language', request.LANGUAGE_CODE)
        return response
