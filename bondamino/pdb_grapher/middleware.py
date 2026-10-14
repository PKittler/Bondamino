from django.conf import settings
from django.utils import translation
from django.utils.cache import patch_vary_headers
from django.utils.translation.trans_real import parse_accept_lang_header

from .models import Language, UserSettings


def select_language(request, codes, preferred=''):
    """Sprache aus den Einstellungen des Users (`preferred`), sonst die gewählte (Cookie), sonst
    die passendste laut Browser, sonst die Standardsprache. Berücksichtigt werden nur die im
    Admin angelegten Sprachen `codes`."""
    if not codes:
        return settings.LANGUAGE_CODE
    for chosen in (preferred, request.COOKIES.get(settings.LANGUAGE_COOKIE_NAME)):
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
    sind. Die Sprachen und die Einstellungen des Users merkt sich der Request für die Navbar.
    Steht nach der AuthenticationMiddleware, weil die Einstellungen am User hängen."""

    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        request.available_languages = list(Language.objects.all())
        request.user_settings = UserSettings.of(request.user) if request.user.is_authenticated else None
        preferred = request.user_settings.language if request.user_settings else ''
        language = select_language(request, [lang.code for lang in request.available_languages], preferred)
        translation.activate(language)
        request.LANGUAGE_CODE = translation.get_language()

        response = self.get_response(request)
        patch_vary_headers(response, ('Accept-Language', 'Cookie'))
        response.headers.setdefault('Content-Language', request.LANGUAGE_CODE)
        return response
