"""Template context processors (registered in ``TEMPLATES`` in settings.py)."""


def languages(request):
    """Provides the language selector and the user's settings to all templates.

    The values are set on the request by ``middleware.LanguageMiddleware``.

    Returns:
        A dict with ``available_languages`` (languages created in the admin, for the navbar language selector)
        and ``user_settings`` (the user's ``UserSettings``, e.g. for the theme; None if not logged in).
    """
    return {
        'available_languages': getattr(request, 'available_languages', []),
        'user_settings': getattr(request, 'user_settings', None),
    }
