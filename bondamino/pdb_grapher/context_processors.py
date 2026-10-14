def languages(request):
    """Im Admin angelegte Sprachen für die Sprachauswahl in der Navbar und die Einstellungen des
    Users (Darstellung); ohne Anmeldung None."""
    return {
        'available_languages': getattr(request, 'available_languages', []),
        'user_settings': getattr(request, 'user_settings', None),
    }
