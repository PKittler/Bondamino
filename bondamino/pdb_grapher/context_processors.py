def languages(request):
    """Im Admin angelegte Sprachen für die Sprachauswahl in der Navbar."""
    return {'available_languages': getattr(request, 'available_languages', [])}
