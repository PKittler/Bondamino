import os

from django import template
from django.contrib.staticfiles import finders
from django.templatetags.static import static as static_url

register = template.Library()


@register.simple_tag
def static_v(path):
    """Wie {% static %}, hängt aber die Änderungszeit der Datei als Query-Parameter an,
    damit Browser (insbesondere Firefox) nach Code-Änderungen keine gecachte alte
    Version mehr ausliefern."""
    url = static_url(path)
    abs_path = finders.find(path)
    if abs_path:
        return f'{url}?v={int(os.path.getmtime(abs_path))}'
    return url
