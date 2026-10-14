"""Template tag ``static_v`` for cache-busting URLs of static files (``{% load versioned_static %}``)."""
import os

from django import template
from django.contrib.staticfiles import finders
from django.templatetags.static import static as static_url

register = template.Library()


@register.simple_tag
def static_v(path):
    """Like ``{% static %}``, but appends the file's modification time as a query parameter.

    This way browsers (Firefox in particular) do not keep serving a cached old version after code changes.

    Args:
        path: Path of the static file, as for ``{% static %}``.

    Returns:
        The static URL with ``?v=<mtime>``, or the plain static URL if the file cannot be found.
    """
    url = static_url(path)
    abs_path = finders.find(path)
    if abs_path:
        return f'{url}?v={int(os.path.getmtime(abs_path))}'
    return url
