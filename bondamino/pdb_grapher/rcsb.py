"""Download of structure files from the RCSB Protein Data Bank.

Lets users load a structure by its four-character PDB ID instead of uploading a file. The
only public entry point is `fetch_pdb()`, which returns the raw file content as bytes; the
caller (the upload form) stores it like an uploaded file. All failures are raised as
`RCSBError` whose message is already translated and meant to be shown to the user.

Only the legacy PDB format is requested, because the parser in `utils.py` reads PDB
records; very large structures that RCSB only offers as mmCIF are therefore not available.
"""
import re
import urllib.error
import urllib.request

from django.utils.translation import gettext as _

# Download URL for a structure in PDB format; `{pdb_id}` is the upper-case four-character ID.
DOWNLOAD_URL = 'https://files.rcsb.org/download/{pdb_id}.pdb'
# Classic PDB ID: a digit 1-9 followed by three alphanumeric characters (e.g. 1LYZ).
PDB_ID_RE = re.compile(r'^[1-9][A-Za-z0-9]{3}$')
TIMEOUT = 20  # Network timeout in seconds.
# Maximum file size in bytes; also used by the upload form (forms.py) for uploaded files.
MAX_SIZE = 20 * 1024 * 1024


class RCSBError(Exception):
    """Error while fetching a structure; the message is translated and meant for the user."""


def normalize_pdb_id(pdb_id):
    """Validates a PDB ID and normalizes it to upper case.

    Args:
        pdb_id: User input, e.g. ' 1lyz '. Surrounding whitespace is ignored.

    Returns:
        The ID in upper case, e.g. '1LYZ'.

    Raises:
        RCSBError: If the ID does not have the format of a PDB ID.
    """
    pdb_id = pdb_id.strip().upper()
    if not PDB_ID_RE.match(pdb_id):
        raise RCSBError(_('Ungültige PDB-ID. Eine PDB-ID besteht aus vier Zeichen, z. B. 1LYZ.'))
    return pdb_id


def fetch_pdb(pdb_id):
    """Downloads a structure in PDB format from RCSB.

    Args:
        pdb_id: PDB ID of the structure; normalized with `normalize_pdb_id()`.

    Returns:
        The file content as bytes (at most `MAX_SIZE` bytes).

    Raises:
        RCSBError: If the ID is invalid, RCSB has no PDB-format file for it (HTTP 404),
            RCSB rejects the request, the server is unreachable or times out, or the file
            is larger than `MAX_SIZE`.
    """
    pdb_id = normalize_pdb_id(pdb_id)
    request = urllib.request.Request(DOWNLOAD_URL.format(pdb_id=pdb_id), headers={'User-Agent': 'Bondamino'})
    try:
        with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
            # Read one byte more than allowed so that an oversized file can be detected
            # without loading all of it.
            content = response.read(MAX_SIZE + 1)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            # Unknown ID, or a very large structure that RCSB only provides as mmCIF.
            raise RCSBError(_('Für %(id)s gibt es bei RCSB keine Datei im PDB-Format.') % {'id': pdb_id}) from e
        raise RCSBError(_('RCSB hat den Abruf von %(id)s abgelehnt (HTTP %(code)s).')
                        % {'id': pdb_id, 'code': e.code}) from e
    except (urllib.error.URLError, TimeoutError) as e:
        raise RCSBError(_('RCSB ist nicht erreichbar. Bitte später erneut versuchen.')) from e
    if len(content) > MAX_SIZE:
        raise RCSBError(_('Die Datei für %(id)s ist größer als 20 MB.') % {'id': pdb_id})
    return content
