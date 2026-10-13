"""Abruf von Strukturdateien aus der RCSB Protein Data Bank."""
import re
import urllib.error
import urllib.request

from django.utils.translation import gettext as _

DOWNLOAD_URL = 'https://files.rcsb.org/download/{pdb_id}.pdb'
PDB_ID_RE = re.compile(r'^[1-9][A-Za-z0-9]{3}$')
TIMEOUT = 20  # Sekunden
MAX_SIZE = 20 * 1024 * 1024  # wie beim Datei-Upload


class RCSBError(Exception):
    """Fehler beim Abruf; die Meldung ist für den User bestimmt."""


def normalize_pdb_id(pdb_id):
    """Gibt die ID in Großbuchstaben zurück oder wirft RCSBError bei ungültigem Format."""
    pdb_id = pdb_id.strip().upper()
    if not PDB_ID_RE.match(pdb_id):
        raise RCSBError(_('Ungültige PDB-ID. Eine PDB-ID besteht aus vier Zeichen, z. B. 1LYZ.'))
    return pdb_id


def fetch_pdb(pdb_id):
    """Lädt die Struktur `pdb_id` im PDB-Format von RCSB und gibt den Inhalt als Bytes zurück."""
    pdb_id = normalize_pdb_id(pdb_id)
    request = urllib.request.Request(DOWNLOAD_URL.format(pdb_id=pdb_id), headers={'User-Agent': 'Bondamino'})
    try:
        with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
            content = response.read(MAX_SIZE + 1)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            # Sehr große Strukturen gibt es bei RCSB nur als mmCIF, nicht im PDB-Format.
            raise RCSBError(_('Für %(id)s gibt es bei RCSB keine Datei im PDB-Format.') % {'id': pdb_id}) from e
        raise RCSBError(_('RCSB hat den Abruf von %(id)s abgelehnt (HTTP %(code)s).')
                        % {'id': pdb_id, 'code': e.code}) from e
    except (urllib.error.URLError, TimeoutError) as e:
        raise RCSBError(_('RCSB ist nicht erreichbar. Bitte später erneut versuchen.')) from e
    if len(content) > MAX_SIZE:
        raise RCSBError(_('Die Datei für %(id)s ist größer als 20 MB.') % {'id': pdb_id})
    return content
