import json
import math
import statistics
from pathlib import Path

BACKBONE_BONDS = [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('C', 'OXT')]

# Längster C–N-Abstand (Å), der noch als Peptidbindung gilt (Soll ≈ 1,33 Å). Größere
# Abstände zwischen aufeinanderfolgenden Residuen sind Kettenlücken, z. B. durch in der
# Struktur fehlende Residuen, und keine Bindung.
MAX_PEPTIDE_BOND_LENGTH = 2.0

# Längster SG–SG-Abstand (Å), der ohne SSBOND-Eintrag als Disulfidbrücke gilt (typisch ≈ 2,0 Å).
# Mit SSBOND-Eintrag gilt die Brücke unabhängig vom Abstand.
MAX_DISULFIDE_LENGTH = 2.5

SIDECHAIN_BONDS = {
    'ALA': [('CA', 'CB')],
    'ARG': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD'), ('CD', 'NE'), ('NE', 'CZ'), ('CZ', 'NH1'), ('CZ', 'NH2')],
    'ASN': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'OD1'), ('CG', 'ND2')],
    'ASP': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'OD1'), ('CG', 'OD2')],
    'CYS': [('CA', 'CB'), ('CB', 'SG')],
    'GLN': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD'), ('CD', 'OE1'), ('CD', 'NE2')],
    'GLU': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD'), ('CD', 'OE1'), ('CD', 'OE2')],
    'GLY': [],
    'HIS': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'ND1'), ('CG', 'CD2'), ('ND1', 'CE1'), ('CD2', 'NE2'), ('CE1', 'NE2')],
    'ILE': [('CA', 'CB'), ('CB', 'CG1'), ('CB', 'CG2'), ('CG1', 'CD1')],
    'LEU': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD1'), ('CG', 'CD2')],
    'LYS': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD'), ('CD', 'CE'), ('CE', 'NZ')],
    'MET': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'SD'), ('SD', 'CE')],
    'PHE': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD1'), ('CG', 'CD2'), ('CD1', 'CE1'), ('CD2', 'CE2'), ('CE1', 'CZ'), ('CE2', 'CZ')],
    'PRO': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD'), ('CD', 'N')],
    'SER': [('CA', 'CB'), ('CB', 'OG')],
    'THR': [('CA', 'CB'), ('CB', 'OG1'), ('CB', 'CG2')],
    'TRP': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD1'), ('CG', 'CD2'), ('CD1', 'NE1'), ('CD2', 'CE2'), ('CD2', 'CE3'), ('NE1', 'CE2'), ('CE2', 'CZ2'), ('CE3', 'CZ3'), ('CZ2', 'CH2'), ('CZ3', 'CH2')],
    'TYR': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD1'), ('CG', 'CD2'), ('CD1', 'CE1'), ('CD2', 'CE2'), ('CE1', 'CZ'), ('CE2', 'CZ'), ('CZ', 'OH')],
    'VAL': [('CA', 'CB'), ('CB', 'CG1'), ('CB', 'CG2')],
}

# Referenzwerte (Mittelwert, Standardabweichung) für Bindungslängen (Å) und Bindungswinkel (°).
# Rückgrat: Engh & Huber (2001) mit den Zuordnungsregeln aus Dangle (Richardson-Labor,
# BSD-3-Clause); Seitenketten: CCP4 Monomer Library (LGPL-3.0). Herkunft und Lizenzen
# stehen in reference_data/README.md.
REFERENCE_SOURCES = {
    'backbone': 'Engh & Huber (2001), International Tables for Crystallography F, 382–392',
    'sidechain': 'CCP4 Monomer Library',
}
REFERENCE_DATA_DIR = Path(__file__).resolve().parent / 'reference_data'
BACKBONE_ATOMS = {'N', 'CA', 'C', 'O', 'OXT'}

# Klassen: das Residuum selbst ('GLY', 'PRO'), über die Peptidbindung das folgende
# ('X-GLY', 'X-PRO', 'cis-X-PRO' bei cis-Peptidbindung) oder das vorige ('GLY-X');
# 'other' gilt sonst. Die erste passende Klasse in der Liste gilt (Vorrang wie in Dangle).
BACKBONE_BOND_REFERENCES = {
    ('N', 'CA'): [('GLY', 1.456, 0.015), ('PRO', 1.468, 0.017), ('other', 1.459, 0.020)],
    ('CA', 'C'): [('GLY', 1.514, 0.016), ('PRO', 1.524, 0.020), ('other', 1.525, 0.026)],
    ('C', 'O'): [('GLY', 1.232, 0.016), ('PRO', 1.228, 0.020), ('other', 1.229, 0.019)],
    ('C', 'N'): [('X-GLY', 1.326, 0.018), ('X-PRO', 1.338, 0.019), ('other', 1.336, 0.023)],
}
# Atome in Reihenfolge vorheriges Residuum -> folgendes (C–N ist die Peptidbindung)
BACKBONE_ANGLE_REFERENCES = {
    ('N', 'CA', 'C'): [('GLY', 113.1, 2.5), ('PRO', 112.1, 2.6), ('other', 111.0, 2.7)],
    ('CA', 'C', 'O'): [('GLY', 120.6, 1.8), ('PRO', 120.2, 2.4), ('other', 120.1, 2.1)],
    ('CA', 'C', 'N'): [('GLY-X', 116.2, 2.0), ('X-PRO', 117.1, 2.8), ('other', 117.2, 2.2)],
    ('O', 'C', 'N'): [('X-PRO', 121.1, 1.9), ('GLY-X', 123.2, 1.7), ('other', 122.7, 1.6)],
    ('C', 'N', 'CA'): [('cis-X-PRO', 127.0, 2.4), ('X-PRO', 119.3, 1.5), ('X-GLY', 122.3, 2.1), ('other', 121.7, 2.5)],
}
CIS_MAX_OMEGA = 30  # |ω| < 30°: cis-Peptidbindung (wie im Browser)


def _load_sidechain_references():
    """Seitenketten-Referenzen je Aminosäure (eigener Residuenname, auch modifizierte):
    Bindungen {(Residuum, {Atom, Atom}): [...]}, Winkel {(Residuum, Mitte, {Ende, Ende}): [...]}."""
    data = json.loads((REFERENCE_DATA_DIR / 'ccp4_sidechains.json').read_text(encoding='utf-8'))['residues']
    bonds, angles = {}, {}
    for res, entry in data.items():
        for a, b, mean, sigma in entry['bonds']:
            bonds[(res, frozenset((a, b)))] = [mean, sigma, f'{res}:{a}-{b}']
        for a, center, c, mean, sigma in entry['angles']:
            angles[(res, center, frozenset((a, c)))] = [mean, sigma, f'{res}:{a}-{center}-{c}']
    return bonds, angles


SIDECHAIN_BOND_REFERENCES, SIDECHAIN_ANGLE_REFERENCES = _load_sidechain_references()


def _classes(own=None, prev=None, next_=None, cis=False):
    """Zutreffende Klassen (siehe BACKBONE_BOND_REFERENCES) für ein Residuum bzw. eine Peptidbindung."""
    classes = {'other'}
    if own in ('GLY', 'PRO'):
        classes.add(own)
    if prev == 'GLY':
        classes.add('GLY-X')
    if next_ in ('GLY', 'PRO'):
        classes.add(f'X-{next_}')
        if next_ == 'PRO' and cis:
            classes.add('cis-X-PRO')
    return classes


def _backbone_reference(table, atoms, classes):
    """Erster passender Rückgrat-Referenzwert als [Mittelwert, Standardabweichung, Gruppe] oder None."""
    for cls, mean, sigma in table.get(atoms, []):
        if cls in classes:
            return [mean, sigma, f"backbone:{'-'.join(atoms)}:{cls}"]
    return None


def _reference(atom1, atom2, bond_type, res_name, res_class, next_class=None):
    """Referenz einer Bindung als [Mittelwert, Standardabweichung, Gruppe] oder None.
    `res_name` ist der Name des Residuums, `res_class` der seiner (Ausgangs-)Aminosäure,
    `next_class` der des folgenden Residuums (nur für die Peptidbindung)."""
    if bond_type == 'peptide':
        if res_class not in SIDECHAIN_BONDS or next_class not in SIDECHAIN_BONDS:
            return None
        return _backbone_reference(BACKBONE_BOND_REFERENCES, ('C', 'N'), _classes(next_=next_class))
    if {atom1, atom2} - BACKBONE_ATOMS:
        return SIDECHAIN_BOND_REFERENCES.get((res_name, frozenset((atom1, atom2))))
    if res_class not in SIDECHAIN_BONDS:
        return None  # keine Aminosäure oder modifizierte ohne bekannte Ausgangs-Aminosäure
    atoms = (atom1, atom2) if (atom1, atom2) in BACKBONE_BOND_REFERENCES else (atom2, atom1)
    return _backbone_reference(BACKBONE_BOND_REFERENCES, atoms, _classes(own=res_class))


def parse_header(file_path):
    """Strukturmetadaten aus dem PDB-Header: Methode (EXPDTA), Auflösung (REMARK 2) und
    R-Werte (REMARK 3). Fehlende Angaben sind None."""
    info = {'id_code': None, 'method': None, 'resolution': None, 'r_work': None, 'r_free': None}
    with open(file_path, 'r') as f:
        for line in f:
            record = line[:6].strip()
            if record in ('ATOM', 'HETATM', 'MODEL'):
                break  # Header ist zu Ende
            if record == 'HEADER':
                info['id_code'] = line[62:66].strip() or None
            elif record == 'EXPDTA':
                method = line[10:].strip()
                info['method'] = f"{info['method']}; {method}" if info['method'] else method
            elif line.startswith('REMARK   2 RESOLUTION.'):
                value = line[22:].split()
                info['resolution'] = _float_or_none(value[0]) if value else None
            elif line.startswith('REMARK   3') and ':' in line:
                label, _, value = line[10:].partition(':')
                label = ' '.join(label.split())
                if label in ('R VALUE (WORKING SET)', 'R VALUE (WORKING + TEST SET)') and info['r_work'] is None:
                    info['r_work'] = _float_or_none(value.strip())
                elif label == 'FREE R VALUE' and info['r_free'] is None:
                    info['r_free'] = _float_or_none(value.strip())
    return info


# Ausgangs-Aminosäure gängiger modifizierter Aminosäuren; gilt nur, wenn die Datei keine
# MODRES-Einträge dazu hat (Dateien von RCSB bringen sie mit).
MODIFIED_PARENTS = {
    'MSE': 'MET', 'FME': 'MET',
    'SEP': 'SER', 'TPO': 'THR', 'PTR': 'TYR',
    'CSO': 'CYS', 'CME': 'CYS', 'OCS': 'CYS', 'CSD': 'CYS',
    'HYP': 'PRO',
    'MLY': 'LYS', 'M3L': 'LYS', 'ALY': 'LYS', 'KCX': 'LYS', 'LLP': 'LYS',
    'CGU': 'GLU', 'NEP': 'HIS',
}

# Kovalenzradien in Å (Cordero et al. 2008) für die Bindungserkennung über Abstände
COVALENT_RADII = {
    'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05, 'SE': 1.20, 'P': 1.07,
    'F': 0.57, 'CL': 1.02, 'BR': 1.20, 'I': 1.39, 'B': 0.84,
}
DEFAULT_COVALENT_RADIUS = 0.80
# Zwei Atome gelten als gebunden, wenn ihr Abstand höchstens die Summe der Radien plus
# diese Toleranz beträgt (übliche Werte 0,4–0,45 Å)
BOND_TOLERANCE = 0.45

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


def _float_or_none(text):
    try:
        return float(text)
    except ValueError:
        return None


def _int_or_none(text):
    try:
        return int(text)
    except ValueError:
        return None


def list_models(file_path):
    """Gibt die im PDB-File vorkommenden MODEL-Nummern in Auftrittsreihenfolge zurück.

    Enthält die Datei keine MODEL-Records, wird eine leere Liste zurückgegeben
    (die Datei wird dann als einzelnes, implizites Modell behandelt).
    """
    models = []
    with open(file_path, 'r') as f:
        for line in f:
            if line[:6].strip() != 'MODEL':
                continue
            try:
                model_num = int(line[10:14].strip())
            except ValueError:
                model_num = len(models) + 1
            models.append(model_num)
    return models


def list_altlocs(file_path):
    """Gibt die in der Datei vorkommenden Kennungen alternativer Konformationen
    (altLoc, z. B. ['A', 'B']) sortiert zurück; leer, wenn es keine gibt."""
    altlocs = set()
    with open(file_path, 'r') as f:
        for line in f:
            if line[:6].strip() in ('ATOM', 'HETATM') and len(line) > 16 and line[16].strip():
                altlocs.add(line[16])
    return sorted(altlocs)


def _choose_altloc(atoms, preferred):
    """Konformation eines Residuums: `preferred`, falls das Residuum sie hat, sonst die
    mit der höchsten mittleren Belegung (bei Gleichstand die zuerst auftretende).
    `atoms` sind (Atomname, altLoc, Residuenname, Atomdaten) in Dateireihenfolge.
    Gibt '' zurück, wenn das Residuum keine Alternativen hat."""
    occupancies = {}
    for _, alt_loc, _, atom in atoms:
        if alt_loc:
            occupancies.setdefault(alt_loc, []).append(atom['occupancy'] or 0.0)
    if not occupancies:
        return ''
    if preferred in occupancies:
        return preferred
    # max() liefert bei Gleichstand das erste Element, dict behält die Dateireihenfolge
    return max(occupancies, key=lambda a: statistics.fmean(occupancies[a]))


class Residue(dict):
    """Atome eines Residuums ({Atomname: Atomdaten}) mit Zusatzangaben:

    `conformations` ordnet jeder Konformation des Residuums in der Datei ihren
    Residuennamen zu, z. B. {"A": "PRO", "B": "SER"} (leer ohne Alternativen).
    `parent` ist bei modifizierten Aminosäuren die Ausgangs-Aminosäure (z. B. "MET" für
    MSE), `sidechain_bonds` deren Seitenketten-Bindungen als Atomnamen-Paare.
    """

    def __init__(self, *args, conformations=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.conformations = conformations or {}
        self.parent = None
        self.sidechain_bonds = []
        # (chain, seq, icode) der Cysteine, mit denen laut SSBOND eine Disulfidbrücke besteht
        self.ssbond_partners = set()


def _select_conformation(atoms, preferred):
    """Atome eines Residuums in einer einzigen Konformation: alle ohne altLoc plus die der
    gewählten Alternative. Liefert (Residuenname, {Atomname: Atomdaten}); der Name stammt
    aus der gewählten Alternative, da er sich bei Mikroheterogenität (z. B. A = PRO,
    B = SER an derselben Position) unterscheidet."""
    chosen = _choose_altloc(atoms, preferred)
    selected = [a for a in atoms if a[1] in ('', chosen)]
    res_name = next((a[2] for a in selected if chosen and a[1] == chosen), selected[0][2])
    conformations = {}
    for _, alt_loc, name, _ in atoms:
        if alt_loc:
            conformations.setdefault(alt_loc, name)
    residue = Residue(conformations=dict(sorted(conformations.items())))
    for atom_name, _, _, atom in selected:
        residue.setdefault(atom_name, atom)
    return res_name, residue


def parse_pdb(file_path, model_id=None, altloc=None):
    """Liest Atome aus einer PDB-Datei.

    Wenn `model_id` angegeben ist, werden nur Atome des entsprechenden MODEL-Blocks
    berücksichtigt. Enthält die Datei keine MODEL-Records, hat `model_id` keine
    Auswirkung und alle ATOM/HETATM-Zeilen werden gelesen.

    Residuen werden über Chain, Sequenznummer und Insertion-Code unterschieden (52 und
    52A sind verschiedene Residuen). Je Residuum wird genau eine Konformation verwendet:
    `altloc`, falls das Residuum sie hat, sonst die mit der höchsten Belegung.

    Schlüssel des Ergebnisses: (chain_id, res_seq, i_code, res_name).
    """
    # (chain, seq, icode) -> [(Atomname, altLoc, Residuenname, Atomdaten)] in Dateireihenfolge
    groups = {}
    modres = {}   # Residuenname -> Ausgangs-Aminosäure (MODRES)
    ssbonds = []  # Paare ((chain, seq, icode), (chain, seq, icode)) aus SSBOND
    conect = []   # Atomnummern-Paare aus CONECT
    serials = {}  # Atomnummer -> ((chain, seq, icode), Residuenname, Atomname), erstes Vorkommen
    current_model = None
    with open(file_path, 'r') as f:
        for line in f:
            record = line[:6].strip()
            if record == 'MODRES':
                modres.setdefault(line[12:15].strip(), line[24:27].strip())
                continue
            if record == 'SSBOND':
                # Brücken zu Symmetrie-Kopien (Symmetrie-Operator ≠ 1555) liegen nicht in der Datei
                sym1, sym2 = line[59:65].strip(), line[66:72].strip()
                seq1, seq2 = _int_or_none(line[17:21]), _int_or_none(line[31:35])
                if seq1 is not None and seq2 is not None and sym1 in ('', '1555') and sym2 in ('', '1555'):
                    ssbonds.append(((line[15], seq1, line[21].strip()), (line[29], seq2, line[35].strip())))
                continue
            if record == 'CONECT':
                numbers = [_int_or_none(line[i:i + 5]) for i in range(6, 31, 5)]
                if numbers[0] is not None:
                    conect.extend((numbers[0], n) for n in numbers[1:] if n is not None)
                continue
            if record == 'MODEL':
                try:
                    current_model = int(line[10:14].strip())
                except ValueError:
                    current_model = (current_model or 0) + 1
                continue
            if record == 'ENDMDL':
                current_model = None
                continue
            if record not in ('ATOM', 'HETATM'):
                continue
            if len(line) < 54:
                continue
            try:
                atom_name = line[12:16].strip()
                alt_loc = line[16]
                res_name = line[17:20].strip()
                chain_id = line[21]
                res_seq = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except (ValueError, IndexError):
                continue

            alt_loc = alt_loc.strip()
            i_code = line[26].strip() if len(line) > 26 else ''
            # CONECT bezieht sich auf die Atomnummern des ersten Models
            serials.setdefault(_int_or_none(line[6:11]), ((chain_id, res_seq, i_code), res_name, atom_name))
            if model_id is not None and current_model != model_id:
                continue
            element = line[76:78].strip() if len(line) > 76 else ''
            if not element:
                element = ''.join(c for c in atom_name if c.isalpha())[:1]
            groups.setdefault((chain_id, res_seq, i_code), []).append((atom_name, alt_loc, res_name, {
                'x': x, 'y': y, 'z': z,
                'element': element,
                'record': record,
                'serial': _int_or_none(line[6:11]),
                'alt_loc': alt_loc,
                'i_code': i_code,
                'occupancy': _float_or_none(line[54:60]),
                'temp_factor': _float_or_none(line[60:66]),
                'charge': line[78:80].strip() if len(line) > 78 else '',
            }))

    residues = {}
    for (chain_id, res_seq, i_code), atoms in groups.items():
        # Verschiedene Residuennamen ohne altLoc sind keine Mikroheterogenität, sondern eine
        # doppelt vergebene Nummer (z. B. Ligand mit der Nummer eines Protein-Residuums):
        # dann jeden Namen als eigenes Residuum behandeln.
        if len({a[2] for a in atoms if not a[1]}) > 1:
            parts = {}
            for atom in atoms:
                parts.setdefault(atom[2], []).append(atom)
            parts = parts.values()
        else:
            parts = [atoms]
        for part in parts:
            res_name, residue = _select_conformation(part, altloc)
            residues[(chain_id, res_seq, i_code, res_name)] = residue

    by_position = {}
    for key, residue in residues.items():
        by_position.setdefault(key[:3], []).append(residue)
    for a, b in ssbonds:
        for residue in by_position.get(a, []):
            residue.ssbond_partners.add(b)
        for residue in by_position.get(b, []):
            residue.ssbond_partners.add(a)

    topology = _conect_topology(conect, serials)
    for (_, _, _, res_name), residue in residues.items():
        if res_name in SIDECHAIN_BONDS or not {'N', 'CA', 'C'} <= residue.keys():
            continue  # Standard-Aminosäure oder keine Aminosäure (Ligand, Wasser, Ion)
        residue.parent = modres.get(res_name) or MODIFIED_PARENTS.get(res_name)
        heavy = [name for name, atom in residue.items() if atom['element'].upper() not in ('H', 'D')]
        if res_name in topology:
            bonds = [(a, b) for a, b in topology[res_name] if a in heavy and b in heavy]
        else:
            bonds = _bonds_by_distance(residue, heavy)
        backbone = {frozenset(pair) for pair in BACKBONE_BONDS}
        residue.sidechain_bonds = [pair for pair in bonds if frozenset(pair) not in backbone]
    return residues


def _conect_topology(conect, serials):
    """Bindungen innerhalb von Residuen aus den CONECT-Einträgen, je Residuenname als
    Atomnamen-Paare (in Reihenfolge der Atome in der Datei). Dateien von RCSB enthalten
    sie für alle HETATM-Gruppen, also auch für modifizierte Aminosäuren."""
    topology = {}
    for a, b in conect:
        if a >= b or a not in serials or b not in serials:
            continue  # jede Bindung steht in beiden Richtungen in CONECT
        (res_a, name_a, atom_a), (res_b, name_b, atom_b) = serials[a], serials[b]
        if res_a == res_b and name_a == name_b:
            pairs = topology.setdefault(name_a, {})
            pairs.setdefault((atom_a, atom_b), None)
    return {name: list(pairs) for name, pairs in topology.items()}


def _bonds_by_distance(residue, names):
    """Bindungen zwischen den Atomen `names` eines Residuums über Kovalenzradien."""
    def radius(name):
        return COVALENT_RADII.get(residue[name]['element'].upper(), DEFAULT_COVALENT_RADIUS)

    return [
        (a, b)
        for i, a in enumerate(names)
        for b in names[i + 1:]
        if _distance(residue[a], residue[b]) <= radius(a) + radius(b) + BOND_TOLERANCE
    ]


def residue_pos(res_seq, i_code):
    """Position eines Residuums wie in der PDB-Datei, z. B. "52" oder "52A"."""
    return f'{res_seq}{i_code}'


def atom_id(chain_id, res_seq, i_code, atom_name):
    """Eindeutige Atom-ID, z. B. "A_52A_CA" (genutzt als Knoten-ID im Graphen)."""
    return f'{chain_id}_{residue_pos(res_seq, i_code)}_{atom_name}'


def _distance(a1, a2):
    return math.sqrt(
        (a1['x'] - a2['x']) ** 2 +
        (a1['y'] - a2['y']) ** 2 +
        (a1['z'] - a2['z']) ** 2
    )


def _disulfide_bonds(residues):
    """Disulfidbrücken als (source_id, target_id, ungerundete Länge, Herkunft) mit Herkunft
    "ssbond" (SSBOND-Eintrag) oder "distance" (SG–SG ≤ MAX_DISULFIDE_LENGTH). Beteiligt sind
    Cysteine und davon abgeleitete modifizierte Aminosäuren mit einem SG-Atom."""
    cysteines = [
        (key, residue) for key, residue in residues.items()
        if 'SG' in residue and (key[3] == 'CYS' or getattr(residue, 'parent', None) == 'CYS')
    ]
    by_position = {key[:3]: (i, key) for i, (key, _) in enumerate(cysteines)}
    found = {}
    for i, (key, residue) in enumerate(cysteines):
        for partner in getattr(residue, 'ssbond_partners', ()):
            if partner in by_position and partner != key[:3]:
                pair = tuple(sorted((i, by_position[partner][0])))
                found.setdefault(pair, 'ssbond')
    for i, (_, a) in enumerate(cysteines):
        for j in range(i + 1, len(cysteines)):
            if (i, j) not in found and _distance(a['SG'], cysteines[j][1]['SG']) <= MAX_DISULFIDE_LENGTH:
                found[(i, j)] = 'distance'
    bonds = []
    for (i, j), via in sorted(found.items()):
        (key_a, res_a), (key_b, res_b) = cysteines[i], cysteines[j]
        bonds.append((atom_id(*key_a[:3], 'SG'), atom_id(*key_b[:3], 'SG'), _distance(res_a['SG'], res_b['SG']), via))
    return bonds


def _residue_class(residues, key):
    """Name der Aminosäure für Referenzwerte: der eigene oder bei modifizierten die Ausgangs-Aminosäure."""
    return getattr(residues[key], 'parent', None) or key[3]


def _iter_bonds(residues):
    """Liefert alle Bindungen als (source_id, target_id, typ, ungerundete Länge, Referenz)."""
    ordered_keys = list(residues.keys())
    for i, key in enumerate(ordered_keys):
        chain_id, res_seq, i_code, res_name = key
        atoms = residues[key]
        res_class = _residue_class(residues, key)

        sidechain = SIDECHAIN_BONDS.get(res_name) or getattr(atoms, 'sidechain_bonds', [])
        for bond_type, bonds in (('backbone', BACKBONE_BONDS), ('sidechain', sidechain)):
            for a1n, a2n in bonds:
                if a1n in atoms and a2n in atoms:
                    yield (atom_id(chain_id, res_seq, i_code, a1n), atom_id(chain_id, res_seq, i_code, a2n),
                           bond_type, _distance(atoms[a1n], atoms[a2n]),
                           _reference(a1n, a2n, bond_type, res_name, res_class))

        if i + 1 < len(ordered_keys):
            next_key = ordered_keys[i + 1]
            next_chain, next_seq, next_i_code, _ = next_key
            if next_chain == chain_id and 'C' in atoms and 'N' in residues[next_key]:
                length = _distance(atoms['C'], residues[next_key]['N'])
                if length <= MAX_PEPTIDE_BOND_LENGTH:
                    yield (atom_id(chain_id, res_seq, i_code, 'C'), atom_id(next_chain, next_seq, next_i_code, 'N'),
                           'peptide', length,
                           _reference('C', 'N', 'peptide', res_name, res_class, _residue_class(residues, next_key)))

    for source, target, length, _ in _disulfide_bonds(residues):
        yield source, target, 'disulfide', length, None


def _dihedral(p0, p1, p2, p3):
    """Torsionswinkel p0–p1–p2–p3 in Grad (−180 … 180); hier zählt nur der Betrag (cis/trans)."""
    sub = lambda a, b: [a['x'] - b['x'], a['y'] - b['y'], a['z'] - b['z']]  # noqa: E731
    cross = lambda u, v: [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]  # noqa: E731
    dot = lambda u, v: sum(a * b for a, b in zip(u, v))  # noqa: E731
    b0, b1, b2 = sub(p1, p0), sub(p2, p1), sub(p3, p2)
    n1, n2 = cross(b0, b1), cross(b1, b2)
    m = cross(n1, [c / math.sqrt(dot(b1, b1)) for c in b1])
    return math.degrees(math.atan2(dot(m, n2), dot(n1, n2)))


def _angle_reference(residues, keys, ends, center):
    """Referenz des Bindungswinkels ends[0]–center–ends[1] oder None. Atome sind
    (Residuum-Index, Atomname); `keys` sind die Residuen-Schlüssel in Dateireihenfolge."""
    (ia, a), (ic, c) = sorted(ends)  # Atom im früheren Residuum zuerst
    i_center, name = center
    if ia == ic == i_center:
        key = keys[i_center]
        if {a, name, c} - BACKBONE_ATOMS:
            return SIDECHAIN_ANGLE_REFERENCES.get((key[3], name, frozenset((a, c))))
        res_class = _residue_class(residues, key)
        if res_class not in SIDECHAIN_BONDS:
            return None
        atoms = (a, name, c) if (a, name, c) in BACKBONE_ANGLE_REFERENCES else (c, name, a)
        return _backbone_reference(BACKBONE_ANGLE_REFERENCES, atoms, _classes(own=res_class))
    # Über die Peptidbindung: CA–C–N und O–C–N (Mitte C) oder C–N–CA (Mitte N)
    prev_key, next_key = keys[ia], keys[ic]
    prev_class, next_class = _residue_class(residues, prev_key), _residue_class(residues, next_key)
    if prev_class not in SIDECHAIN_BONDS or next_class not in SIDECHAIN_BONDS:
        return None
    prev_atoms, next_atoms = residues[prev_key], residues[next_key]
    cis = (all(n in prev_atoms for n in ('CA', 'C')) and all(n in next_atoms for n in ('N', 'CA'))
           and abs(_dihedral(prev_atoms['CA'], prev_atoms['C'], next_atoms['N'], next_atoms['CA'])) < CIS_MAX_OMEGA)
    return _backbone_reference(BACKBONE_ANGLE_REFERENCES, (a, name, c),
                               _classes(prev=prev_class, next_=next_class, cis=cis))


def _angle_references(residues, bonds):
    """Referenzen der Bindungswinkel je Atom in der Mitte: {Atom-ID: {"Ende|Ende": Referenz}}
    mit den alphabetisch sortierten Namen der Endatome. Ohne Disulfidbrücken."""
    keys = list(residues.keys())
    where = {
        atom_id(*key[:3], name): (i, name)
        for i, key in enumerate(keys) for name in residues[key]
    }
    neighbors = {}
    for source, target, bond_type, _, _ in bonds:
        if bond_type != 'disulfide':
            neighbors.setdefault(source, []).append(target)
            neighbors.setdefault(target, []).append(source)
    references = {}
    for center_id, ids in neighbors.items():
        for i, first in enumerate(ids):
            for second in ids[i + 1:]:
                ref = _angle_reference(residues, keys, (where[first], where[second]), where[center_id])
                if ref:
                    label = '|'.join(sorted((where[first][1], where[second][1])))
                    references.setdefault(center_id, {})[label] = ref
    return references


def build_graph(residues):
    nodes = []
    # Reihenfolge, in der die Residuen in der PDB-Datei auftreten, beibehalten
    # (nicht nach Chain/Sequenznummer neu sortieren) — HETATM-Einträge können
    # eine Chain-ID wiederverwenden, obwohl sie physisch nach anderen Chains
    # in der Datei stehen.
    ordered_keys = list(residues.keys())

    for i, key in enumerate(ordered_keys):
        chain_id, res_seq, i_code, res_name = key
        atoms = residues[key]

        for atom_name, atom_data in atoms.items():
            nodes.append({
                'id': atom_id(chain_id, res_seq, i_code, atom_name),
                'label': atom_name,
                'residue': res_name,
                'residue_seq': res_seq,
                'chain': chain_id,
                'residue_index': i,
                'element': atom_data.get('element', ''),
                'hetatm': atom_data.get('record') == 'HETATM',
                # Rohdaten aus der ATOM/HETATM-Zeile (für die Datentabelle)
                'record': atom_data.get('record'),
                'serial': atom_data.get('serial'),
                'alt_loc': atom_data.get('alt_loc', ''),
                'i_code': atom_data.get('i_code', ''),
                'x': atom_data['x'],
                'y': atom_data['y'],
                'z': atom_data['z'],
                'occupancy': atom_data.get('occupancy'),
                'temp_factor': atom_data.get('temp_factor'),
                'charge': atom_data.get('charge', ''),
            })

    bonds = list(_iter_bonds(residues))
    edges = [
        # ref: [Mittelwert, Standardabweichung, Gruppe] nach REFERENCE_SOURCES oder None
        {'source': source, 'target': target, 'length': round(length, 3), 'type': bond_type, 'ref': ref}
        for source, target, bond_type, length, ref in bonds
    ]
    # Referenzen der Bindungswinkel am Atom in der Mitte, siehe _angle_references
    angle_refs = _angle_references(residues, bonds)
    for node in nodes:
        if node['id'] in angle_refs:
            node['angle_refs'] = angle_refs[node['id']]
    # Herkunft der Disulfidbrücken: "ssbond" oder "distance"
    origin = {(source, target): via for source, target, _, via in _disulfide_bonds(residues)}
    for edge in edges:
        if edge['type'] == 'disulfide':
            edge['via'] = origin[(edge['source'], edge['target'])]

    residue_list = [
        {
            'index': i,
            'name': res_name,
            'one_letter': THREE_TO_ONE.get(res_name) or THREE_TO_ONE.get(getattr(residues[key], 'parent', None), '?'),
            # Ausgangs-Aminosäure, falls es eine modifizierte ist (z. B. "MET" für MSE)
            'parent': getattr(residues[key], 'parent', None),
            'seq': res_seq,
            'i_code': i_code,
            # Position wie in der Datei ("52A"); Teil der Atom-IDs
            'pos': residue_pos(res_seq, i_code),
            'chain': chain_id,
            # Konformationen des Residuums in der Datei mit ihrer Aminosäure, z. B. {"A": "PRO", "B": "SER"}
            'conformations': getattr(residues[key], 'conformations', {}),
        }
        for i, key in enumerate(ordered_keys)
        for chain_id, res_seq, i_code, res_name in [key]
    ]

    return {'nodes': nodes, 'edges': edges, 'residues': residue_list}


# Über die Models gemittelte Atomfelder mit ihrer Anzahl Nachkommastellen (wie im PDB-Format)
AVERAGED_ATOM_FIELDS = {'x': 3, 'y': 3, 'z': 3, 'occupancy': 2, 'temp_factor': 2}


def build_average_graph(file_path, models, altloc=None):
    """Baut den Graphen eines Ensembles, dessen Bindungslängen der Mittelwert der
    jeweiligen Bindung über alle angegebenen Models sind.

    Atome und Bindungen stammen aus dem ersten Model; gemittelt wird je Bindung über
    alle Models, in denen beide Atome vorkommen. Gemittelt werden die ungerundeten
    Abstände, gerundet wird erst das Ergebnis. Ebenso werden Koordinaten, Besetzung
    und B-Faktor jedes Atoms über die Models gemittelt, in denen es vorkommt.
    """
    if not models:
        return build_graph(parse_pdb(file_path, altloc=altloc))

    base = build_graph(parse_pdb(file_path, model_id=models[0], altloc=altloc))
    sums = [0.0] * len(base['edges'])
    counts = [0] * len(base['edges'])
    atom_sums = {node['id']: {f: [0.0, 0] for f in AVERAGED_ATOM_FIELDS} for node in base['nodes']}
    for model_id in models:
        residues = parse_pdb(file_path, model_id=model_id, altloc=altloc)
        atoms = {
            atom_id(chain_id, res_seq, i_code, atom_name): atom
            for (chain_id, res_seq, i_code, _), res_atoms in residues.items()
            for atom_name, atom in res_atoms.items()
        }
        for i, edge in enumerate(base['edges']):
            a1, a2 = atoms.get(edge['source']), atoms.get(edge['target'])
            if a1 and a2:
                length = _distance(a1, a2)
                # Wie in _iter_bonds: In diesem Model ist hier eine Kettenlücke, keine Bindung
                if edge['type'] == 'peptide' and length > MAX_PEPTIDE_BOND_LENGTH:
                    continue
                sums[i] += length
                counts[i] += 1
        for node_id, fields in atom_sums.items():
            atom = atoms.get(node_id)
            if not atom:
                continue
            for field, acc in fields.items():
                if atom.get(field) is not None:
                    acc[0] += atom[field]
                    acc[1] += 1

    for i, edge in enumerate(base['edges']):
        if counts[i]:
            edge['length'] = round(sums[i] / counts[i], 3)
    for node in base['nodes']:
        for field, (total, count) in atom_sums[node['id']].items():
            if count:
                node[field] = round(total / count, AVERAGED_ATOM_FIELDS[field])
    return base


def _bond_atom(residues, node_id):
    """Angaben zu einem Atom einer Bindung (gleiche Schlüssel wie die Graph-Knoten)."""
    chain_id, pos, atom_name = node_id.split('_', 2)
    for (c, seq, i_code, res_name), atoms in residues.items():
        if c == chain_id and residue_pos(seq, i_code) == pos and atom_name in atoms:
            atom = atoms[atom_name]
            return {
                'label': atom_name, 'serial': atom.get('serial'), 'residue': res_name,
                'residue_seq': seq, 'i_code': atom.get('i_code', ''), 'chain': chain_id,
            }
    return None


def build_bond_ensemble(file_path, models, altloc=None):
    """Bindungslängen über alle Models: Mittelwert und Median je Bindung sowie die
    kürzeste und längste Bindung aller Models (inkl. beteiligter Atome und Model).

    Mittelwert und Median werden aus den ungerundeten Abständen berechnet, damit sie
    mit dem Model "Durchschnitt" übereinstimmen. Dateien ohne MODEL-Records werden
    als ein einzelnes Model behandelt.
    """
    lengths = {}
    extremes = {'min': None, 'max': None}
    for model_id in models or [None]:
        residues = parse_pdb(file_path, model_id=model_id, altloc=altloc)
        model_min = model_max = None
        for source, target, _, length, _ in _iter_bonds(residues):
            lengths.setdefault(f'{source}|{target}', []).append(length)
            if model_min is None or length < model_min[2]:
                model_min = (source, target, length)
            if model_max is None or length > model_max[2]:
                model_max = (source, target, length)
        for kind, candidate, better in (('min', model_min, lambda a, b: a < b),
                                        ('max', model_max, lambda a, b: a > b)):
            if candidate and (extremes[kind] is None or better(candidate[2], extremes[kind]['_raw'])):
                source, target, length = candidate
                extremes[kind] = {
                    'source': source, 'target': target, 'length': round(length, 3), '_raw': length,
                    'model': model_id,
                    'atoms': [_bond_atom(residues, source), _bond_atom(residues, target)],
                }
    for extreme in extremes.values():
        if extreme:
            del extreme['_raw']

    return {
        'models': models,
        # Bindung "source|target" -> [Mittelwert, Median, Anzahl Models]
        'bonds': {
            key: [round(statistics.fmean(v), 3), round(statistics.median(v), 3), len(v)]
            for key, v in lengths.items()
        },
        'min': extremes['min'],
        'max': extremes['max'],
    }
