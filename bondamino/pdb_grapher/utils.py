"""Scientific core of Bondamino: PDB parsing, bond graph, reference values and statistics.

This module is pure Python (no Django imports) and turns a PDB file into the JSON-ready
data the front-end draws as interactive graphs. All lengths are in Å, all angles in degrees.

Typical flow (see views.py)::

    info = parse_header(path)                           # method, resolution, R values
    models = list_models(path)                          # [] if the file has no MODEL records
    residues = parse_pdb(path, model_id=models[0] if models else None, altloc=None)
    graph = build_graph(residues)                       # nodes, edges, residues (+ Ramachandran)
    graph = build_average_graph(path, models)           # alternative: mean over all models
    ensemble = build_bond_ensemble(path, models)        # per-bond mean/median, min/max bond

Main entry points:
    parse_header, list_models, list_altlocs: Cheap scans of a file for metadata.
    parse_pdb: Reads ATOM/HETATM records of one model into `Residue` objects, choosing a
        single alternate conformation per residue and deriving bonds of ligands and
        modified amino acids (from CONECT records or covalent radii).
    build_graph: Builds the graph dict for the front-end from `parse_pdb()` output.
    build_average_graph: Same, with bond lengths and atom fields averaged over models.
    build_bond_ensemble: Bond length statistics across all models.
    evaluate_ramachandran, ramachandran_value, ramachandran_region: Top8000-based
        Ramachandran scoring as in MolProbity.

Graph dict produced by `build_graph()` / `build_average_graph()`::

    {
      'nodes': [{                      # one per atom, in file order
          'id': 'A_52A_CA',            # see atom_id(); referenced by edges
          'label': 'CA', 'residue': 'ALA', 'residue_seq': 52, 'chain': 'A',
          'residue_index': 0,          # index into 'residues'
          'element': 'C', 'hetatm': False,
          'record', 'serial', 'alt_loc', 'i_code', 'x', 'y', 'z', 'occupancy',
          'temp_factor', 'charge',     # raw ATOM/HETATM fields for the data table
          'angle_refs': {'C|N': [mean, sigma, group], ...},   # optional, only for atoms
      }, ...],                         # that are the centre of an angle with a reference
      'edges': [{                      # one per bond
          'source': node id, 'target': node id,
          'length': 1.458,             # Å, rounded to 3 decimals
          'type': 'backbone' | 'sidechain' | 'peptide' | 'disulfide' | 'ligand',
          'ref': [mean, sigma, group] | None,   # reference value, see REFERENCE_SOURCES
          'via': 'ssbond' | 'distance',          # only for disulfide bonds
      }, ...],
      'residues': [{                   # one per residue, in file order
          'index', 'name', 'one_letter', 'parent', 'seq', 'i_code', 'pos', 'chain',
          'conformations': {'A': 'PRO', 'B': 'SER'},
          'rama': {'phi', 'psi', 'type', 'value', 'region'} | None,
      }, ...],
    }

Reference values are lists `[mean, sigma, group]`, where `group` is a string such as
'backbone:N-CA:GLY' or 'SER:CB-OG' that identifies which reference applied (used to group
statistics). Backbone values follow Engh & Huber (2001), side-chain values the CCP4
Monomer Library, Ramachandran regions the Top8000 data set; the data files and their
licenses are described in reference_data/README.md.
"""
import gzip
import json
import math
import statistics
from pathlib import Path

# Bonds within the backbone of every amino acid, as atom name pairs. OXT only exists at
# the C terminus; pairs whose atoms are missing are skipped.
BACKBONE_BONDS = [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('C', 'OXT')]

# Longest C–N distance (Å) still treated as a peptide bond (ideal ≈ 1.33 Å). Larger
# distances between consecutive residues are chain breaks, e.g. caused by residues
# missing from the model, not bonds.
MAX_PEPTIDE_BOND_LENGTH = 2.0

# Longest SG–SG distance (Å) treated as a disulfide bond without an SSBOND record (typically
# ≈ 2.0 Å). With an SSBOND record the bond is drawn regardless of its length.
MAX_DISULFIDE_LENGTH = 2.5

# Side-chain bonds (heavy atoms only, PDB atom names) of the 20 standard amino acids.
# CA–CB is counted as side chain; for PRO the ring closure CD–N is included as well.
# Modified amino acids get their side-chain bonds at parse time instead (Residue.sidechain_bonds).
# Membership in this dict is also used throughout as the test "is a standard amino acid".
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

# Reference values (mean, standard deviation) for bond lengths (Å) and bond angles (°).
# Backbone: Engh & Huber (2001) with the assignment rules from Dangle (Richardson Lab,
# BSD-3-Clause); side chains: CCP4 Monomer Library (LGPL-3.0). Provenance and licenses
# are documented in reference_data/README.md. The citations below are shown in the UI.
REFERENCE_SOURCES = {
    'backbone': 'Engh & Huber (2001), International Tables for Crystallography F, 382–392',
    'sidechain': 'CCP4 Monomer Library',
}
REFERENCE_DATA_DIR = Path(__file__).resolve().parent / 'reference_data'
# Atom names that belong to the backbone; a bond or angle involving any other atom is
# looked up in the side-chain references.
BACKBONE_ATOMS = {'N', 'CA', 'C', 'O', 'OXT'}

# Backbone references per atom tuple: [(class, mean, sigma), ...]. Classes describe the
# residue itself ('GLY', 'PRO'), across the peptide bond the following residue ('X-GLY',
# 'X-PRO', 'cis-X-PRO' for a cis peptide bond) or the preceding one ('GLY-X'); 'other'
# applies otherwise. The first class in the list that applies wins (precedence as in Dangle).
BACKBONE_BOND_REFERENCES = {
    ('N', 'CA'): [('GLY', 1.456, 0.015), ('PRO', 1.468, 0.017), ('other', 1.459, 0.020)],
    ('CA', 'C'): [('GLY', 1.514, 0.016), ('PRO', 1.524, 0.020), ('other', 1.525, 0.026)],
    ('C', 'O'): [('GLY', 1.232, 0.016), ('PRO', 1.228, 0.020), ('other', 1.229, 0.019)],
    ('C', 'N'): [('X-GLY', 1.326, 0.018), ('X-PRO', 1.338, 0.019), ('other', 1.336, 0.023)],
}
# Angles as (end, centre, end). Angles across the peptide bond (C–N) list their atoms in
# order preceding residue -> following residue; within a residue both orders are looked up.
BACKBONE_ANGLE_REFERENCES = {
    ('N', 'CA', 'C'): [('GLY', 113.1, 2.5), ('PRO', 112.1, 2.6), ('other', 111.0, 2.7)],
    ('CA', 'C', 'O'): [('GLY', 120.6, 1.8), ('PRO', 120.2, 2.4), ('other', 120.1, 2.1)],
    ('CA', 'C', 'N'): [('GLY-X', 116.2, 2.0), ('X-PRO', 117.1, 2.8), ('other', 117.2, 2.2)],
    ('O', 'C', 'N'): [('X-PRO', 121.1, 1.9), ('GLY-X', 123.2, 1.7), ('other', 122.7, 1.6)],
    ('C', 'N', 'CA'): [('cis-X-PRO', 127.0, 2.4), ('X-PRO', 119.3, 1.5), ('X-GLY', 122.3, 2.1), ('other', 121.7, 2.5)],
}
CIS_MAX_OMEGA = 30  # |ω| < 30°: cis peptide bond (same threshold as in the browser).


def _load_sidechain_references():
    """Loads the side-chain references from reference_data/ccp4_sidechains.json.

    Keys use the residue's own name, so modified amino acids (e.g. MSE) have their own
    entries. Atom pairs and angle end atoms are frozensets so that the lookup does not
    depend on atom order.

    Returns:
        A tuple (bonds, angles) of dicts mapping to [mean, sigma, group]:
        bonds {(res_name, frozenset({atom1, atom2})): ...} and
        angles {(res_name, center_atom, frozenset({end1, end2})): ...}.
    """
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
    """Returns the reference classes (see BACKBONE_BOND_REFERENCES) that apply.

    Args:
        own: Amino acid of the residue itself (for bonds/angles within a residue).
        prev: Amino acid before the peptide bond.
        next_: Amino acid after the peptide bond.
        cis: Whether the peptide bond is cis.

    Returns:
        A set of class names, always including 'other'.
    """
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
    """Returns the first matching backbone reference as [mean, sigma, group], or None.

    Args:
        table: BACKBONE_BOND_REFERENCES or BACKBONE_ANGLE_REFERENCES.
        atoms: Atom name tuple in the order used as key in `table`.
        classes: Applicable classes from `_classes()`.
    """
    for cls, mean, sigma in table.get(atoms, []):
        if cls in classes:
            return [mean, sigma, f"backbone:{'-'.join(atoms)}:{cls}"]
    return None


def _reference(atom1, atom2, bond_type, res_name, res_class, next_class=None):
    """Returns the reference value of a bond.

    Side-chain bonds are looked up by the residue's own name (CCP4 has entries for the
    supported modified amino acids); backbone and peptide bonds by the class of the
    (parent) amino acid. Residues that are not amino acids get no reference.

    Args:
        atom1: Name of the first atom.
        atom2: Name of the second atom.
        bond_type: Edge type, e.g. 'backbone', 'sidechain' or 'peptide'.
        res_name: Residue name as in the file (e.g. 'MSE').
        res_class: Name of the residue's (parent) amino acid (e.g. 'MET' for MSE).
        next_class: Class of the following residue; only used for the peptide bond.

    Returns:
        [mean, sigma, group] or None if no reference applies.
    """
    if bond_type == 'peptide':
        if res_class not in SIDECHAIN_BONDS or next_class not in SIDECHAIN_BONDS:
            return None
        return _backbone_reference(BACKBONE_BOND_REFERENCES, ('C', 'N'), _classes(next_=next_class))
    if {atom1, atom2} - BACKBONE_ATOMS:
        return SIDECHAIN_BOND_REFERENCES.get((res_name, frozenset((atom1, atom2))))
    if res_class not in SIDECHAIN_BONDS:
        return None  # Not an amino acid, or a modified one with unknown parent amino acid.
    atoms = (atom1, atom2) if (atom1, atom2) in BACKBONE_BOND_REFERENCES else (atom2, atom1)
    return _backbone_reference(BACKBONE_BOND_REFERENCES, atoms, _classes(own=res_class))


def parse_header(file_path):
    """Reads structure metadata from the PDB header.

    Reading stops at the first ATOM, HETATM or MODEL record.

    Args:
        file_path: Path of the PDB file.

    Returns:
        A dict with the keys 'id_code' (HEADER), 'method' (EXPDTA; several methods are
        joined with '; '), 'resolution' (REMARK 2, Å), 'r_work' and 'r_free' (first values
        found in REMARK 3). Missing values are None.
    """
    info = {'id_code': None, 'method': None, 'resolution': None, 'r_work': None, 'r_free': None}
    with open(file_path, 'r') as f:
        for line in f:
            record = line[:6].strip()
            if record in ('ATOM', 'HETATM', 'MODEL'):
                break  # End of the header.
            if record == 'HEADER':
                info['id_code'] = line[62:66].strip() or None  # Columns 63-66: PDB ID.
            elif record == 'EXPDTA':
                method = line[10:].strip()
                info['method'] = f"{info['method']}; {method}" if info['method'] else method
            elif line.startswith('REMARK   2 RESOLUTION.'):
                value = line[22:].split()
                info['resolution'] = _float_or_none(value[0]) if value else None
            elif line.startswith('REMARK   3') and ':' in line:
                # "REMARK   3   R VALUE     (WORKING SET) : 0.190"; whitespace in labels varies.
                label, _, value = line[10:].partition(':')
                label = ' '.join(label.split())
                if label in ('R VALUE (WORKING SET)', 'R VALUE (WORKING + TEST SET)') and info['r_work'] is None:
                    info['r_work'] = _float_or_none(value.strip())
                elif label == 'FREE R VALUE' and info['r_free'] is None:
                    info['r_free'] = _float_or_none(value.strip())
    return info


# Parent amino acid of common modified amino acids; only used if the file has no MODRES
# record for the residue (files from RCSB include them). Each of these also needs side-chain
# references in ccp4_sidechains.json (see reference_data/build_ccp4_sidechains.py).
MODIFIED_PARENTS = {
    'MSE': 'MET', 'FME': 'MET',
    'SEP': 'SER', 'TPO': 'THR', 'PTR': 'TYR',
    'CSO': 'CYS', 'CME': 'CYS', 'OCS': 'CYS', 'CSD': 'CYS',
    'HYP': 'PRO',
    'MLY': 'LYS', 'M3L': 'LYS', 'ALY': 'LYS', 'KCX': 'LYS', 'LLP': 'LYS',
    'CGU': 'GLU', 'NEP': 'HIS',
}

# Covalent radii in Å (Cordero et al. 2008) per element symbol (upper case), used to
# detect bonds by distance when a residue has no CONECT records. 'CA' is calcium here,
# not the C-alpha atom name.
COVALENT_RADII = {
    'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05, 'SE': 1.20, 'P': 1.07,
    'F': 0.57, 'CL': 1.02, 'BR': 1.20, 'I': 1.39, 'B': 0.84,
    # Metals in ligands (e.g. Fe in haem, Mg in chlorophyll); Cordero et al. (2008).
    'FE': 1.32, 'ZN': 1.22, 'MG': 1.41, 'CU': 1.32, 'MN': 1.39, 'CO': 1.26, 'NI': 1.24,
    'CA': 1.76, 'MO': 1.54, 'W': 1.62, 'V': 1.53,
}
# Water gets no bonds (only one heavy atom).
WATER_NAMES = {'HOH', 'DOD', 'WAT', 'H2O'}
DEFAULT_COVALENT_RADIUS = 0.80  # Å, for elements missing from COVALENT_RADII.
# Two atoms count as bonded if their distance is at most the sum of their radii plus this
# tolerance in Å (common values are 0.4–0.45 Å).
BOND_TOLERANCE = 0.45

# One-letter codes of the standard amino acids (modified ones use their parent's code).
THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


def _float_or_none(text):
    """Parses a float from a PDB field; returns None for empty or invalid text."""
    try:
        return float(text)
    except ValueError:
        return None


def _int_or_none(text):
    """Parses an int from a PDB field; returns None for empty or invalid text."""
    try:
        return int(text)
    except ValueError:
        return None


def list_models(file_path):
    """Lists the MODEL numbers of a PDB file in order of appearance.

    Args:
        file_path: Path of the PDB file.

    Returns:
        A list of model numbers. Empty if the file has no MODEL records; the file is then
        treated as a single, implicit model (pass `model_id=None` to `parse_pdb()`).
    """
    models = []
    with open(file_path, 'r') as f:
        for line in f:
            if line[:6].strip() != 'MODEL':
                continue
            try:
                model_num = int(line[10:14].strip())  # Columns 11-14: model serial number.
            except ValueError:
                # Missing or invalid number: use the position of the model in the file.
                model_num = len(models) + 1
            models.append(model_num)
    return models


def list_altlocs(file_path):
    """Lists the alternate location (altLoc) identifiers used in a PDB file.

    Args:
        file_path: Path of the PDB file.

    Returns:
        The sorted identifiers, e.g. ['A', 'B']; empty if there are none.
    """
    altlocs = set()
    with open(file_path, 'r') as f:
        for line in f:
            # Column 17: altLoc indicator.
            if line[:6].strip() in ('ATOM', 'HETATM') and len(line) > 16 and line[16].strip():
                altlocs.add(line[16])
    return sorted(altlocs)


def _choose_altloc(atoms, preferred):
    """Chooses the alternate conformation to use for one residue.

    Args:
        atoms: Tuples (atom_name, alt_loc, res_name, atom_data) in file order.
        preferred: altLoc requested by the user, or None.

    Returns:
        `preferred` if the residue has it, otherwise the altLoc with the highest mean
        occupancy (on a tie the first one in the file); '' if the residue has no
        alternate conformations.
    """
    occupancies = {}
    for _, alt_loc, _, atom in atoms:
        if alt_loc:
            occupancies.setdefault(alt_loc, []).append(atom['occupancy'] or 0.0)
    if not occupancies:
        return ''
    if preferred in occupancies:
        return preferred
    # On a tie max() returns the first element, and the dict keeps file order.
    return max(occupancies, key=lambda a: statistics.fmean(occupancies[a]))


class Residue(dict):
    """Atoms of one residue as {atom_name: atom_data}, with extra attributes.

    `atom_data` is a dict with the keys 'x', 'y', 'z' (Å), 'element', 'record' ('ATOM' or
    'HETATM'), 'serial', 'alt_loc', 'i_code', 'occupancy', 'temp_factor' and 'charge'.

    Attributes:
        conformations: Maps every alternate conformation of the residue in the file to its
            residue name, e.g. {'A': 'PRO', 'B': 'SER'}; empty without alternatives.
        parent: For modified amino acids the parent amino acid (e.g. 'MET' for MSE), or
            None if it is unknown or the residue is a standard amino acid or ligand.
        sidechain_bonds: For modified amino acids their side-chain bonds as atom name
            pairs (all non-backbone bonds); empty otherwise.
        ligand_bonds: For ligands the bonds between their atoms as atom name pairs
            (empty for single-atom ions); None for amino acids and water.
        ssbond_partners: Positions (chain, seq, icode) of the cysteines this residue
            forms a disulfide bond with according to SSBOND records.
    """

    def __init__(self, *args, conformations=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.conformations = conformations or {}
        self.parent = None
        self.ligand_bonds = None
        self.sidechain_bonds = []
        # (chain, seq, icode) of the cysteines this residue is disulfide-bonded to per SSBOND.
        self.ssbond_partners = set()


def _select_conformation(atoms, preferred):
    """Reduces the atoms of one residue to a single conformation.

    Keeps all atoms without altLoc plus those of the chosen alternative (see
    `_choose_altloc()`); if an atom name occurs more than once, the first one wins.

    Args:
        atoms: Tuples (atom_name, alt_loc, res_name, atom_data) in file order.
        preferred: altLoc requested by the user, or None.

    Returns:
        A tuple (res_name, residue) with a `Residue`. The name is taken from the chosen
        alternative, because it differs in case of microheterogeneity (e.g. A = PRO,
        B = SER at the same position).
    """
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
    """Reads the atoms of one model from a PDB file.

    Residues are distinguished by chain, sequence number and insertion code (52 and 52A
    are different residues). Exactly one conformation is used per residue: `altloc` if
    the residue has it, otherwise the one with the highest occupancy. Besides ATOM/HETATM
    the records MODRES (parent amino acids), SSBOND (disulfide bonds) and CONECT (bonds
    within HETATM groups) are evaluated; residues that are neither standard amino acids
    nor water get their bonds from CONECT or, failing that, from covalent radii.

    Args:
        file_path: Path of the PDB file.
        model_id: MODEL number to read. If given, only atoms inside that MODEL block are
            read. Use None for files without MODEL records (otherwise no atoms are found).
            With None in a multi-model file all models are read, and for each atom its
            first occurrence (i.e. the first model) is kept.
        altloc: Preferred alternate conformation (e.g. 'A'), or None.

    Returns:
        A dict {(chain_id, res_seq, i_code, res_name): Residue} in file order.
    """
    # (chain, seq, icode) -> [(atom_name, alt_loc, res_name, atom_data)] in file order
    groups = {}
    modres = {}   # Residue name -> parent amino acid (MODRES).
    ssbonds = []  # Pairs ((chain, seq, icode), (chain, seq, icode)) from SSBOND.
    conect = []   # Pairs of atom serial numbers from CONECT.
    serials = {}  # Atom serial -> ((chain, seq, icode), res_name, atom_name), first occurrence.
    current_model = None
    with open(file_path, 'r') as f:
        for line in f:
            record = line[:6].strip()
            if record == 'MODRES':
                # Columns 13-15: residue name, 25-27: standard (parent) residue name.
                modres.setdefault(line[12:15].strip(), line[24:27].strip())
                continue
            if record == 'SSBOND':
                # Columns 16/30: chain IDs, 18-21/32-35: residue numbers, 22/36: insertion
                # codes, 60-65/67-72: symmetry operators. Bonds to symmetry copies
                # (operator other than 1555) involve atoms that are not in the file.
                sym1, sym2 = line[59:65].strip(), line[66:72].strip()
                seq1, seq2 = _int_or_none(line[17:21]), _int_or_none(line[31:35])
                if seq1 is not None and seq2 is not None and sym1 in ('', '1555') and sym2 in ('', '1555'):
                    ssbonds.append(((line[15], seq1, line[21].strip()), (line[29], seq2, line[35].strip())))
                continue
            if record == 'CONECT':
                # Columns 7-11: atom serial, 12-31: up to four bonded atom serials.
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
                continue  # Too short to contain coordinates.
            # Fixed PDB columns (1-based): 13-16 atom name, 17 altLoc, 18-20 residue name,
            # 22 chain, 23-26 residue number, 27 insertion code, 31-54 x/y/z (Å),
            # 55-60 occupancy, 61-66 B factor, 77-78 element, 79-80 charge.
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
            # CONECT refers to the atom serial numbers of the first model, so record them
            # before filtering by model.
            serials.setdefault(_int_or_none(line[6:11]), ((chain_id, res_seq, i_code), res_name, atom_name))
            if model_id is not None and current_model != model_id:
                continue
            element = line[76:78].strip() if len(line) > 76 else ''
            if not element:
                # No element column (older files): guess from the first letter of the atom name.
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
        # Different residue names without altLoc are not microheterogeneity but a number
        # used twice (e.g. a ligand with the number of a protein residue): treat each name
        # as a residue of its own.
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

    # Attach SSBOND partners to both residues (all residues at that position).
    by_position = {}
    for key, residue in residues.items():
        by_position.setdefault(key[:3], []).append(residue)
    for a, b in ssbonds:
        for residue in by_position.get(a, []):
            residue.ssbond_partners.add(b)
        for residue in by_position.get(b, []):
            residue.ssbond_partners.add(a)

    # Bonds of everything that is neither a standard amino acid nor water (hydrogens ignored).
    topology = _conect_topology(conect, serials)
    for (_, _, _, res_name), residue in residues.items():
        if res_name in SIDECHAIN_BONDS or res_name in WATER_NAMES:
            continue
        heavy = [name for name, atom in residue.items() if atom['element'].upper() not in ('H', 'D')]
        if res_name in topology:
            bonds = [(a, b) for a, b in topology[res_name] if a in heavy and b in heavy]
        else:
            bonds = _bonds_by_distance(residue, heavy)
        if {'N', 'CA', 'C'} <= residue.keys():
            # Modified amino acid: backbone as usual, everything else is its side chain.
            residue.parent = modres.get(res_name) or MODIFIED_PARENTS.get(res_name)
            backbone = {frozenset(pair) for pair in BACKBONE_BONDS}
            residue.sidechain_bonds = [pair for pair in bonds if frozenset(pair) not in backbone]
        else:
            residue.ligand_bonds = bonds  # Ligand (ions have only one atom and thus no bonds).
    return residues


def _conect_topology(conect, serials):
    """Derives intra-residue bonds from CONECT records.

    Files from RCSB contain CONECT records for all HETATM groups, i.e. also for modified
    amino acids. Bonds between different residues (e.g. to metal ions) are ignored. The
    topology is collected per residue name, so all residues of the same name share it.

    Args:
        conect: Pairs of atom serial numbers from CONECT records.
        serials: Atom serial -> ((chain, seq, icode), res_name, atom_name).

    Returns:
        A dict {res_name: [(atom_name, atom_name), ...]}, without duplicates and in order
        of first appearance in the CONECT records.
    """
    topology = {}
    for a, b in conect:
        if a >= b or a not in serials or b not in serials:
            continue  # Every bond appears in both directions in CONECT; keep one.
        (res_a, name_a, atom_a), (res_b, name_b, atom_b) = serials[a], serials[b]
        if res_a == res_b and name_a == name_b:
            pairs = topology.setdefault(name_a, {})
            pairs.setdefault((atom_a, atom_b), None)  # dict as an insertion-ordered set
    return {name: list(pairs) for name, pairs in topology.items()}


def _bonds_by_distance(residue, names):
    """Detects bonds between atoms of a residue by distance.

    Two atoms are bonded if their distance is at most the sum of their covalent radii
    plus BOND_TOLERANCE.

    Args:
        residue: The `Residue`.
        names: Names of the atoms to consider (normally the heavy atoms).

    Returns:
        A list of bonded atom name pairs.
    """
    def radius(name):
        return COVALENT_RADII.get(residue[name]['element'].upper(), DEFAULT_COVALENT_RADIUS)

    return [
        (a, b)
        for i, a in enumerate(names)
        for b in names[i + 1:]
        if _distance(residue[a], residue[b]) <= radius(a) + radius(b) + BOND_TOLERANCE
    ]


def residue_pos(res_seq, i_code):
    """Returns the residue position as written in the PDB file, e.g. "52" or "52A"."""
    return f'{res_seq}{i_code}'


def atom_id(chain_id, res_seq, i_code, atom_name):
    """Returns the unique atom ID, e.g. "A_52A_CA", used as node ID in the graph.

    The residue name is not part of the ID; `evaluate_ramachandran()` rebuilds IDs
    from chain, position and atom name.
    """
    return f'{chain_id}_{residue_pos(res_seq, i_code)}_{atom_name}'


def _distance(a1, a2):
    """Returns the Euclidean distance in Å between two atom data dicts."""
    return math.sqrt(
        (a1['x'] - a2['x']) ** 2 +
        (a1['y'] - a2['y']) ** 2 +
        (a1['z'] - a2['z']) ** 2
    )


def _disulfide_bonds(residues):
    """Finds the disulfide bonds between cysteines.

    Candidates are cysteines and modified amino acids derived from cysteine that have an
    SG atom. A pair is bonded if an SSBOND record says so or, failing that, if the SG–SG
    distance is at most MAX_DISULFIDE_LENGTH.

    Args:
        residues: Output of `parse_pdb()`.

    Returns:
        A list of (source_id, target_id, unrounded length in Å, origin) sorted by residue
        order, where origin is 'ssbond' (SSBOND record) or 'distance'.
    """
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
    """Returns the amino acid name used for references: the parent for modified residues, else its own."""
    return getattr(residues[key], 'parent', None) or key[3]


def _iter_bonds(residues):
    """Yields all bonds of a structure.

    Within each residue: backbone and side-chain bonds of amino acids, or the own bonds
    of ligands. Between consecutive residues of the same chain: the peptide bond C–N if
    it is not longer than MAX_PEPTIDE_BOND_LENGTH. Finally all disulfide bonds.

    Args:
        residues: Output of `parse_pdb()`, in file order.

    Yields:
        Tuples (source_id, target_id, type, unrounded length in Å, reference), where type
        is 'backbone', 'sidechain', 'ligand', 'peptide' or 'disulfide' and reference is
        [mean, sigma, group] or None.
    """
    ordered_keys = list(residues.keys())
    for i, key in enumerate(ordered_keys):
        chain_id, res_seq, i_code, res_name = key
        atoms = residues[key]
        res_class = _residue_class(residues, key)

        sidechain = SIDECHAIN_BONDS.get(res_name) or getattr(atoms, 'sidechain_bonds', [])
        # Ligands only get their own bonds (atoms named e.g. C and O are not backbone there).
        ligand = getattr(atoms, 'ligand_bonds', None)
        groups = [('ligand', ligand)] if ligand is not None else [('backbone', BACKBONE_BONDS), ('sidechain', sidechain)]
        for bond_type, bonds in groups:
            for a1n, a2n in bonds:
                if a1n in atoms and a2n in atoms:
                    yield (atom_id(chain_id, res_seq, i_code, a1n), atom_id(chain_id, res_seq, i_code, a2n),
                           bond_type, _distance(atoms[a1n], atoms[a2n]),
                           _reference(a1n, a2n, bond_type, res_name, res_class))

        # Peptide bond to the next residue in file order if it is in the same chain.
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
    """Returns the torsion angle p0–p1–p2–p3 in degrees.

    Range −180 … 180, sign per IUPAC convention (clockwise positive), same as in the browser.

    Args:
        p0: Atom data dict (with 'x', 'y', 'z'); likewise p1, p2, p3.
    """
    sub = lambda a, b: [a['x'] - b['x'], a['y'] - b['y'], a['z'] - b['z']]  # noqa: E731
    cross = lambda u, v: [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]  # noqa: E731
    dot = lambda u, v: sum(a * b for a, b in zip(u, v))  # noqa: E731
    b0, b1, b2 = sub(p1, p0), sub(p2, p1), sub(p3, p2)
    n1, n2 = cross(b0, b1), cross(b1, b2)
    m = cross(n1, [c / math.sqrt(dot(b1, b1)) for c in b1])
    # atan2 of the normals n1, n2 in the frame (n1, m); negated to get the IUPAC sign.
    return -math.degrees(math.atan2(dot(m, n2), dot(n1, n2)))


def _angle_reference(residues, keys, ends, center):
    """Returns the reference of the bond angle ends[0]–center–ends[1].

    Angles within one residue use the side-chain references if a non-backbone atom is
    involved, else the backbone references. Angles across the peptide bond use the
    backbone references with the classes of both residues and the cis/trans state of
    the peptide bond (ω = CA–C–N–CA).

    Args:
        residues: Output of `parse_pdb()`.
        keys: The residue keys of `residues` in file order.
        ends: The two end atoms as (residue index, atom name).
        center: The centre atom as (residue index, atom name).

    Returns:
        [mean, sigma, group] or None if no reference applies.
    """
    (ia, a), (ic, c) = sorted(ends)  # Atom in the earlier residue first.
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
    # Across the peptide bond: CA–C–N and O–C–N (centre C) or C–N–CA (centre N).
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
    """Collects the bond angle references for every centre atom.

    Every pair of bonded neighbours of an atom forms an angle. Disulfide bonds are not
    considered.

    Args:
        residues: Output of `parse_pdb()`.
        bonds: Bonds as yielded by `_iter_bonds()`.

    Returns:
        A dict {centre atom ID: {"end|end": [mean, sigma, group]}}, where "end|end" are
        the alphabetically sorted names of the two end atoms. Angles without reference
        are left out.
    """
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


# Ramachandran scoring as in MolProbity (cctbx mmtbx/validation/ramalyze.py) with the
# Top8000 grids of the Richardson Lab (CC BY 4.0, see reference_data/README.md).
RAMACHANDRAN_SOURCE = 'Top8000 (Richardson Lab, Duke University), Williams et al. (2018), Protein Science 27, 293–315'
# Region thresholds on the interpolated grid value: favored from 2 %, allowed from the
# per-category value, outlier below. Must match reference_data/build_top8000_rama.py.
RAMA_FAVORED = 0.02
RAMA_ALLOWED = {'general': 0.0005, 'glycine': 0.001, 'cis-proline': 0.002, 'trans-proline': 0.001,
                'pre-proline': 0.001, 'ile-val': 0.001}
RAMA_CIS_MAX_OMEGA = 90  # Proline counts as cis up to |ω| < 90° (twisted peptides go to the nearer form).
RAMA_BINS, RAMA_BIN_WIDTH = 180, 2.0  # Grid size per axis and bin width in degrees.
_rama_grids = None  # Cache {category: flat grid}, loaded lazily by _rama_grid().


def _rama_grid(rama_type):
    """Returns the grid values of a category (180 × 180, φ-major), loaded on first use."""
    global _rama_grids
    if _rama_grids is None:
        with gzip.open(REFERENCE_DATA_DIR / 'top8000_rama.json.gz', 'rt', encoding='utf-8') as f:
            # Stored as integers in millionths.
            _rama_grids = {name: [v / 1e6 for v in values] for name, values in json.load(f)['grids'].items()}
    return _rama_grids[rama_type]


def ramachandran_value(rama_type, phi, psi):
    """Returns the value of the Top8000 distribution at (φ, ψ).

    Bilinear interpolation between the four nearest bin centres, wrapping around at
    ±180° (as NDimTable.valueAt in MolProbity).

    Args:
        rama_type: Category, one of the keys of RAMA_ALLOWED.
        phi: φ in degrees (−180 … 180).
        psi: ψ in degrees (−180 … 180).

    Returns:
        The interpolated density value (a fraction, compared with RAMA_FAVORED/RAMA_ALLOWED).
    """
    grid = _rama_grid(rama_type)
    home, neighbor, contrib = [], [], []
    for angle in (phi, psi):
        # Bin containing the angle, the neighbouring bin towards the angle, and the
        # fractional distance from the home bin centre (0 … 0.5).
        i = min(int((angle + 180) // RAMA_BIN_WIDTH), RAMA_BINS - 1)
        center = -180 + (i + 0.5) * RAMA_BIN_WIDTH
        home.append(i)
        neighbor.append((i - 1 if angle < center else i + 1) % RAMA_BINS)
        contrib.append(abs(angle - center) / RAMA_BIN_WIDTH)
    value = 0.0
    for use_i in (False, True):
        for use_j in (False, True):
            i = neighbor[0] if use_i else home[0]
            j = neighbor[1] if use_j else home[1]
            weight = (contrib[0] if use_i else 1 - contrib[0]) * (contrib[1] if use_j else 1 - contrib[1])
            value += weight * grid[i * RAMA_BINS + j]
    return value


def ramachandran_region(rama_type, value):
    """Classifies a Ramachandran value.

    Args:
        rama_type: Category, one of the keys of RAMA_ALLOWED.
        value: Value from `ramachandran_value()`.

    Returns:
        'favored', 'allowed' or 'outlier'.
    """
    if value >= RAMA_FAVORED:
        return 'favored'
    return 'allowed' if value >= RAMA_ALLOWED[rama_type] else 'outlier'


def evaluate_ramachandran(graph):
    """Adds Ramachandran data to every residue of a graph.

    Works on the graph's nodes (not on `parse_pdb()` output), so that it also applies to
    the averaged coordinates of `build_average_graph()`. φ and ψ require peptide bonds to
    both neighbours. The category follows MolProbity: 'glycine', 'cis-proline' or
    'trans-proline', 'pre-proline' (residue before a proline), 'ile-val' or 'general';
    modified amino acids are categorised by their parent.

    Args:
        graph: Graph dict as built by `build_graph()`; modified in place.

    Returns:
        The same graph, where each residue has 'rama' set to {'phi', 'psi' (degrees,
        1 decimal), 'type' (category), 'value', 'region' ('favored', 'allowed' or
        'outlier')}, or None if φ or ψ is undefined (chain end, chain break, missing atoms).
    """
    nodes = {n['id']: n for n in graph['nodes']}
    peptides = {(e['source'], e['target']) for e in graph['edges'] if e['type'] == 'peptide'}
    residues = graph['residues']
    atom = lambda res, name: nodes.get(f"{res['chain']}_{res['pos']}_{name}")  # noqa: E731

    def linked(a, b):
        """Whether residue `a` is joined to residue `b` by a peptide bond C(a)–N(b)."""
        c, n = atom(a, 'C'), atom(b, 'N')
        return bool(c and n and (c['id'], n['id']) in peptides)

    for i, res in enumerate(residues):
        res['rama'] = None
        prev = residues[i - 1] if i > 0 else None
        next_ = residues[i + 1] if i + 1 < len(residues) else None
        n, ca, c = atom(res, 'N'), atom(res, 'CA'), atom(res, 'C')
        if not (n and ca and c and prev and next_ and linked(prev, res) and linked(res, next_)):
            continue
        phi = _dihedral(atom(prev, 'C'), n, ca, c)
        psi = _dihedral(n, ca, c, atom(next_, 'N'))
        name = res['parent'] or res['name']
        if name == 'GLY':
            rama_type = 'glycine'
        elif name == 'PRO':
            prev_ca = atom(prev, 'CA')
            cis = prev_ca is not None and abs(_dihedral(prev_ca, atom(prev, 'C'), n, ca)) < RAMA_CIS_MAX_OMEGA
            rama_type = 'cis-proline' if cis else 'trans-proline'
        elif (next_['parent'] or next_['name']) == 'PRO':
            rama_type = 'pre-proline'
        elif name in ('ILE', 'VAL'):
            rama_type = 'ile-val'
        else:
            rama_type = 'general'
        value = ramachandran_value(rama_type, phi, psi)
        res['rama'] = {
            'phi': round(phi, 1), 'psi': round(psi, 1), 'type': rama_type,
            'value': round(value, 6), 'region': ramachandran_region(rama_type, value),
        }
    return graph


def build_graph(residues):
    """Builds the graph dict consumed by the front-end.

    Args:
        residues: Output of `parse_pdb()`.

    Returns:
        A dict with 'nodes' (atoms), 'edges' (bonds with lengths and reference values)
        and 'residues' (incl. Ramachandran data); the exact shape is described in the
        module docstring.
    """
    nodes = []
    # Keep the order in which residues appear in the PDB file (do not re-sort by
    # chain/sequence number): HETATM records may reuse a chain ID although they
    # are located after other chains in the file.
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
                # Raw data from the ATOM/HETATM line (for the data table).
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
        # ref: [mean, sigma, group] according to REFERENCE_SOURCES, or None.
        {'source': source, 'target': target, 'length': round(length, 3), 'type': bond_type, 'ref': ref}
        for source, target, bond_type, length, ref in bonds
    ]
    # Bond angle references at the centre atom, see _angle_references().
    angle_refs = _angle_references(residues, bonds)
    for node in nodes:
        if node['id'] in angle_refs:
            node['angle_refs'] = angle_refs[node['id']]
    # Origin of the disulfide bonds: "ssbond" or "distance".
    origin = {(source, target): via for source, target, _, via in _disulfide_bonds(residues)}
    for edge in edges:
        if edge['type'] == 'disulfide':
            edge['via'] = origin[(edge['source'], edge['target'])]

    residue_list = [
        {
            'index': i,
            'name': res_name,
            'one_letter': THREE_TO_ONE.get(res_name) or THREE_TO_ONE.get(getattr(residues[key], 'parent', None), '?'),
            # Parent amino acid if this is a modified one (e.g. "MET" for MSE).
            'parent': getattr(residues[key], 'parent', None),
            'seq': res_seq,
            'i_code': i_code,
            # Position as in the file ("52A"); part of the atom IDs.
            'pos': residue_pos(res_seq, i_code),
            'chain': chain_id,
            # Conformations of the residue in the file with their residue name, e.g. {"A": "PRO", "B": "SER"}.
            'conformations': getattr(residues[key], 'conformations', {}),
        }
        for i, key in enumerate(ordered_keys)
        for chain_id, res_seq, i_code, res_name in [key]  # Unpacks the key within the comprehension.
    ]

    return evaluate_ramachandran({'nodes': nodes, 'edges': edges, 'residues': residue_list})


# Atom fields averaged over the models, with their number of decimals (as in the PDB format).
AVERAGED_ATOM_FIELDS = {'x': 3, 'y': 3, 'z': 3, 'occupancy': 2, 'temp_factor': 2}


def build_average_graph(file_path, models, altloc=None):
    """Builds the graph of an ensemble whose bond lengths are averaged over the models.

    Atoms and bonds are taken from the first model; each bond is averaged over all
    models in which both atoms occur (a peptide bond is skipped in models where it is a
    chain break). The unrounded distances are averaged and only the result is rounded.
    Likewise, coordinates, occupancy and B factor of each atom are averaged over the
    models in which it occurs, and φ/ψ are computed from the averaged coordinates. Bond
    angle references and residue data are those of the first model.

    Args:
        file_path: Path of the PDB file.
        models: Model numbers from `list_models()`. If empty, the graph of the whole
            file is returned unchanged.
        altloc: Preferred alternate conformation, passed to `parse_pdb()`.

    Returns:
        A graph dict as returned by `build_graph()`.
    """
    if not models:
        return build_graph(parse_pdb(file_path, altloc=altloc))

    base = build_graph(parse_pdb(file_path, model_id=models[0], altloc=altloc))
    # Running sums and counts per edge (same index as base['edges']) and per atom field.
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
                # As in _iter_bonds(): in this model there is a chain break here, not a bond.
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
    # Recompute φ/ψ from the averaged coordinates.
    return evaluate_ramachandran(base)


def _bond_atom(residues, node_id):
    """Returns details of a bond atom (same keys as the graph nodes), or None if not found.

    Args:
        residues: Output of `parse_pdb()`.
        node_id: Atom ID as built by `atom_id()`.
    """
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
    """Computes bond length statistics across all models.

    Mean and median are computed from the unrounded distances so that they agree with
    the "average" model (`build_average_graph()`). Files without MODEL records are
    treated as a single model.

    Args:
        file_path: Path of the PDB file.
        models: Model numbers from `list_models()` (may be empty).
        altloc: Preferred alternate conformation, passed to `parse_pdb()`.

    Returns:
        A dict::

            {
              'models': models,
              'bonds': {'source|target': [mean, median, number of models], ...},
              'min': {'source', 'target', 'length', 'model', 'atoms': [atom, atom]} | None,
              'max': same as 'min',
            }

        'min'/'max' are the shortest and longest bond over all models, where 'model' is
        None for files without MODEL records and 'atoms' come from `_bond_atom()`.
        Lengths are in Å, rounded to 3 decimals.
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
                    # '_raw' keeps the unrounded length for comparisons; removed below.
                    'source': source, 'target': target, 'length': round(length, 3), '_raw': length,
                    'model': model_id,
                    'atoms': [_bond_atom(residues, source), _bond_atom(residues, target)],
                }
    for extreme in extremes.values():
        if extreme:
            del extreme['_raw']

    return {
        'models': models,
        # Bond "source|target" -> [mean, median, number of models containing the bond].
        'bonds': {
            key: [round(statistics.fmean(v), 3), round(statistics.median(v), 3), len(v)]
            for key, v in lengths.items()
        },
        'min': extremes['min'],
        'max': extremes['max'],
    }
