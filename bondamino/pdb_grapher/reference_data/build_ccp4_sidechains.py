"""Builds the side-chain reference values from the CCP4 Monomer Library as compact JSON.

One-off maintenance script, not used at runtime; `utils.py` only reads its output.

Usage: python -I build_ccp4_sidechains.py <folder with the .cif files> ccp4_sidechains.json

The input folder must contain one `<RES>.cif` file per residue in `RESIDUES`, taken from
the monomer library repository at commit `COMMIT` (see README.md). Bond lengths and bond
angles are kept only if they involve no hydrogen atom and at least one atom outside the
backbone (N, CA, C, O, OXT); backbone values come from Engh & Huber instead, because the
library describes the amino acids as free molecules.

Output format (``residues`` part)::

    {"ALA": {"bonds": [[atom1, atom2, mean_A, esd_A], ...],
             "angles": [[atom1, center_atom, atom3, mean_deg, esd_deg], ...]}, ...}

plus `_source`, `_commit`, `_license` and `_content` keys describing provenance.
"""
import json
import re
import sys
from pathlib import Path

# Commit of https://github.com/MonomerLibrary/monomers the .cif files were taken from.
COMMIT = '1c97ea23a7df06373ed01cd7c3156a1f8e9b27e3'
# The 20 standard amino acids followed by the modified amino acids of utils.MODIFIED_PARENTS.
RESIDUES = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
            'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            'MSE', 'FME', 'SEP', 'TPO', 'PTR', 'CSO', 'CME', 'OCS', 'CSD', 'HYP', 'MLY', 'M3L', 'ALY', 'KCX',
            'LLP', 'CGU', 'NEP']
# CIF token: enclosed in ' or " (the quote character may appear inside the value if it is
# not followed by whitespace), or unquoted. Exactly one of the three groups matches.
TOKEN = re.compile(r"'((?:[^']|'(?=\S))*)'(?=\s|$)|\"((?:[^\"]|\"(?=\S))*)\"(?=\s|$)|(\S+)")
# Backbone atom names; restraints between these atoms only are skipped.
BACKBONE = {'N', 'CA', 'C', 'O', 'OXT'}


def loops(path):
    """Parses the loop_ tables of a monomer library CIF file.

    A minimal CIF reader sufficient for the monomer library: it reads every loop_ in the
    data blocks other than `data_comp_list` (i.e. in `data_comp_<RES>`), including
    multi-line text fields delimited by lines starting with ';'. Single key-value items
    outside loops are ignored.

    Args:
        path: `pathlib.Path` of the .cif file.

    Returns:
        A dict {category: [row, ...]}, e.g. {'_chem_comp_bond': [{'atom_id_1': 'N', ...}]},
        where each row maps the item names (without category prefix) to string values.
    """
    tables, lines = {}, path.read_text(encoding='utf-8').splitlines()
    in_comp, i = False, 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('data_'):
            in_comp = line != 'data_comp_list'
        if in_comp and line == 'loop_':
            i += 1
            fields = []
            while i < len(lines) and lines[i].strip().startswith('_'):
                fields.append(lines[i].strip())
                i += 1
            category = fields[0].split('.')[0]
            names = [f.split('.')[1] for f in fields]
            rows, tokens = [], []
            while i < len(lines) and lines[i].strip() and not lines[i].strip().startswith(('loop_', '_', 'data_', '#')):
                if lines[i].startswith(';'):  # Multi-line text field up to the next ';' line.
                    text = [lines[i][1:]]
                    i += 1
                    while not lines[i].startswith(';'):
                        text.append(lines[i])
                        i += 1
                    tokens.append('\n'.join(text))
                else:
                    tokens += [next(g for g in m.groups() if g is not None) for m in TOKEN.finditer(lines[i])]
                # A row may span several lines, so collect tokens until a full row is available.
                while len(tokens) >= len(names):
                    rows.append(dict(zip(names, tokens[:len(names)])))
                    tokens = tokens[len(names):]
                i += 1
            tables[category] = rows
            continue
        i += 1
    return tables


def main(src, out):
    """Extracts the side-chain restraints and writes them as JSON.

    Args:
        src: Folder containing `<RES>.cif` for every residue in `RESIDUES`.
        out: Path of the JSON file to write (normally `ccp4_sidechains.json`).
    """
    data = {}
    for res in RESIDUES:
        t = loops(Path(src) / f'{res}.cif')
        hydrogen = {a['atom_id'] for a in t['_chem_comp_atom'] if a['type_symbol'].upper() in ('H', 'D')}
        # Keep a restraint if it has no hydrogen and at least one non-backbone atom.
        keep = lambda atoms: not hydrogen & set(atoms) and set(atoms) - BACKBONE  # noqa: E731
        # Lengths and their esd are rounded to 0.001 Å, angles to 0.1° and their esd to 0.01°.
        bonds = [[b['atom_id_1'], b['atom_id_2'], round(float(b['value_dist']), 3), round(float(b['value_dist_esd']), 3)]
                 for b in t['_chem_comp_bond'] if keep((b['atom_id_1'], b['atom_id_2']))]
        angles = [[a['atom_id_1'], a['atom_id_2'], a['atom_id_3'], round(float(a['value_angle']), 1),
                   round(float(a['value_angle_esd']), 2)]
                  for a in t.get('_chem_comp_angle', []) if keep((a['atom_id_1'], a['atom_id_2'], a['atom_id_3']))]
        data[res] = {'bonds': bonds, 'angles': angles}
    doc = {
        '_source': 'CCP4 Monomer Library, https://github.com/MonomerLibrary/monomers',
        '_commit': COMMIT,
        '_license': 'LGPL-3.0, see LICENSE-monomers.txt and GPL-3.0.txt',
        '_content': 'value_dist/value_dist_esd (Å) and value_angle/value_angle_esd (°) without hydrogens, '
                    'only restraints involving at least one atom outside N, CA, C, O, OXT',
        'residues': data,
    }
    # Hand-formatted JSON with one line per bond or angle, so that changes stay readable in diffs.
    dump = lambda v: json.dumps(v, ensure_ascii=False)  # noqa: E731
    lines = ['{'] + [f' {dump(k)}: {dump(v)},' for k, v in doc.items() if k != 'residues'] + [' "residues": {']
    for n, (res, entry) in enumerate(data.items()):
        lines.append(f'  {dump(res)}: {{')
        for key in ('bonds', 'angles'):
            rows = ',\n'.join(f'   {dump(r)}' for r in entry[key])
            lines.append(f'   "{key}": [\n{rows}\n   ]' + (',' if key == 'bonds' else '') if rows else f'   "{key}": []' + (',' if key == 'bonds' else ''))
        lines.append('  }' + (',' if n < len(data) - 1 else ''))
    lines += [' }', '}']
    text = '\n'.join(lines) + '\n'
    assert json.loads(text) == doc  # The manual formatting must yield exactly `doc`.
    Path(out).write_text(text, encoding='utf-8')
    print(sum(len(v['bonds']) for v in data.values()), 'bonds,', sum(len(v['angles']) for v in data.values()), 'angles')


if __name__ == '__main__':
    main(*sys.argv[1:])
