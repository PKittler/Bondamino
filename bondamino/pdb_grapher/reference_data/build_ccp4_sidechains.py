"""Seitenketten-Referenzwerte aus der CCP4 Monomer Library als kompaktes JSON.

Aufruf: python -I build_ccp4_sidechains.py <ordner mit den .cif-dateien> ccp4_sidechains.json
Übernommen werden Bindungen und Winkel ohne Wasserstoff, an denen mindestens ein Atom
außerhalb des Rückgrats (N, CA, C, O, OXT) beteiligt ist.
"""
import json
import re
import sys
from pathlib import Path

COMMIT = '1c97ea23a7df06373ed01cd7c3156a1f8e9b27e3'
RESIDUES = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
            'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            'MSE', 'FME', 'SEP', 'TPO', 'PTR', 'CSO', 'CME', 'OCS', 'CSD', 'HYP', 'MLY', 'M3L', 'ALY', 'KCX',
            'LLP', 'CGU', 'NEP']
# CIF-Token: in ' oder " (Anführungszeichen im Wert erlaubt, wenn kein Leerraum folgt) oder ohne
TOKEN = re.compile(r"'((?:[^']|'(?=\S))*)'(?=\s|$)|\"((?:[^\"]|\"(?=\S))*)\"(?=\s|$)|(\S+)")
BACKBONE = {'N', 'CA', 'C', 'O', 'OXT'}


def loops(path):
    """Alle loop_-Blöcke des Blocks data_comp_<RES> als {Kategorie: [Zeilen als dict]}."""
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
                if lines[i].startswith(';'):  # mehrzeiliges Textfeld bis zur nächsten ';'-Zeile
                    text = [lines[i][1:]]
                    i += 1
                    while not lines[i].startswith(';'):
                        text.append(lines[i])
                        i += 1
                    tokens.append('\n'.join(text))
                else:
                    tokens += [next(g for g in m.groups() if g is not None) for m in TOKEN.finditer(lines[i])]
                while len(tokens) >= len(names):
                    rows.append(dict(zip(names, tokens[:len(names)])))
                    tokens = tokens[len(names):]
                i += 1
            tables[category] = rows
            continue
        i += 1
    return tables


def main(src, out):
    data = {}
    for res in RESIDUES:
        t = loops(Path(src) / f'{res}.cif')
        hydrogen = {a['atom_id'] for a in t['_chem_comp_atom'] if a['type_symbol'].upper() in ('H', 'D')}
        keep = lambda atoms: not hydrogen & set(atoms) and set(atoms) - BACKBONE  # noqa: E731
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
    # Eine Zeile je Bindung bzw. Winkel, damit Änderungen im Diff lesbar bleiben
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
    assert json.loads(text) == doc
    Path(out).write_text(text, encoding='utf-8')
    print(sum(len(v['bonds']) for v in data.values()), 'bonds,', sum(len(v['angles']) for v in data.values()), 'angles')


if __name__ == '__main__':
    main(*sys.argv[1:])
