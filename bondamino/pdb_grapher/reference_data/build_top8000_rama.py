"""Ramachandran-Raster (Top8000, Richardson-Labor) als kompakte Dateien für Server und Browser.

Aufruf: python -I build_top8000_rama.py <ordner mit den rama8000-*.data> top8000_rama.json.gz <regions.json>

- top8000_rama.json.gz: je Kategorie die 180 × 180 Rasterwerte (φ außen, ψ innen, Mitten
  bei −179°, −177°, … 179°) als ganze Zahlen in Millionsteln; für die Bewertung im Server.
- <regions.json>: je Kategorie eine Zeichenkette mit 0 (Ausreißer), 1 (erlaubt) oder
  2 (bevorzugt) je Rasterfeld; für die Bereiche im Ramachandran-Plot.
"""
import gzip
import json
import sys
from pathlib import Path

COMMIT = '5ee6875fc29eccc3c9dc7cbe705ff9fd1b505d7d'
FILES = {
    'general': 'rama8000-general-noGPIVpreP.data',
    'glycine': 'rama8000-gly-sym.data',
    'cis-proline': 'rama8000-cispro.data',
    'trans-proline': 'rama8000-transpro.data',
    'pre-proline': 'rama8000-prepro-noGP.data',
    'ile-val': 'rama8000-ileval-nopreP.data',
}
# Wie in MolProbity (cctbx mmtbx/validation/ramalyze.py)
FAVORED = 0.02
ALLOWED = {'general': 0.0005, 'cis-proline': 0.002}
ALLOWED_DEFAULT = 0.001
BINS, WIDTH = 180, 2.0


def read_grid(path):
    values = [0.0] * (BINS * BINS)  # Nicht aufgeführte Rasterpunkte sind 0
    for line in path.read_text(encoding='utf-8').splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        phi, psi, value = (float(v) for v in line.split())
        i, j = int((phi + 180) // WIDTH), int((psi + 180) // WIDTH)
        values[i * BINS + j] = value
    return values


def main(src, server_out, regions_out):
    grids = {name: read_grid(Path(src) / file) for name, file in FILES.items()}
    server = {
        '_source': 'Top8000 Ramachandran contour grids, Richardson Lab, Duke University, '
                   'https://github.com/rlabduke/reference_data',
        '_commit': COMMIT,
        '_license': 'CC BY 4.0, see LICENSE-top8000.txt',
        '_content': 'values × 10^6, 180 × 180 bins of 2°, phi-major, bin centres at -179 … 179',
        'grids': {name: [round(v * 1e6) for v in values] for name, values in grids.items()},
    }
    with gzip.open(server_out, 'wt', encoding='utf-8', compresslevel=9) as f:
        json.dump(server, f, separators=(',', ':'))
    regions = {
        'source': server['_source'], 'license': 'CC BY 4.0', 'bins': BINS,
        'regions': {
            name: ''.join('2' if v >= FAVORED else '1' if v >= ALLOWED.get(name, ALLOWED_DEFAULT) else '0'
                          for v in values)
            for name, values in grids.items()
        },
    }
    Path(regions_out).write_text(json.dumps(regions, separators=(',', ':')) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main(*sys.argv[1:])
