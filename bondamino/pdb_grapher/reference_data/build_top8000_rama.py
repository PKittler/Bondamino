"""Builds compact Ramachandran grid files (Top8000, Richardson Lab) for server and browser.

One-off maintenance script, not used at runtime.

Usage: python -I build_top8000_rama.py <folder with the rama8000-*.data files> top8000_rama.json.gz <regions.json>

The input files are the Top8000 contour grids at commit `COMMIT` (see README.md), one per
MolProbity category (`FILES`). Outputs:

- top8000_rama.json.gz: per category the 180 × 180 grid values (φ outer, ψ inner index,
  bin centres at −179°, −177°, … 179°) as integers in millionths; read by
  `utils.ramachandran_value()` to score residues on the server.
- <regions.json> (normally static/pdb_grapher/data/rama8000_regions.json): per category a
  string with one character per grid bin, 0 (outlier), 1 (allowed) or 2 (favored); used by
  the browser to draw the regions in the Ramachandran plot.
"""
import gzip
import json
import sys
from pathlib import Path

# Commit of https://github.com/rlabduke/reference_data the .data files were taken from.
COMMIT = '5ee6875fc29eccc3c9dc7cbe705ff9fd1b505d7d'
# Ramachandran category (as used in utils.evaluate_ramachandran) -> Top8000 grid file.
FILES = {
    'general': 'rama8000-general-noGPIVpreP.data',
    'glycine': 'rama8000-gly-sym.data',
    'cis-proline': 'rama8000-cispro.data',
    'trans-proline': 'rama8000-transpro.data',
    'pre-proline': 'rama8000-prepro-noGP.data',
    'ile-val': 'rama8000-ileval-nopreP.data',
}
# Region thresholds on the grid value, as in MolProbity (cctbx mmtbx/validation/ramalyze.py);
# must match RAMA_FAVORED and RAMA_ALLOWED in utils.py.
FAVORED = 0.02
ALLOWED = {'general': 0.0005, 'cis-proline': 0.002}
ALLOWED_DEFAULT = 0.001
# Number of bins per axis and bin width in degrees (180 × 2° = 360°).
BINS, WIDTH = 180, 2.0


def read_grid(path):
    """Reads a Top8000 grid file.

    Args:
        path: `pathlib.Path` of a rama8000-*.data file with lines "phi psi value";
            lines starting with '#' are comments.

    Returns:
        A flat list of BINS × BINS values, index `i * BINS + j` for φ bin i and ψ bin j.
    """
    values = [0.0] * (BINS * BINS)  # Grid points not listed in the file are 0.
    for line in path.read_text(encoding='utf-8').splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        phi, psi, value = (float(v) for v in line.split())
        i, j = int((phi + 180) // WIDTH), int((psi + 180) // WIDTH)
        values[i * BINS + j] = value
    return values


def main(src, server_out, regions_out):
    """Converts the grid files and writes both output files.

    Args:
        src: Folder containing the files listed in `FILES`.
        server_out: Path of the gzipped JSON with the grid values (top8000_rama.json.gz).
        regions_out: Path of the JSON with the region strings for the browser.
    """
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
