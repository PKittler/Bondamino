# Reference values for bond lengths, bond angles and the Ramachandran plot

This folder contains the reference data Bondamino compares structures with, the scripts that
produced it from the original sources, and the licences of those sources. See also
[LICENSE.md](../../../LICENSE.md) and [CITATION.md](../../../CITATION.md) in the repository root.

## Backbone: Engh & Huber (2001)

Engh, R. A. & Huber, R. (2001). Structure quality and target parameters.
*International Tables for Crystallography*, Vol. F, ch. 18.3, 382–392.
https://doi.org/10.1107/97809553602060000695

The values are stored as `BACKBONE_BOND_REFERENCES` and `BACKBONE_ANGLE_REFERENCES` in
`utils.py`. They were taken, together with the assignment rules (Gly, Pro, cis-Pro, residue
before or after the peptide bond), from the file `EnghHuber_IntlTblsF_1999.txt` of the program
Dangle (Richardson Lab, Duke University):

- https://github.com/rlabduke/javadev/blob/7a10fd2f4be7fa2f5162765d35bb8e975e9bce5b/dangle/resource/dangle/EnghHuber_IntlTblsF_1999.txt
- Licence: BSD-3-Clause, see `LICENSE-javadev.txt`

## Side chains: CCP4 Monomer Library

`ccp4_sidechains.json` contains the bond lengths (`value_dist`, `value_dist_esd`) and bond
angles (`value_angle`, `value_angle_esd`) of the 20 standard amino acids and common modified
amino acids, without hydrogen atoms. Only entries involving at least one atom outside the
backbone (N, CA, C, O, OXT) are included, because the library describes the amino acids as
free molecules.

- https://github.com/MonomerLibrary/monomers, commit `1c97ea23a7df06373ed01cd7c3156a1f8e9b27e3`
- Licence: LGPL-3.0, see `LICENSE-monomers.txt` and `GPL-3.0.txt`

To regenerate: download the `.cif` files of the residues listed in `build_ccp4_sidechains.py`
from the repository (folders `a/` to `z/`) into one folder, then run

    python -I build_ccp4_sidechains.py <folder> ccp4_sidechains.json

## Ramachandran regions: Top8000

Williams, C. J. et al. (2018). MolProbity: More and better reference data for improved
all-atom structure validation. *Protein Science* 27, 293–315. https://doi.org/10.1002/pro.3330

Grids of the Top8000 distributions (Richardson Lab, Duke University) for the six MolProbity
categories: general, glycine, cis-proline, trans-proline, pre-proline, and isoleucine/valine.
The evaluation follows MolProbity (`mmtbx/validation/ramalyze.py`): linear interpolation in the
grid, favoured from 2 %, allowed from 0.05 % (general), 0.2 % (cis-proline) or 0.1 % (all
other categories).

- https://github.com/rlabduke/reference_data/tree/5ee6875fc29eccc3c9dc7cbe705ff9fd1b505d7d/Top8000/Top8000_ramachandran_pct_contour_grids
- Licence: CC BY 4.0, see `LICENSE-top8000.txt`. Changes: the grids were converted to JSON,
  the values rounded to millionths, and the favoured/allowed/outlier regions derived from them
  for the plot.

Generated with `build_top8000_rama.py`:

- `top8000_rama.json.gz`: grid values for the evaluation on the server
- `static/pdb_grapher/data/rama8000_regions.json`: regions per grid cell for the plot

      python -I build_top8000_rama.py <folder with the rama8000-*.data files> top8000_rama.json.gz ../static/pdb_grapher/data/rama8000_regions.json
