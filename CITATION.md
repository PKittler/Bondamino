# How to cite Bondamino

If Bondamino contributed to work that you publish (article, thesis, poster, report), please
cite Bondamino itself **and** the data sources behind the results you report. This is good
scientific practice; for the data sources it is also how their authors ask to be credited.

- [Citing Bondamino](#citing-bondamino)
- [Citing the data sources](#citing-the-data-sources)
- [Software packages](#software-packages)
- [Example text for a methods section](#example-text-for-a-methods-section)
- [References](#references)

## Citing Bondamino

> Kittler, P. (2026). *Bondamino: interactive analysis of bond geometry in protein
> structures* [Computer software]. https://gitlab.com/d11m-group/bondamino

Please also state **which version** you used. Bondamino has no numbered releases yet, so give
the date and, if you know it, the commit (shown on GitLab next to the latest change, or
with `git rev-parse --short HEAD`), e.g. "version of 14 October 2026, commit `0096f02`".

BibTeX:

```bibtex
@software{bondamino,
  author  = {Kittler, Philipp},
  title   = {Bondamino: interactive analysis of bond geometry in protein structures},
  year    = {2026},
  url     = {https://gitlab.com/d11m-group/bondamino},
  license = {AGPL-3.0-or-later},
  note    = {Version of <date>, commit <hash>}
}
```

The repository also contains a [`CITATION.cff`](CITATION.cff) file with the same information in
a machine-readable form (Citation File Format). Reference managers such as Zotero and archives
such as Zenodo can read it; `cffconvert` turns it into BibTeX, APA and other formats.

## Citing the data sources

Bondamino compares structures with published reference data. Cite the sources of the values
you report:

| You report … | Cite |
|---|---|
| Bond lengths or bond angles of the **backbone** compared with reference values (Z-scores, RMSZ, outliers) | Engh & Huber (2001) |
| Bond lengths or bond angles of **side chains** or modified amino acids compared with reference values | Vagin et al. (2004), for the CCP4 Monomer Library |
| The **Ramachandran evaluation** (favoured / allowed / outlier) or the plot with its regions | Williams et al. (2018) for the Top8000 data; optionally Chen et al. (2010) for the MolProbity method and Ramachandran et al. (1963) for the plot itself |
| Results for **structures from the Protein Data Bank** | The structure's primary publication, the entry's DOI (`https://doi.org/10.2210/pdbXXXX/pdb`, e.g. `10.2210/pdb1lyz/pdb` for 1LYZ), and Berman et al. (2000) for the PDB |

The Top8000 data are licensed under CC BY 4.0, which requires attribution. Citing Williams et
al. (2018) satisfies this in a publication. If you redistribute the data files themselves,
keep their licence file and say what was changed (see [LICENSE.md](LICENSE.md)).

## Software packages

Bondamino is built on Django, D3.js, Bootstrap, jsPDF and other packages listed in
[LICENSE.md](LICENSE.md). You do **not** need to cite them in a publication; their licence
conditions are met by Bondamino itself. If you want to acknowledge them, an acknowledgement
such as "Bondamino uses Django, D3.js and jsPDF" is sufficient.

How Bondamino was developed, including the use of AI, is described in
[TRANSPARENCY.md](TRANSPARENCY.md). Bondamino computes all results with its own code; no AI is
involved when you use it.

## Example text for a methods section

> Bond lengths, bond angles and backbone dihedral angles were analysed with Bondamino
> (Kittler, 2026; version of 14 October 2026). Backbone geometry was compared with the
> reference values of Engh & Huber (2001), side-chain geometry with the CCP4 Monomer Library
> (Vagin et al., 2004); deviations were expressed as Z-scores and root-mean-square Z-scores
> (RMSZ), with |Z| > 4 counted as outliers. Backbone conformations were evaluated against the
> Top8000 Ramachandran distributions (Williams et al., 2018). The structure 1LYZ was obtained
> from the Protein Data Bank (Berman et al., 2000).

Adapt the text to what you actually used. The definitions of Z-score, RMSZ and the thresholds
are given in [DOCUMENTATION.md](DOCUMENTATION.md).

## References

Berman, H. M., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T. N., Weissig, H.,
Shindyalov, I. N. & Bourne, P. E. (2000). The Protein Data Bank. *Nucleic Acids Research*,
28, 235–242. https://doi.org/10.1093/nar/28.1.235

Chen, V. B., Arendall, W. B., Headd, J. J., Keedy, D. A., Immormino, R. M., Kapral, G. J.,
Murray, L. W., Richardson, J. S. & Richardson, D. C. (2010). MolProbity: all-atom structure
validation for macromolecular crystallography. *Acta Crystallographica Section D*, 66, 12–21.
https://doi.org/10.1107/S0907444909042073

Engh, R. A. & Huber, R. (2001). Structure quality and target parameters. In M. G. Rossmann &
E. Arnold (Eds.), *International Tables for Crystallography, Vol. F: Crystallography of
Biological Macromolecules* (ch. 18.3, pp. 382–392). Kluwer Academic Publishers.
https://doi.org/10.1107/97809553602060000695

Kittler, P. (2026). *Bondamino: interactive analysis of bond geometry in protein structures*
[Computer software]. https://gitlab.com/d11m-group/bondamino

Ramachandran, G. N., Ramakrishnan, C. & Sasisekharan, V. (1963). Stereochemistry of
polypeptide chain configurations. *Journal of Molecular Biology*, 7, 95–99.
https://doi.org/10.1016/S0022-2836(63)80023-6

Vagin, A. A., Steiner, R. A., Lebedev, A. A., Potterton, L., McNicholas, S., Long, F. &
Murshudov, G. N. (2004). REFMAC5 dictionary: organization of prior chemical knowledge and
guidelines for its use. *Acta Crystallographica Section D*, 60, 2184–2195.
https://doi.org/10.1107/S0907444904023510

Williams, C. J., Headd, J. J., Moriarty, N. W., Prisant, M. G., Videau, L. L., Deis, L. N.,
Verma, V., Keedy, D. A., Hintze, B. J., Chen, V. B., Jain, S., Lewis, S. M., Arendall, W. B.,
Snoeyink, J., Adams, P. D., Lovell, S. C., Richardson, J. S. & Richardson, D. C. (2018).
MolProbity: More and better reference data for improved all-atom structure validation.
*Protein Science*, 27, 293–315. https://doi.org/10.1002/pro.3330
