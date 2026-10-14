# Licensing

This file explains under which terms Bondamino may be used, changed and passed on, which
third-party software and data Bondamino contains or needs, and what that means in practice.

> This is an explanation in plain language, not legal advice. If anything here seems to differ
> from the licence texts, the licence texts apply: [`LICENSE`](LICENSE) for Bondamino and the
> licence files listed in the tables below for the third-party components.

- [1. Bondamino's licence](#1-bondaminos-licence)
- [2. Additional terms](#2-additional-terms)
- [3. Third-party software and data](#3-third-party-software-and-data)
- [4. What this means for you](#4-what-this-means-for-you)
- [5. Contributions](#5-contributions)
- [6. Requests (not conditions)](#6-requests-not-conditions)

## 1. Bondamino's licence

Copyright © 2025–2026 Philipp Kittler

Bondamino is free software: you can redistribute it and/or modify it under the terms of the
**GNU Affero General Public License** as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version (SPDX identifier
`AGPL-3.0-or-later`), with the additional terms in [section 2](#2-additional-terms).

Bondamino is distributed in the hope that it will be useful, but **without any warranty**;
without even the implied warranty of merchantability or fitness for a particular purpose. See
the GNU Affero General Public License for more details. The full licence text is in
[`LICENSE`](LICENSE).

**In short, everyone may**

- use Bondamino for any purpose, including commercial use and use in research and teaching,
- study how it works and change it,
- pass it on, changed or unchanged, free of charge or for a fee,
- run it as a web service for others.

**In return**

- the copyright notices and licences must be kept and shown (see [section 2](#2-additional-terms)),
- changed versions must again be licensed under the AGPL-3.0 (version 3 or later) as a whole,
- whoever passes Bondamino on, or lets other people use a *changed* version over a network,
  must give them the complete source code of that version under the same licence,
- no one may add further restrictions to the licence.

The AGPL-3.0 is the GPL-3.0 plus one rule for web applications (section 13): if you run a
changed version on a server, the people who use it through their browser must be able to get
the source code of exactly that version. This keeps improvements to Bondamino open, even when
they are only offered as an online service.

Results you produce with Bondamino (graphs, tables, PDF reports, exported data) are yours. The
licence covers the software, not your analyses. Cite the data sources as described in
[CITATION.md](CITATION.md).

## 2. Additional terms

As permitted by section 7 of the AGPL-3.0, the following additional terms apply to Bondamino.
They are also binding for everyone who passes on Bondamino or a work based on it:

1. **Licences must be shown** (section 7 b). Every copy and every changed version must keep
   - the copyright notices,
   - the files `LICENSE` and `LICENSE.md` (an updated `LICENSE.md` is fine as long as no
     component is left out),
   - the licence files of the third-party software and data sources listed in
     [section 3](#3-third-party-software-and-data), in a place where users can find them.
2. **Changed versions must be marked** (section 7 c). A changed version must say in a
   reasonable way that it has been changed (for example in its README and in the commit
   history) and must not be presented as the original Bondamino by Philipp Kittler.

Section 7 of the AGPL-3.0 allows no other additional conditions. Everything else follows from
the AGPL-3.0 itself and from the licences of the components below. Where a component has its
own licence, that licence applies to the component, and nothing in Bondamino's terms
restricts the rights it grants.

## 3. Third-party software and data

### 3.1 Included in this repository

These files are part of Bondamino's source code and of every copy of it. Their licences are
all compatible with the AGPL-3.0.

| Component | Version | Used for | Licence | Licence text |
|---|---|---|---|---|
| [Bootstrap](https://getbootstrap.com/) | 5.3.8 | Page layout, buttons, dialogs | MIT | `bondamino/pdb_grapher/static/pdb_grapher/vendor/bootstrap/LICENSE` |
| [D3.js](https://d3js.org/) | 7.9.0 | Interactive graphs and plots | ISC | `.../vendor/d3/LICENSE` |
| [jsPDF](https://github.com/parallax/jsPDF) | 4.2.1 | PDF export (contains further small libraries; their notices are kept in the file) | MIT | `.../vendor/jspdf/LICENSE-jspdf.txt` |
| [jsPDF-AutoTable](https://github.com/simonbengtsson/jsPDF-AutoTable) | 5.0.8 | Tables in the PDF export | MIT | `.../vendor/jspdf/LICENSE-jspdf-autotable.txt` |
| [Roboto, Roboto Mono](https://fonts.google.com/specimen/Roboto) (via [Fontsource](https://fontsource.org/)) | 5.3.0 | Fonts of the user interface | SIL Open Font License 1.1 | `.../vendor/fonts/LICENSE-roboto.txt`, `LICENSE-roboto-mono.txt` |

**Data sources** (folder `bondamino/pdb_grapher/reference_data/`, details in its
[README](bondamino/pdb_grapher/reference_data/README.md)):

| Data | Used for | Origin | Licence | Licence text |
|---|---|---|---|---|
| Backbone reference values after Engh & Huber (2001), in `utils.py` | Bond lengths and angles of the backbone | File `EnghHuber_IntlTblsF_1999.txt` of Dangle (javadev), Richardson Lab, Duke University | BSD-3-Clause | `reference_data/LICENSE-javadev.txt` |
| `ccp4_sidechains.json` (extract) | Bond lengths and angles of the side chains | [CCP4 Monomer Library](https://github.com/MonomerLibrary/monomers) | LGPL-3.0 | `reference_data/LICENSE-monomers.txt`, `GPL-3.0.txt` |
| `top8000_rama.json.gz`, `static/pdb_grapher/data/rama8000_regions.json` (converted) | Ramachandran evaluation and plot regions | Top8000 contour grids, [Richardson Lab](https://github.com/rlabduke/reference_data) | CC BY 4.0 | `reference_data/LICENSE-top8000.txt` |

The scientific values themselves (numbers published in the literature) should be cited as
described in [CITATION.md](CITATION.md).

### 3.2 Installed with pip (not part of this repository)

These packages are downloaded when Bondamino is installed (`requirements*.txt`). They are not
contained in the source code, but they are contained in the Docker image.

| Package | Used for | Licence | Installed with |
|---|---|---|---|
| [Django](https://www.djangoproject.com/) | Web framework | BSD-3-Clause | all |
| asgiref, sqlparse | Parts of Django | BSD-3-Clause | all |
| tzdata | Time zone database (for the time zone setting) | Apache-2.0 | all |
| [WhiteNoise](https://whitenoise.readthedocs.io/) | Delivers CSS, JavaScript and fonts | MIT | all |
| [gunicorn](https://gunicorn.org/) | Web server (Linux/macOS, Docker) | MIT | all, not on Windows |
| [waitress](https://github.com/Pylons/waitress) | Web server on Windows (optional) | ZPL 2.1 | manually, see [INSTALL.md](INSTALL.md) |
| [psycopg 3](https://www.psycopg.org/) (incl. `psycopg-binary`) | PostgreSQL driver; the binary package contains the PostgreSQL client library and further libraries under their own licences | LGPL-3.0 | `requirements-postgresql.txt` |
| [mysqlclient](https://github.com/PyMySQL/mysqlclient) | MariaDB/MySQL driver | GPL-2.0-or-later | `requirements-mariadb.txt` |

### 3.3 Additionally contained in the Docker image

The image built from the `Dockerfile` contains, besides Bondamino and the packages from 3.2,
Python (PSF License 2.0), a minimal Debian system (various free licences, listed in the image
under `/usr/share/doc/*/copyright`) and the MariaDB client library `libmariadb3`
(LGPL-2.1-or-later). The database servers started by the Compose files (`postgres`, `mariadb`)
are separate, unchanged official images with their own licences (PostgreSQL License,
GPL-2.0); they are not part of Bondamino.

### 3.4 Data loaded at run time

Structures imported by PDB ID are downloaded from the [RCSB Protein Data
Bank](https://www.rcsb.org/). PDB archive data is available under CC0 1.0 (public domain);
the PDB asks to cite the structure's primary publication and the PDB (see
[CITATION.md](CITATION.md)).

## 4. What this means for you

**You use Bondamino in your lab, unchanged.** Nothing to do. Your results are yours. Please
cite as described in [CITATION.md](CITATION.md).

**You change Bondamino for your own group and run it on your institute's server.** The people
who use it through the browser must be able to get the source code of your changed version
(AGPL section 13), e.g. through a public repository or a link in the application. Mark the
version as changed. Keep all licence files.

**You pass on Bondamino or a changed version** (a fork on GitLab or elsewhere, a ZIP file, a Docker
image). Include `LICENSE`, `LICENSE.md` and all licence files from section 3, license your
version under the AGPL-3.0 (version 3 or later), and provide the complete source code.
If you distribute a Docker image, the components from 3.2 and 3.3 are passed on with it;
their source code must be obtainable as well. For the GPL- and LGPL-licensed packages, a
common way is to name the exact versions and point to their unchanged source releases (e.g.
on PyPI and in the Debian archive); you remain responsible that these sources stay available.

**You reuse parts of Bondamino's code in another project.** The other project must be
compatible with the AGPL-3.0: as a whole it must be licensed under the AGPL-3.0. (A GPL-3.0
project may also combine its code with Bondamino's, as section 13 of both licences allows;
Bondamino's parts then stay under the AGPL-3.0.) You may not use Bondamino's code in
closed-source software.

**You only use the reference data** (outside Bondamino). Only the licence of that data set
applies: BSD-3-Clause for the Engh & Huber file from Dangle, LGPL-3.0 for the CCP4 extract,
CC BY 4.0 (attribution, mark changes) for the Top8000 data.

**You choose a licence for your own version.** It must be the AGPL-3.0 (version 3 or any later
version published by the Free Software Foundation) for the work as a whole. You may license
your own additions additionally under a more permissive, AGPL-compatible licence (for example
MIT), but the combined work stays under the AGPL. A licence that forbids changes, commercial
use or passing on, or that removes the requirement to show the licences, is not allowed.

## 5. Contributions

Contributions to this repository (pull requests, patches) are accepted under the same licence
as Bondamino: AGPL-3.0-or-later with the additional terms in section 2. By submitting a
contribution you confirm that you have the right to do so.

## 6. Requests (not conditions)

The following are requests. They are not licence conditions, because the AGPL-3.0 does not
allow additional conditions of this kind:

- If you use Bondamino for work that you publish, please cite it as described in
  [CITATION.md](CITATION.md).
- If you improve Bondamino, please consider contributing the changes back to the original
  repository.

How Bondamino was developed, including the use of AI, is described in
[TRANSPARENCY.md](TRANSPARENCY.md).
