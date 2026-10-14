# Bondamino

**Bond geometry of protein structures at a glance.**

Bondamino is a web application for structural biologists. It reads a protein structure in PDB
format and shows its bond lengths, bond angles and backbone dihedral angles as an interactive
graph, compares them with published reference values and points out unusual geometry. A work
group can keep its structures in one place, look at them together and discuss them in notes.

Bondamino runs in the browser. It is installed once on a computer or server, for example at
your institute, and everyone in the group signs in with their own account.

[![License: AGPL-3.0-or-later](https://img.shields.io/badge/license-AGPL--3.0--or--later-blue)](LICENSE.md)

## Features

- **Load structures** by uploading a `.pdb` file or by entering a PDB ID (e.g. `1LYZ`), which is
  fetched from the RCSB Protein Data Bank.
- **Interactive graph** of every residue with its atoms and bonds, coloured by bond length,
  by deviation from the reference value or by B-factor; zoom, tooltips and a sequence bar for
  navigation.
- **Comparison with reference values**: backbone after Engh & Huber (2001), side chains after
  the CCP4 Monomer Library; Z-scores, RMSZ and a list of outliers.
- **Ramachandran plot** with the Top8000 regions (MolProbity method) and an evaluation of
  every residue as favoured, allowed or outlier.
- **Real-world structures**: multiple models (e.g. NMR ensembles, including an averaged
  model), alternate conformations, chain breaks, disulfide bridges, ligands and modified amino
  acids.
- **Exports**: PDF report, CSV tables for spreadsheets, JSON data, and the Ramachandran plot
  and histograms as images (SVG, PNG, JPEG) for figures.
- **Work groups and notes**: share structures with your group, write notes visible to you or
  to selected groups.
- **For the whole group**: user accounts with invitation emails, English and German user
  interface, light and dark mode.

## Getting started

How you install Bondamino depends on where it should run. The quickest way to try it on your
own computer is Docker:

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/) and start it.
2. Download Bondamino (on this page **Code → Download source code → zip**, then unpack it and
   rename the folder to `Bondamino`) or clone it with Git.
3. In the Bondamino folder, copy the file `.env.example` to `.env`, fill in
   `ADMIN_USERNAME`, `ADMIN_EMAIL` and `ADMIN_PASSWORD` (your first administrator account) and
   empty the line `EMAIL_BACKEND=` (emails are then written to the log instead of being sent).
4. Open a terminal in the Bondamino folder and run

   ```sh
   docker compose -f docker-compose.sqlite.yml up -d --build
   ```

5. Open http://localhost:8000 in your browser and sign in.

This setup is meant for trying Bondamino out. For a server that your group uses every day, with
PostgreSQL or MariaDB, HTTPS, email and backups, and for installation without Docker (with
gunicorn) or in Portainer, follow [INSTALL.md](INSTALL.md).

## Documentation

| Document | Contents |
|---|---|
| [INSTALL.md](INSTALL.md) | Installation step by step: on your own computer, on a server with gunicorn, with Docker, in Portainer; configuration, updates, backups |
| [TUTORIAL.md](TUTORIAL.md) | A complete analysis from the first sign-in to the export, and guides for individual features |
| [DOCUMENTATION.md](DOCUMENTATION.md) | Reference of every page, option and export, and the scientific definitions (Z-score, RMSZ, thresholds) |
| [CITATION.md](CITATION.md) | How to cite Bondamino and the data sources in publications |
| [LICENSE.md](LICENSE.md) | Licence of Bondamino and of the software and data it uses, and what that means for you |
| [TRANSPARENCY.md](TRANSPARENCY.md) | How Bondamino was developed with AI (Claude) |

## Citing

If Bondamino contributed to your publication, please cite it together with the reference data
you used; see [CITATION.md](CITATION.md). Machine-readable citation data is in
[CITATION.cff](CITATION.cff).

## Licence

Bondamino is free software under the GNU Affero General Public License v3.0 or later, with
additional terms that require the licences of Bondamino and of its components to be kept. You
may use, change and pass it on; changed versions must stay under the same licence. Details,
including all third-party software and data sources, are in [LICENSE.md](LICENSE.md).

## Development with AI

Bondamino was developed with the help of Claude (Anthropic). What this means and how the
results were checked is described in [TRANSPARENCY.md](TRANSPARENCY.md).
