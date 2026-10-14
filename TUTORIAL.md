# Bondamino tutorial

This tutorial shows you how to use Bondamino for your own research. It is written for
scientists, not for software developers. You do not need to install anything: you only need a
web browser and an account on a Bondamino server.

The tutorial has two parts:

1. A **standard workflow** that takes you from your first sign-in to a finished analysis that
   you can use in a paper, a thesis or a group meeting.
2. **Task guides** ("How do I …?") for individual features, with hints on when they are useful.

Other documents:

- [README.md](README.md): what Bondamino is.
- [INSTALL.md](INSTALL.md): installing and running a Bondamino server (for administrators).
- [DOCUMENTATION.md](DOCUMENTATION.md): complete feature reference and the scientific
  background, with exact definitions and thresholds.
- [CITATION.md](CITATION.md): how to cite Bondamino and the reference data.

Button and menu names are written in **bold**, exactly as the English (UK) interface shows them.
If you use Bondamino in another language, the positions of the buttons are the same.

## Contents

- [Before you start](#before-you-start)
  - [What you need](#what-you-need)
  - [Terms used in this tutorial](#terms-used-in-this-tutorial)
  - [Example structures](#example-structures)
- [Standard workflow](#standard-workflow)
  - [Step 1: Set your password and sign in](#step-1-set-your-password-and-sign-in)
  - [Step 2: Load a structure](#step-2-load-a-structure)
  - [Step 3: Read the structure information](#step-3-read-the-structure-information)
  - [Step 4: Explore the graph](#step-4-explore-the-graph)
  - [Step 5: Check the bond statistics](#step-5-check-the-bond-statistics)
  - [Step 6: Look at the outliers](#step-6-look-at-the-outliers)
  - [Step 7: Check the angles and the Ramachandran plot](#step-7-check-the-angles-and-the-ramachandran-plot)
  - [Step 8: Look at the atom data](#step-8-look-at-the-atom-data)
  - [Step 9: Export your results](#step-9-export-your-results)
  - [Step 10: Share and discuss with your work group](#step-10-share-and-discuss-with-your-work-group)
  - [Step 11: Cite Bondamino](#step-11-cite-bondamino)
  - [Workflow checklist](#workflow-checklist)
- [Task guides](#task-guides)
  - [How do I compare the models of an NMR ensemble?](#how-do-i-compare-the-models-of-an-nmr-ensemble)
  - [How do I look at alternate conformations?](#how-do-i-look-at-alternate-conformations)
  - [How do I find chain breaks?](#how-do-i-find-chain-breaks)
  - [How do I check disulfide bridges?](#how-do-i-check-disulfide-bridges)
  - [How do I find cis peptide bonds?](#how-do-i-find-cis-peptide-bonds)
  - [How do I look at ligands and modified amino acids?](#how-do-i-look-at-ligands-and-modified-amino-acids)
  - [How do I find flexible regions with the B-factor?](#how-do-i-find-flexible-regions-with-the-b-factor)
  - [How do I work with a large structure?](#how-do-i-work-with-a-large-structure)
  - [How do I prepare a figure?](#how-do-i-prepare-a-figure)
  - [How do I prepare data for a spreadsheet?](#how-do-i-prepare-data-for-a-spreadsheet)
  - [How do I change who can see a file, or remove it?](#how-do-i-change-who-can-see-a-file-or-remove-it)
  - [What can a group leader do?](#what-can-a-group-leader-do)
  - [How do I write and manage notes?](#how-do-i-write-and-manage-notes)
  - [How do I change my personal settings?](#how-do-i-change-my-personal-settings)
  - [What if a structure cannot be loaded?](#what-if-a-structure-cannot-be-loaded)
- [Troubleshooting](#troubleshooting)

---

## Before you start

### What you need

- **An account.** You cannot register yourself. Your Bondamino administrator creates your
  account and adds you to your work group(s). You then receive an email with the subject
  "Your access to Bondamino". It contains your username and a link to set your password.
- **A web browser**, for example Firefox, Chrome, Edge or Safari, in a current version.
- **The address of your Bondamino server.** Ask your administrator. If you want to run a server
  yourself, see [INSTALL.md](INSTALL.md).
- **A structure.** Either a PDB ID (four characters, for example `1LYZ`) or a structure file in
  PDB format (file name ending in `.pdb`).

### Terms used in this tutorial

| Term | Meaning |
|---|---|
| PDB file | A text file with the atom coordinates of a structure, in the classic format of the Protein Data Bank. |
| PDB ID | The four-character code of a structure in the Protein Data Bank, for example `1LYZ`. |
| Model | One complete set of coordinates. NMR structures usually contain several models (an "ensemble"); crystal structures usually contain one. |
| Chain | One polypeptide (or a group of other molecules) with its own letter, for example chain `A`. |
| Residue | One building block of a chain: an amino acid, or a ligand, ion or water molecule. |
| Alternate conformation (altLoc) | Some atoms are modelled in two or more positions, labelled `A`, `B`, … in the file. |
| Occupancy | The fraction of molecules in the crystal that have an atom in this position (1.00 = always). |
| B-factor | A number per atom that describes how much its position is smeared out (motion or disorder). Higher means less well defined. |
| Reference value | The expected bond length or bond angle, with a standard deviation σ, taken from carefully determined small-molecule data. |
| Z-score | How far a value is from its reference, in standard deviations: Z = (value − reference) / σ. |
| RMSZ | Root mean square of all Z-scores. A single number for the overall agreement with the reference values. |
| Outlier | A bond or angle with \|Z\| > 4, or an amino acid outside the allowed Ramachandran region. |
| φ, ψ, ω | Backbone torsion (dihedral) angles. φ and ψ define the Ramachandran plot; ω describes the peptide bond (trans ≈ 180°, cis ≈ 0°). |
| Work group | A group of users (for example a lab) who share files and notes. Every group has a group leader. |

The exact definitions and thresholds are in [DOCUMENTATION.md](DOCUMENTATION.md).

### Example structures

You can follow every step of this tutorial with public structures. Bondamino loads them
directly from the Protein Data Bank (RCSB).

| PDB ID | What it is | Good for |
|---|---|---|
| `1LYZ` | Hen egg-white lysozyme, X-ray, 2.0 Å, one of the first protein structures (1970s) | The standard workflow |
| `4LZT` | Hen egg-white lysozyme, X-ray, 0.95 Å | Comparison with a modern lysozyme structure |
| `1D3Z` | Ubiquitin, solution NMR, 10 models | Models and the averaged model |
| `1EJG` | Crambin, X-ray, 0.54 Å, alternate conformations `A`, `B` and `C` | Alternate conformations |
| `1ATP` | Protein kinase A with ATP, two Mn²⁺ ions and an inhibitor peptide | Ligands and modified amino acids |
| `1HZH` | Complete IgG1 antibody, X-ray, 2.7 Å, about 10,600 atoms | Large structures, chain breaks, disulfide bridges |

The numbers quoted below for these structures were checked with the files from RCSB at the
time of writing. If RCSB updates an entry, small differences are possible.

---

## Standard workflow

This workflow uses `1LYZ`. It takes about 30 minutes the first time.

### Step 1: Set your password and sign in

1. Open the email "Your access to Bondamino" and click the link.
   You see the page **Set password** with your username.
2. Enter a new password twice and click **Save password**.
   The password must not be too short, too common, entirely numeric or too similar to your
   username. If it is not accepted, a red message under the field tells you why.
3. You are now signed in and see the start page with the list of files.
   A green message says "Your password has been saved."

The link is valid for 24 hours and works only once. If it has expired, see
[Troubleshooting](#troubleshooting).

Next time, open the address of your server, enter your **Username** and **Password** and click
**Sign in**.

> Tip: The language and the light/dark display can already be changed on the sign-in page (top
> right). After signing in, both are in the account menu, see
> [How do I change my personal settings?](#how-do-i-change-my-personal-settings).

### Step 2: Load a structure

The start page shows all files you may see. The columns are **File**, **Visible to**,
**Uploaded** and **Latest note**.

1. Click the blue button with the document-and-plus icon at the top left (tooltip: "Upload a
   new PDB file or load one from RCSB").
   The dialog **New PDB file** opens.
2. Choose **one** source:
   - **From the Protein Data Bank:** type `1LYZ` into the field **Load PDB ID from RCSB**.
   - **From your computer:** drag a `.pdb` file into the dashed area, or click **Choose file**.
     The file name appears in the area. Files may be up to 20 MB.

   If you fill in both, Bondamino uses the one you chose last and clears the other.
3. Under **Available to**, choose who can see the file:
   - **For me personally**: you can always see your own uploads.
   - **Working group …**: one checkbox per work group you belong to. All members of the
     selected groups can see the file.

   At least one box must be ticked. For your own first tests, **For me personally** alone is
   fine. You can change this later (see
   [How do I change who can see a file, or remove it?](#how-do-i-change-who-can-see-a-file-or-remove-it)).
4. Click **Analyse protein**.
   A spinning wheel appears while Bondamino downloads and checks the file. This can take a few
   seconds.
5. The analysis page of the structure opens.

What you see: five tabs at the top (**Graph**, **Ramachandran**, **Data**, **Statistics**,
**Notes**). On the right of the tabs are the file name (`1LYZ.pdb`), a download button for the
original file (**Download PDB file**) and two selection boxes for the model and the
conformation. For `1LYZ` both boxes are grey and show **Model 1** and **Single conformation**,
because the file has only one model and no alternate conformations.

If something goes wrong, a red message in the dialog tells you why. See
[What if a structure cannot be loaded?](#what-if-a-structure-cannot-be-loaded).

### Step 3: Read the structure information

Before you interpret any numbers, check what kind of structure you have.

1. Click the tab **Statistics**.
2. Look at the section **General**. Next to the table **Elements** you find:
   - **Models**: the number of models in the file.
   - **Method**: for example *X-ray crystallography* or *Solution NMR*.
   - **Resolution**: for crystal and cryo-EM structures, in Å.
   - **R work** and **R free** (shown with subscripts): the agreement of a crystal structure
     with its diffraction data. Lower is better. R free is the more honest of the two, because
     it uses data that were left out of the refinement.

   Fields that are missing from the file header are not shown.
3. Look at the tables **Elements** and **Amino acids**. They show what the file contains.

For `1LYZ` you see **Models** 1, **Method** X-ray crystallography and **Resolution** 2.00 Å.
There are no R values, because this early entry does not report them. The protein has 129
amino acids.

Why this matters: resolution and method tell you how much detail you can expect. A 3 Å crystal
structure or an NMR ensemble is expected to deviate more from ideal geometry than a 1 Å
crystal structure.

### Step 4: Explore the graph

Click the tab **Graph**. Each amino acid is one horizontal band, from the N terminus at the top
to the C terminus at the bottom:

- On the left: the one-letter code, the three-letter code and the position, for example
  **K / LYS / #1**.
- A vertical column with the backbone atoms N, CA and C; the carbonyl O sits to the left.
- The side chain grows to the right as a tree. Rings (Phe, Tyr, Trp, His, Pro) are drawn as
  polygons.
- Every bond is labelled with its length in Å.
- Solid lines are backbone bonds, dashed lines side-chain bonds, dash-dotted lines ligand
  bonds. Disulfide bridges are drawn as brackets on the right. The **Legend** on the right
  explains the line types and the element colours.

Atoms without any bond (hydrogen atoms, water, single metal ions) are not drawn in the graph.
You find them in the **Data** tab and in the statistics.

**Moving around**

| Action | Effect |
|---|---|
| Mouse wheel | Scroll up and down |
| Shift + mouse wheel | Scroll left and right |
| Ctrl + mouse wheel, or pinch on a touchpad | Zoom in and out |
| Drag with the mouse | Move the view |
| Round arrow button in the toolbar (**Reset zoom**) | Back to the start view |

**Tooltips:** hold the mouse pointer over a bond. A small box shows the two atoms, the length,
the bond type and the reference value with the Z-score. Over an atom you see its name, residue,
chain, B-factor and all bond angles with this atom at the centre.

**Sequence bar:** above the graph, the sequence is shown in one-letter code, with numbers
above it (chain start and every tenth position). Click a letter and the graph jumps to that
amino acid. While you scroll, the residues currently in view are highlighted. Ligands, ions and
water appear as `?`.

**Bond colours** (first group of buttons in the toolbar):

| Button | What the colour means | When to use it |
|---|---|---|
| **By length** | Red = shortest, green = longest bond in this structure | To see the range of lengths at a glance. The scale is relative: red does not mean "wrong". |
| **Uniform** | All bonds the same colour | For a calm picture, for example for teaching |
| **Reference** | Green \|Z\| ≤ 2, yellow 2–4, red \|Z\| > 4 (outlier), grey no reference | To find problem spots quickly. Usually the most useful view. |
| **B-factor** | Atoms and bonds from blue (low B) to red (high B) | To find flexible or poorly defined regions |

**Layout** (second group of buttons):

- **Schematic**: a fixed grid. Every amino acid looks the same, so the graph is easy to read.
- **Proportional**: the drawn length of each bond is proportional to its real length in Å.
  Unusually long or short bonds stand out. This holds for the backbone, the side chains and the
  peptide bonds; ligand bonds and the disulfide brackets are not drawn to scale. Note: this is
  still a 2D drawing, not the 3D structure.

Bondamino remembers your choice of colour and layout in this browser.

Try it: click **Reference**. In `1LYZ` many bonds turn yellow or red. Step 5 explains why.

### Step 5: Check the bond statistics

1. Click **Statistics** and scroll to the section **Bonds**.
2. Read the first block, **Deviation from reference values**:
   - **Assessed bonds**: how many bonds have a reference value (for `1LYZ`: 1,020 of 1,025).
     Some bonds, for example those in ligands and the disulfide bridges, have no reference
     value and are not rated.
   - **RMSZ**: the overall agreement with the references.
   - **Outliers (|Z| > 4)**: the number of clearly deviating bonds.
3. The table below lists the RMSZ and outliers per section: first **Backbone**, then the side
   chains of each amino acid. Click a row to expand it and see the individual atom pairs, with
   reference value, mean length, mean deviation and RMSZ.
4. The source of the reference values is printed under the table: Engh & Huber (2001) for the
   backbone, CCP4 Monomer Library for the side chains.

**How to read the RMSZ.** As a rough guide, values around 1 or below are typical for
structures refined with modern geometry restraints. Clearly higher values mean that the
geometry deviates more than expected. The value depends on resolution and refinement method,
so compare structures of similar quality. Details are in
[DOCUMENTATION.md](DOCUMENTATION.md).

**Example.** `1LYZ` has an RMSZ of about 12.8 and several hundred outliers. This does not mean
that Bondamino is wrong. `1LYZ` is one of the very first protein structures; its coordinates
were not refined with today's geometry restraints. Load `4LZT` (lysozyme at 0.95 Å) for
comparison: its RMSZ is about 1.3 with only a few outliers. This comparison is a good first
exercise to get a feeling for the numbers.

Further down in the **Bonds** section:

- **Peptide bonds**: how many of the theoretically possible peptide bonds exist. Missing ones
  are chain breaks (see [How do I find chain breaks?](#how-do-i-find-chain-breaks)).
- **Disulfide bridges**: see [How do I check disulfide bridges?](#how-do-i-check-disulfide-bridges).
- **Shortest and longest bond**: in all models and in the selected model, compared with the
  mean and median of the same bond over all models.
- **Distribution of bond lengths**: a histogram with 0.01 Å bins.
- **Bonds by atoms involved**: all bonds grouped by bond type and atom pair, for example
  "C–N (peptide bond)" or "CA–CB (side chain)", with count, mean, median, minimum and maximum.
  A C–O in a ligand is therefore counted separately from the backbone C–O. Click a row to see
  every single bond.

### Step 6: Look at the outliers

1. In **Statistics → Bonds**, find the table **Outliers (|Z| > 4)**. It lists every outlier with
   **Atoms**, **Chain**, **Pos.**, **Amino acid**, length, reference and **Z**, sorted from the
   largest deviation down. If there are no outliers, the table is not shown.
2. Click the row of an outlier you want to inspect, for example chain A, position 1 (or move
   to the row with the Tab key and press Enter). Bondamino switches to the tab **Graph** and
   scrolls to this amino acid, which is now selected.
3. Click **Reference** in the toolbar, so that the outlier bond is drawn in red.
4. Hold the mouse over the red bond to see its length, reference and Z-score.

The other outlier tables in **Statistics** (bond-angle outliers, Ramachandran outliers and
peptide bonds that are not trans) work the same way. You can also click a letter in the
sequence bar to jump to any amino acid.

**What can an outlier mean?** Be careful with conclusions. An outlier is a hint to look
closer, not a proof of an error. Common causes are:

- **An error in the model**, for example a wrongly placed atom, a mix-up of atom names or a
  poorly defined region. Check the electron density or the experimental data if you can.
- **Low resolution or old refinement.** Many outliers spread over the whole structure (as in
  `1LYZ`) point to the overall quality of the model rather than to individual errors.
- **Genuine strain.** A bond or angle can really be distorted, for example in an active site,
  near a bound metal, or in a covalently modified residue. Strong evidence is needed for such a
  claim: good resolution, clear density and ideally an independent structure that shows the
  same.
- **A wrong reference class.** The reference value may not fit the chemistry, for example for
  an unusual protonation state, a modification that the file does not declare, or a residue
  that is not a standard amino acid.
- **Alternate conformations.** If atoms of different conformations are combined, distances can
  look wrong. See [How do I look at alternate conformations?](#how-do-i-look-at-alternate-conformations).

### Step 7: Check the angles and the Ramachandran plot

**The Ramachandran plot**

1. Click the tab **Ramachandran**. You see φ (horizontal) against ψ (vertical) for every
   amino acid that has neighbours on both sides.
2. The **Legend** on the right shows the point colours: **Other amino acids**, **Glycine**,
   **Proline**, each with its count, and the **Outliers** (points with a red ring).
3. Hold the mouse over a point to see the residue, φ, ψ, the category and the region
   (**Favoured**, **Allowed** or **Outliers**).
4. Click a point. The residue is selected. When you go back to **Graph** or **Data**, the view
   jumps to it.
5. Use **Category** to show one category only: **General**, **Glycine**, **cis-proline**,
   **trans-proline**, **Pre-proline** or **Ile and Val**. The favoured and allowed regions of
   this category are then drawn as shaded areas. Each category has its own regions, so this is
   the correct way to judge a single point. Choose **All** to return.

The regions come from the Top8000 data set (Richardson Lab), the same basis as in MolProbity.

**The angle statistics**

1. Click **Statistics** and scroll to the section **Angles**.
2. **Ramachandran evaluation**: the number of **Assessed amino acids** and how many are
   **Favoured**, **Allowed** and **Outliers**, also per category. A list of all Ramachandran
   outliers follows; click a row to see the amino acid in the graph.
3. **cis and trans peptide bonds**: the number of trans, cis and twisted peptide bonds. Bonds
   that are not trans are listed. A cis bond that is not before a proline is highlighted in
   red, because it is rare.
4. **Deviation from reference values**: the same kind of rating as for bonds (step 5), now for
   bond angles.
5. **Bond angles**: all angles grouped by type (for example N–CA–C) with count, mean, standard
   deviation, minimum and maximum. Click a row to see every single angle.

For `1LYZ` 127 amino acids are assessed: 94 favoured, 22 allowed and 11 outliers.

**Interpretation.** A few Ramachandran outliers in a large protein are not unusual, and
glycines are often found in regions that are rare for other amino acids. Outliers in
well-defined parts of a high-resolution structure deserve a closer look. Clusters of outliers
in one region often point to a problem in that region of the model.

### Step 8: Look at the atom data

1. Click the tab **Data**. You see one row per atom of the selected model, in file order, with
   the columns **Record**, **No.**, **Atom**, **Alt**, **Residue**, **Chain**, **Pos.**,
   **iCode**, **x**, **y**, **z**, **Occupancy**, **B-factor**, **Element** and **Charge**.
   Hold the mouse over a column heading to see what it means.
2. The sequence bar is shown here too. Click a letter to scroll the table to that amino acid.
   An amino acid selected in the graph or in the Ramachandran plot is highlighted.

This view is useful to check the raw values behind a suspicious bond, for example a low
occupancy or a very high B-factor.

### Step 9: Export your results

All exports describe the **currently selected model and conformation**. The file name says
which one, for example `1LYZ_statistics_model-1.pdf`.

**A report for your records or a supplement (PDF)**

1. Click **Statistics**.
2. In the **Download** column on the right, decide whether to tick **Individual bonds**. With
   it, the PDF lists every single bond per atom pair. This is useful for small proteins, but
   for large structures it adds many pages.
3. Click **PDF**. The PDF contains the statistics and the Ramachandran tab: composition, bond
   and angle ratings, outlier lists, peptide bonds, disulfide bridges, the histogram and the
   Ramachandran plot.

**Tables for a spreadsheet (CSV)**

In the same column under **CSV**: **Atoms**, **Bonds**, **Bond angles** and **Dihedral angles**
(φ, ψ, ω per amino acid). See
[How do I prepare data for a spreadsheet?](#how-do-i-prepare-data-for-a-spreadsheet).

**All numbers in one machine-readable file (JSON)**

Click **JSON**. This file contains everything that the PDF contains, always including the
individual bonds. It is meant for your own scripts.

**Images**

The bond length histogram (**Statistics → Download → Histogram**) and the Ramachandran plot
(**Ramachandran → Download**) can be saved as PNG, JPEG or SVG. See
[How do I prepare a figure?](#how-do-i-prepare-a-figure).

### Step 10: Share and discuss with your work group

**Share the file**

If you did not share the file during upload: on the start page, click the row of the file
(not the file name) and then the pencil button (**Edit availability of the selected file**).
Tick your work group and click **Save**. All members of the group now see the file in their
list.

**Discuss with notes**

1. Open the file and click the tab **Notes**.
2. Write your note in the text field, for example: "Chain A, Lys 1: CA–CB is 2.02 Å (Z ≈ 49).
   Probably a coordinate error in this early model."
3. Under **Visible to**, choose **Only me** for a private note, or **For working groups** and
   then tick the group(s) that may read it. Only your groups in which the file is shared are
   offered, because nobody else could open the file. If there are none, you see "The file is not
   shared with any of your working groups."; share the file with your group first.
4. Click **Save note**. The note appears in the list below with your name and the date.

Your colleagues see the note when they open the file. On the start page, the column
**Latest note** shows the date of the newest note they may read, so new discussions are easy
to spot.

Notes are plain text. Write chain, position and atom names in the note so that others can find
the spot with the sequence bar. More in
[How do I write and manage notes?](#how-do-i-write-and-manage-notes).

### Step 11: Cite Bondamino

If you use results from Bondamino in a publication, a thesis or a poster, cite Bondamino and
the reference data it uses (Engh & Huber 2001, CCP4 Monomer Library, Top8000). The exact
citations are in [CITATION.md](CITATION.md).

In the methods section, also state:

- the PDB ID or the origin of your file,
- the model and conformation you analysed (they are in the export file name and at the top of
  the PDF),
- that bonds and angles were rated with |Z| > 4 as the outlier threshold, and
- that the Ramachandran regions are those of Top8000.

### Workflow checklist

- [ ] Method, resolution and R values checked (**Statistics → General**)
- [ ] Bond RMSZ and outliers checked (**Statistics → Bonds**)
- [ ] Interesting outliers inspected in the graph with **Reference** colouring
- [ ] Ramachandran plot and outliers checked, per **Category** where needed
- [ ] cis peptide bonds and chain breaks checked
- [ ] Right model and conformation selected before exporting
- [ ] PDF, CSV and images exported
- [ ] File shared with the work group, findings written as notes
- [ ] Bondamino and the reference data cited ([CITATION.md](CITATION.md))

---

## Task guides

### How do I compare the models of an NMR ensemble?

**When is this useful?** NMR structures contain several models. Each model is a valid
solution, and the differences between them show how well each part is defined. Example:
`1D3Z` (ubiquitin, 10 models).

1. Load `1D3Z`.
2. Open the model box at the top right (it shows **Model 1**). It lists **Model 1** to
   **Model 10** and **Average**.
3. Choose another model. All tabs are recalculated for this model. The selection boxes are
   grey for a moment while the data load.
4. In **Statistics → Bonds**, compare the RMSZ and the outliers between models. Write the
   values down or export a CSV per model.
5. In **Shortest and longest bond**, the rows under **All models** show the extreme bonds of
   the whole ensemble and the model they come from. The columns **Mean (Å)** and
   **Median (Å)** compare each bond with the same bond in all models; the number in brackets
   is the difference.

**The averaged model.** Choose **Average** to see the geometry averaged over all models: each
bond length, bond angle and torsion angle (φ, ψ, ω) is measured in every model and then
averaged. Use it to judge the typical geometry of the ensemble. Torsion angles are averaged on
the circle, so −179° and 179° give 180°, not 0°. The Ramachandran evaluation of the averaged
model uses the averaged φ and ψ. In `1D3Z`, all 74 assessed amino acids are favoured in the
averaged model, as in most single models (models 7 and 8 each have one amino acid in the
allowed region).

An average hides the spread between the models. In flexible regions, such as the C-terminal
tail of ubiquitin, the angles can differ a lot from model to model even if their mean looks
normal. There, also look at the individual models. The coordinates shown in the **Data** tab
for **Average** are plain means of the model coordinates; Bondamino does not calculate any
length or angle from them.

The B-factor column of NMR files does not describe thermal motion as in crystal structures,
and it is often zero or used for other values. Do not interpret B-factor colouring for NMR
structures in the same way.

### How do I look at alternate conformations?

**When is this useful?** In high-resolution crystal structures, some side chains or even
backbone segments are modelled in two or more positions. Example: `1EJG` (crambin, 0.54 Å),
with conformations `A`, `B` and `C` for 17 residues.

1. Load `1EJG`. The conformation box shows **Conformation: automatic**. In this mode, each
   residue uses the conformation with the highest occupancy.
2. Open the box and choose **Conformation A**, **Conformation B** or **Conformation C**.
   Residues that do not have the chosen conformation keep the one with the highest occupancy.
3. All tabs, including the statistics over all models, are recalculated.
4. In the **Data** tab, the column **Alt** shows which conformation each atom comes from, and
   **Occupancy** shows its share.

**Different amino acids in one position (microheterogeneity).** Sometimes the conformations
are different amino acids. In `1EJG`, position 22 is Pro in conformation A and Ser in B and C,
and position 25 is Leu in A and Ile in B and C. Such residues have a small bar at the bottom
of their letter in the sequence bar. The tooltip says, for example, "Amino acid per
conformation: PRO (A) · SER (B, C)". When you switch the conformation, the graph shows the
other amino acid.

**Interpretation.** Compare the RMSZ between conformations. If one conformation has many more
outliers, it may be less well modelled, often because it has a low occupancy and is therefore
less well defined by the data.

### How do I find chain breaks?

**When is this useful?** Residues that are missing from a model (often flexible loops) leave
gaps. Without knowing them, you might misread the sequence or the statistics.

1. Look at the sequence bar. A chain break is shown as a red vertical line between two
   letters. The tooltip says, for example, "Chain break: no peptide bond between SER 127 and
   THR 137 (distance 10.7 Å)".
2. In **Statistics → Bonds → Peptide bonds**, the row **Missing** (Chain breaks) counts
   all gaps. **Theoretically possible** is the number of pairs of consecutive amino acids in
   the same chain.

Example: `1HZH` has two chain breaks, both in chain K (between 127 and 137, and between 235
and 239).

Bondamino counts a peptide bond as missing if the C–N distance is too long for a bond. A gap
in the numbering without a large distance is not a chain break. Consecutive amino acids that
are directly connected but numbered with a jump are therefore not marked.

### How do I check disulfide bridges?

**When is this useful?** Disulfide bridges stabilise many extracellular proteins and
antibodies. A missing or strained bridge can point to a modelling problem or to reduced
cysteines.

1. In the graph, disulfide bridges are dotted brackets on the right side.
2. In the sequence bar, cysteines in a bridge have a small dot at the top left. The tooltip
   names the partner, for example "Disulfide bridge to CYS 127".
3. **Statistics → Bonds → Disulfide bridges** lists every bridge with chain, positions,
   length and **Source**:
   - **SSBOND record**: the file declares the bridge.
   - **Distance**: the file does not declare it, but the two sulfur atoms are close enough.

Examples: `1LYZ` has 4 bridges, `1EJG` 3, `1HZH` 15 (some of them between chains).

**Interpretation.** S–S bonds are usually close to 2.0 Å. Clearly longer distances for a
declared bridge (in `1LYZ` up to about 2.4 Å) suggest a less accurate model, radiation
damage, or a partly reduced bridge. Check the data before drawing conclusions.

### How do I find cis peptide bonds?

**When is this useful?** Almost all peptide bonds are trans. cis bonds occur mainly before
proline. A cis bond before another amino acid is rare and should be supported by clear data.

1. In the sequence bar, a cis peptide bond has its own marker between two letters. The tooltip
   says, for example, "cis peptide bond between TYR 140 and PRO 141 (ω = 0.0°)".
2. **Statistics → Angles → cis and trans peptide bonds** gives the counts and lists all
   peptide bonds that are not trans, with ω and **Type** (**cis**, **cis (before proline)** or
   **twisted**). A modified proline counts as proline. cis bonds that are not before a proline
   are highlighted. Click a row to see the peptide bond in the graph.

Example: `1HZH` has 10 cis peptide bonds, all before a proline.

### How do I look at ligands and modified amino acids?

**When is this useful?** To see how a cofactor or substrate is modelled, and whether
phosphorylated or otherwise modified amino acids are recognised. Example: `1ATP` (protein
kinase A with ATP, two Mn²⁺ ions, phosphothreonine 197 and phosphoserine 338).

**Ligands**

1. Load `1ATP`. Ligands, ions and water come after the protein chains, as in the file. Click
   one of the `?` cells at the end of the sequence bar (hold the mouse over a cell to see its
   name, for example "ATP 355") or scroll to the end of the graph.
2. ATP is drawn in its own band. Its layout follows the real 3D shape, flattened to 2D. Ligand
   bonds are dash-dotted.
3. Ligand bonds come from the CONECT records of the file or, if there are none, from the
   distances between atoms.
4. Ligand bonds have no reference value. In **Reference** colouring they are grey. Their
   lengths are included in **Bonds by atoms involved** and in the bond CSV.

Single ions such as Mn²⁺ have no bonds within the residue and are therefore not drawn. They
are listed in the **Data** tab. Bonds between a ligand and the protein (for example metal
coordination) are not shown.

**Modified amino acids**

1. In the sequence bar, modified amino acids have a small dot at the top right. The tooltip
   says, for example, "Modified amino acid, derived from Threonine".
2. They are shown with the one-letter code of the parent amino acid.
3. Backbone bonds are rated like the parent amino acid. For the side chain, Bondamino uses the
   reference values of the modified residue where the CCP4 Monomer Library provides them (in
   `1ATP`, TPO and SEP are rated).
4. In the composition tables they appear as, for example, "Threonine (modified)".

### How do I find flexible regions with the B-factor?

**When is this useful?** In crystal structures, high B-factors mark mobile or poorly ordered
regions, for example surface loops or chain ends. Outliers in such regions are less
surprising than in the well-ordered core.

1. Open the **Graph** tab and click **B-factor**.
2. Atoms and bonds are coloured from blue (low) to red (high). The **Legend** shows the scale
   in Å². The scale covers the middle 96 % of the values, so a few extreme atoms do not wash
   out the colours; values beyond the ends are marked with ≤ or ≥.
3. Scroll through the graph or use the sequence bar. Red stretches are the flexible regions.
4. Compare with **Reference** colouring: are the outliers in blue (ordered) or red (flexible)
   regions?

For NMR structures, see the note in
[How do I compare the models of an NMR ensemble?](#how-do-i-compare-the-models-of-an-nmr-ensemble).

### How do I work with a large structure?

**When is this useful?** Antibodies, complexes and other large structures have thousands of
atoms. Example: `1HZH` (IgG1 antibody), with protein chains H, K, L and M, carbohydrate chains
A and B, and about 10,600 atoms.

- Loading and switching models take longer than for small proteins. Wait until the selection
  boxes are active again.
- Use the sequence bar to jump instead of scrolling. Each chain starts with a label such as
  "H 1".
- **Statistics → Elements** and **Amino acids** have tables per chain (**Elements per chain**,
  **Composition per chain**, **Atoms per chain**, **Amino acids per chain**).
- For the PDF, leave **Individual bonds** unticked, otherwise the PDF becomes very long. Use
  the JSON or the bond CSV if you need every bond.
- For structures that are too large for the PDB format, see
  [What if a structure cannot be loaded?](#what-if-a-structure-cannot-be-loaded).

### How do I prepare a figure?

**Ramachandran plot**

1. Click **Ramachandran**.
2. Optional: choose a **Category** to show only this category, with its favoured and allowed
   regions shaded. With **All**, all points are shown without regions.
3. Under **Download**, click **SVG** (for editing in a vector program such as Inkscape or
   Illustrator, and for print), **PNG** (for slides and documents) or **JPEG**.

The image has a white background and fixed colours. It does not depend on the dark mode.

**Bond length histogram**

1. Click **Statistics**.
2. In the **Download** column under **Histogram**, click **SVG**, **PNG** or **JPEG**.

**The graph itself**

There is no export button for the graph. For a figure of the graph, set up the view (colour
mode, layout, zoom) and take a screenshot with the tool of your operating system. For print,
use the light display and zoom in so that the labels are sharp.

Before you publish a figure, check that the right model and conformation are selected, and
state them in the figure legend.

### How do I prepare data for a spreadsheet?

**When is this useful?** For your own plots, for statistics across several structures, or for
the supplementary material of a paper.

1. Select the model and conformation you want.
2. Click **Statistics** and, in the **Download** column under **CSV**, choose:
   - **Atoms**: the table of the **Data** tab.
   - **Bonds**: both atoms (chain, position, amino acid, atom name, atom number), type,
     length, reference, σ and Z.
   - **Bond angles**: all three atoms, angle type, value, reference, σ and Z.
   - **Dihedral angles**: φ, ψ and ω per amino acid with Ramachandran category and region.
3. Open the file in Excel, LibreOffice Calc or a similar program.

**Decimal separator.** The CSV format follows the language of the Bondamino interface:

| Interface language | Decimal separator | Column separator |
|---|---|---|
| English | point (`1.529`) | comma |
| German | comma (`1,529`) | semicolon |

If the numbers end up in one column, or appear as text or dates, your spreadsheet program
expects a different format. Either switch the Bondamino language to match your spreadsheet
program before exporting, or use the import function of the program (for example
**Data → From Text/CSV** in Excel) and set the separators there. The files are UTF-8 encoded,
so Å and Greek letters are kept.

To combine several structures, export one CSV per structure and add a column with the PDB ID
in your spreadsheet.

### How do I change who can see a file, or remove it?

Only the person who uploaded a file can change its availability or delete it.

**Change who can see a file**

1. On the start page, click the row of the file (not the file name). The row is highlighted.
2. Click the pencil button (**Edit availability of the selected file**). If the button stays
   grey, you are not the uploader of this file.
3. Tick or untick **For me personally** and the working groups, then click **Save**.

What **For me personally** means: as the uploader, you always see your own files while they
are shared with a group. The box decides what happens when the file is no longer shared with
any group. If it is ticked, the file stays available to you. If it is not ticked, the file is
deleted at that moment.

**Delete a file**

1. Select the row and click the red bin button (**Delete selected file**).
2. Confirm the question. The file and all notes on it are deleted for everybody, in all
   groups. This cannot be undone.

**Find files**

- The box next to the search field filters the list: **All available**, **My uploads** or one
  of your groups.
- **Search file name …** filters the list while you type.
- The column **Visible to** shows **personal** and the groups a file is shared with.

### What can a group leader do?

Each work group has a leader. Groups and members are managed by the administrator.

As a group leader you can:

- **Remove a file from your group.** On the start page, each group badge in the column
  **Visible to** has a small × (**Remove from group**) for the groups you lead. Click it and
  confirm. You do not need to be the uploader. If the file is then available nowhere (not
  personal and in no other group), it is deleted.
- **Remove a note from your group.** In the **Notes** tab, each group badge of a note has a ×
  for your groups. The note stays available to its author. A note without an author and
  without any group is deleted.

A group leader cannot edit or delete other people's notes, and cannot change the availability
of other people's files beyond removing them from the group.

**When is this useful?** To keep the group's list tidy, for example after a project ends or
when a file was shared with the wrong group.

### How do I write and manage notes?

1. Open a file and click **Notes**.
2. Write a note and choose **Visible to**: **Only me**, or **For working groups** plus the
   groups (only your groups in which the file is shared are offered). Click **Save note**.
3. To change a note you wrote, click its pencil button (**Edit note**), change the text or the
   visibility, and click **Save**. The note is then marked "edited".
4. To delete a note you wrote, click its bin button (**Delete note**) and confirm.
5. Use the box next to the heading **Notes** to filter: **All visible**, **Only me** (your
   notes that are not shared) or one group.

Only the author can edit or delete a note. Notes of other people are shown with their name.

**Ideas for useful notes**

- Your interpretation of an outlier, with chain, position and atoms.
- Which model or conformation you used for a figure.
- Open questions for the next group meeting.
- A comparison with another structure, for example "4LZT shows normal geometry here".

### How do I change my personal settings?

Click the round person icon at the top right. The account menu shows your name and:

- **Account**: opens your account page.
- **Language**: the available languages (only shown if the administrator has set up more than
  one). The page switches immediately.
- **Display**: **Light** or **Dark**. The page switches immediately.
- **Sign out**.

On the **Account** page you find:

- **User account**: your **Username**, **Name**, **Email** and **Working groups**. These are
  maintained by the administrator. Ask them to correct anything.
- **Change password**: enter your current password, then the new one twice, and click
  **Change password**. You stay signed in.
- **Settings**: **Language**, **Display** and **Time zone**. Click **Save**. Dates and times
  (upload time, latest note, notes) are then shown in your time zone. **Same as the server (…)**
  uses the server's default time zone, named in brackets. The **Language** field is missing if
  the administrator has not set up any language.

Your language, display and time zone are saved with your account and apply on every computer
you use.

Remember that the language also changes the number format of the CSV exports (see
[How do I prepare data for a spreadsheet?](#how-do-i-prepare-data-for-a-spreadsheet)).

### What if a structure cannot be loaded?

| Message | Cause | What to do |
|---|---|---|
| "Invalid PDB ID. A PDB ID consists of four characters, e.g. 1LYZ." | Typing error | Check the ID. It starts with a digit from 1 to 9. |
| "RCSB has no file in PDB format for …" | The ID does not exist, or the structure is only available as mmCIF | See below |
| "RCSB cannot be reached. Please try again later." | The server cannot reach the internet, or RCSB is down | Try later, or download the `.pdb` file yourself and upload it |
| "Only .pdb files are allowed." | Wrong file type (for example `.cif`, `.ent` or `.gz`) | Use a file in PDB format and give it the ending `.pdb` |
| "The file must not be larger than 20 MB." / "The file for … is larger than 20 MB." | Very large structure | See below |
| "The file could not be processed." | The file could not be read as a PDB file | Check that it is a complete PDB-format file, for example by opening it in a structure viewer |

**Structures that exist only as mmCIF.** Very large structures (for example ribosomes or
large virus capsids) do not fit into the classic PDB format, so RCSB offers them only as mmCIF.
Bondamino reads only the PDB format. If you need only part of such a structure, for example
one protein chain, you can save that part in PDB format with a structure viewer such as PyMOL
or UCSF ChimeraX and upload the resulting `.pdb` file.

---

## Troubleshooting

**I forgot my password.**
On the sign-in page, click **Forgotten your password?**, enter your email address and click
**Send link**. You receive an email "New password for Bondamino" with a link. The link is
valid for 24 hours and works only once. For privacy reasons, the page always says that a link
was sent, even if the address is not known. If no email arrives, check your spam folder and
ask your administrator which address is stored for you. Only a few emails per hour are sent
per address.

**The link in my email says "This link is invalid or has expired."**
The link was already used, is older than 24 hours, or your password was changed since. Click
**Request a new link** and enter your email address. This also works for an invitation link
that has expired.

**"Too many failed sign-in attempts. Please try again in 15 minutes."**
After five wrong passwords for the same username, sign-in is blocked for 15 minutes. (Many
wrong passwords from the same network address, for any username, have the same effect.) During
this time even the correct password is rejected. Wait 15 minutes, or use
**Forgotten your password?** to set a new password.

**I cannot see a file that a colleague uploaded.**
The file is not shared with any of your groups. Ask the uploader to share it with your group,
or ask the administrator to add you to the group.

**The pencil and bin buttons on the start page stay grey.**
First click a row to select it. If they stay grey, you are not the uploader of this file. Only
the uploader can change or delete it.

**The model or conformation box is grey.**
The file has only one model or no alternate conformations. The boxes are also grey for a
moment while a new model loads.

**Nothing happens when I choose another model, or a download does not start.**
Wait a few seconds; large structures need time. If loading or an export fails, a red message
appears at the top right, for example "The model or conformation could not be loaded." or "The
PDF could not be created.", and the model box returns to the model still shown. Reload the page
and try again. If the problem persists, tell your administrator which file and which model you
used.

**Dates and times are off by some hours.**
Choose your time zone under **Account → Settings → Time zone** and click **Save**.

**Bondamino does not remember my colour and layout choice.**
Your browser blocks local storage for this site (for example in some private windows).
Everything works, but the choice is not remembered. Allow site data for the Bondamino address.

**The CSV opens with all values in one column.**
See the decimal separator table in
[How do I prepare data for a spreadsheet?](#how-do-i-prepare-data-for-a-spreadsheet).

**The upload was rejected.**
See [What if a structure cannot be loaded?](#what-if-a-structure-cannot-be-loaded).

**Something else does not work.**
Contact your Bondamino administrator. Tell them the file name, the tab and what you did. For
server problems, administrators find help in [INSTALL.md](INSTALL.md).
