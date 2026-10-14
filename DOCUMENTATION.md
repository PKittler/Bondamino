# Bondamino documentation

This document describes every page, button, option and export of the Bondamino web interface,
and the scientific definitions Bondamino uses. It is a reference, not a tutorial.

- For an overview of the project, see [README.md](README.md).
- For installation and server configuration, see [INSTALL.md](INSTALL.md).
- For step-by-step workflows, see [TUTORIAL.md](TUTORIAL.md).
- For how to cite Bondamino and its data sources, see [CITATION.md](CITATION.md).
- For licences and transparency information, see [LICENSE.md](LICENSE.md) and [TRANSPARENCY.md](TRANSPARENCY.md).

Labels of buttons, tabs and options are quoted as the British English interface shows them,
in bold, for example **Statistics**. The interface is also available in German.

## Contents

- [Concepts](#concepts)
  - [Users and administrators](#users-and-administrators)
  - [Working groups and group leaders](#working-groups-and-group-leaders)
  - [Who can see and change what](#who-can-see-and-change-what)
  - [Personal and group availability of files](#personal-and-group-availability-of-files)
  - [What happens when something is deleted](#what-happens-when-something-is-deleted)
- [Pages](#pages)
  - [Navigation bar and account menu](#navigation-bar-and-account-menu)
  - [Sign-in page](#sign-in-page)
  - [Forgotten password](#forgotten-password)
  - [Invitation and setting a password](#invitation-and-setting-a-password)
  - [Start page: file list](#start-page-file-list)
  - [Upload dialog](#upload-dialog)
  - [Edit availability page](#edit-availability-page)
  - [Detail page: header](#detail-page-header)
  - [Sequence bar](#sequence-bar)
  - [Graph tab](#graph-tab)
  - [Ramachandran tab](#ramachandran-tab)
  - [Data tab](#data-tab)
  - [Statistics tab](#statistics-tab)
  - [Exports](#exports)
  - [Notes tab](#notes-tab)
  - [Account page](#account-page)
  - [Administration](#administration)
- [Limits and behaviour worth knowing](#limits-and-behaviour-worth-knowing)
- [Scientific background](#scientific-background)
  - [Reading the PDB file](#reading-the-pdb-file)
  - [Models and the averaged model](#models-and-the-averaged-model)
  - [Alternate conformations](#alternate-conformations)
  - [How bonds are determined](#how-bonds-are-determined)
  - [Bond angles](#bond-angles)
  - [Reference values](#reference-values)
  - [Z-score, RMSZ, warnings and outliers](#z-score-rmsz-warnings-and-outliers)
  - [Dihedral angles](#dihedral-angles)
  - [Peptide bond conformation](#peptide-bond-conformation)
  - [Ramachandran evaluation](#ramachandran-evaluation)
  - [Statistics over all models](#statistics-over-all-models)
  - [Colour scales](#colour-scales)
  - [References](#references)
- [Glossary](#glossary)

---

## Concepts

### Users and administrators

- **Users** sign in with a username and a password. There is no self-registration: an
  administrator creates every account. The new user receives an invitation email with a link
  and sets the password there.
- **Administrators** are users with access to the administration area at `/admin/` (in Django
  terms: users with "staff status"; a "superuser" has all permissions). The first administrator
  is created during installation, see [INSTALL.md](INSTALL.md). Administrators see and can change
  all files, notes, users and working groups, regardless of the rules below.

### Working groups and group leaders

- A **working group** is a set of users, for example a lab or a project. Administrators create
  working groups and assign members.
- Every working group has exactly one **group leader**. The leader is always also a member
  (the administration adds the leader automatically when the group is saved).
- Files and notes can be shared with one or several working groups. All members of a group can
  then see them.
- The group leader may remove any file or note from the group, even one the leader did not
  upload or write.

### Who can see and change what

| Action | Who may do it |
|---|---|
| Upload a file | Every signed-in user |
| See a file, open its detail page, download it, use all tabs and exports | The uploader, and every member of a working group the file is shared with |
| Change a file's availability | Only the uploader |
| Delete a file (for everybody) | Only the uploader |
| Remove a file from a working group | Only the leader of that group |
| Write a note on a file | Everybody who can see the file |
| Read a note | The author, and the members of the working groups the note is shared with (they must also be able to see the file) |
| Edit or delete a note | Only the author |
| Remove a note from a working group | Only the leader of that group |
| Everything | Administrators, in the administration area |

Files and notes you are not entitled to see do not appear anywhere. If you open the address of
such a file directly, Bondamino answers "not found", so nobody can find out which files exist.

### Personal and group availability of files

When you upload a file, you choose under **Available to**:

- **For me personally**: the file stays available to you even when it is not shared with any group.
- **Working group …**: the file is shared with that group.

At least one of these must be selected. As long as a file exists, its uploader always sees it.
The **personal** setting matters when the file is no longer shared with any group: a file that
is neither personal nor shared with a group is deleted automatically. This can happen when a
group leader removes the file from the last group it was shared with (Bondamino then reports that
the file was removed and deleted).

### What happens when something is deleted

| Deleted | Consequence |
|---|---|
| A file | The stored file and all notes on it are deleted. |
| A working group | Files that were available only through this group (not personal, no other group) are deleted. Notes that have no author any more and were readable only by this group are deleted. Everything else stays, just no longer shared with the deleted group. |
| A user | Their files that are not shared with any group are deleted. Their files that are shared with groups stay in those groups; they lose the personal setting and no longer have an uploader. Their notes that are not shared with any group are deleted; shared notes stay, shown without author ("–"). |

A file or note without uploader or author can no longer be edited or deleted by normal users;
group leaders can still remove it from their groups (it is then deleted once it is in no group),
and administrators can change or delete it. A user who leads a working group cannot be deleted
until another leader is assigned.

---

## Pages

Every page except sign-in, forgotten password and the password links requires you to be
signed in; otherwise you are sent to the sign-in page.

### Navigation bar and account menu

The bar at the top shows **Bondamino** on the left (a link to the start page).

When you are signed in, a round button with a person icon (**User account**) on the right opens
the account menu:

| Entry | Function |
|---|---|
| Your name and username | Information only |
| **Account** | Opens the [account page](#account-page) |
| **Language** | One entry per available language; the selected one is ticked. Choosing one switches the interface and saves the choice in your settings. Shown only if the administrator has set up more than one language. |
| **Display**: **Light**, **Dark** | Switches between light and dark colours immediately and saves the choice in your settings. |
| **Log out** | Signs you out and returns to the sign-in page. |

When you are not signed in, the bar shows a moon button (**Toggle dark mode**, remembered only in
this browser) and, if more than one language is set up, a language selector.

Bondamino chooses the interface language in this order: your saved setting, the language last
chosen in this browser, your browser's preferred languages, German, the first language set up.

### Sign-in page

Fields **Username** and **Password**, button **Log in**, and the link **Forgotten your password?**.
After signing in you reach the start page (or the page you originally wanted to open).

Failed attempts are limited (see [Limits](#limits-and-behaviour-worth-knowing)). While a lock is
active, even the correct password is rejected with the message "Too many failed sign-in attempts.
Please try again in 15 minutes."

### Forgotten password

Enter your **Email** and click **Send link**. Bondamino always shows the same confirmation ("If
a user account belongs to this email address, a link … has been sent"), whether or not an account
exists, so that nobody can test which addresses are registered. The link in the email is valid
for 24 hours and works once. This also works for invited users who never set a password or whose
invitation has expired. **Back to sign-in** returns to the sign-in page.

### Invitation and setting a password

The link from the invitation email opens the page **Set password**; the link from a
forgotten-password email opens **New password**. Both show your username and two fields,
**New password** and **New password confirmation**, and the button **Save password**. The password
must not be too similar to your user data, must have the minimum length set on the server (by
default 8 characters), must not be a common password and must not consist of digits only. After
saving you are signed in automatically and see the start page with "Your password has been saved."

If the link is older than 24 hours or has already been used (also after the password was changed
or the user signed in by other means since the link was sent), the page
shows "This link is invalid or has expired." with the button **Request a new link**.

### Start page: file list

The start page lists all files you may see, newest upload first.

**Toolbar** (from left to right):

| Element | Function |
|---|---|
| Button with a "+" document icon (**New PDB file**) | Opens the [upload dialog](#upload-dialog). |
| Pencil button (**Edit availability of the selected file**) | Opens the [edit availability page](#edit-availability-page) for the selected file. Active only if you uploaded the selected file. |
| Bin button (**Delete selected file**) | Deletes the selected file permanently, for everybody, after a confirmation ("Permanently delete "…"?"). Active only if you uploaded the selected file. |
| **Filter** selector | **All available** (everything you may see), **My uploads** (files you uploaded), or the name of one of your working groups (files shared with it). |
| **Search file name …** | Hides all rows whose file name does not contain the text (not case-sensitive) while you type. Pressing Enter applies the search on the server and keeps it in the page address. |

Select a row by clicking it (or with the keyboard: Tab to the row, then Enter or Space).
Clicking the file name opens the [detail page](#detail-page-header).

**Columns:**

| Column | Content |
|---|---|
| **File** | File name, link to the detail page |
| **Visible to** | A grey badge **personal** if the file is personal, and one blue badge per working group it is shared with. If you lead that group, the badge has a small "×" (**Remove from group**); after a confirmation the file is removed from the group (and deleted if it is then available nowhere). |
| **Uploaded** | Date and time of the upload |
| **Latest note** | Date and time of the newest note on this file that you may read, or "–" |

If nothing matches, the list shows "No PDB files yet." or "No PDB files for this selection."

### Upload dialog

Title **New PDB file**. You provide the structure in one of two ways, not both:

1. **Local file**: drag a file onto the dashed area ("Drag a file here or") or click
   **Choose file**. The chosen file name replaces the hint **Supports .pdb files**.
2. **Load PDB ID from RCSB**: type a four-character PDB ID, for example `1LYZ`. Bondamino downloads
   the structure in PDB format from the RCSB Protein Data Bank and stores it as `1LYZ.pdb`.

Choosing a file clears the PDB ID field, and typing an ID clears the chosen file.

Under **Available to** select **For me personally** (selected by default) and/or **Working group …**
for each of your groups (see [availability](#personal-and-group-availability-of-files)). If you
are in no working group, the dialog says so.

**Cancel** closes the dialog. **Analyse protein** uploads (or downloads) the file; the button shows
a spinner while this runs. Bondamino then reads the file once as a test. On success the detail
page opens. On failure the dialog opens again with the error, and nothing is stored.

Possible error messages:

| Message | Cause |
|---|---|
| Only .pdb files are allowed. | The file name does not end in `.pdb` |
| The file must not be larger than 20 MB. | Uploaded file too large |
| Please either choose a file or enter a PDB ID, not both. / Please choose a file or enter a PDB ID. | Neither or both sources given |
| Invalid PDB ID. A PDB ID consists of four characters, e.g. 1LYZ. | The ID is not a digit 1–9 followed by three letters or digits |
| RCSB has no file in PDB format for … | Unknown ID, or a large structure that RCSB offers only in mmCIF format |
| RCSB refused the download of … (HTTP …) / RCSB cannot be reached. Please try again later. / The file for … is larger than 20 MB. | Download problems |
| The file must be available to you personally or to at least one working group. … | Nothing selected under **Available to** |
| The file could not be processed. | The test read failed; the file was discarded |

### Edit availability page

Opened with the pencil button on the start page (only for your own uploads). Title
**Edit availability**, followed by the file name. The same **Available to** options as in the upload
dialog. The working groups offered are your current groups plus any group the file is already
shared with (so you can withdraw a file from a group you have left). **Save** stores the change and
returns to the start page; **Cancel** returns without saving. At least one option must remain
selected; to withdraw a file completely, delete it instead.

### Detail page: header

Clicking a file name opens its detail page. Bondamino reads the file anew each time. If it cannot
be processed, you return to the start page with the message "The file "…" could not be processed."

The top row contains:

| Element | Function |
|---|---|
| Tabs **Graph**, **Ramachandran**, **Data**, **Statistics**, **Notes** | Switch between the views described below. |
| **File** and the file name | Information |
| Download button (**Download PDB file**) | Downloads the original file under its original name. |
| **Model** selector | One entry **Model 1**, **Model 2**, … per model in the file, plus **Average** (the [averaged model](#models-and-the-averaged-model)). The first model is shown initially. If the file has only one model, the selector is disabled ("The file contains only one model."). |
| **Conformation** selector | **Conformation: automatic** and one entry per alternate location identifier in the file (**Conformation A**, **Conformation B**, …). See [alternate conformations](#alternate-conformations). If the file has none, the selector shows **Single conformation** and is disabled. |

Changing the model or conformation reloads the data and redraws all tabs. While loading, both
selectors are disabled. Changing the conformation also recomputes the [statistics over all
models](#statistics-over-all-models), because they depend on it. If loading fails, the display
stays unchanged.

The structure information from the file header (method, resolution, R values) is shown in the
**Statistics** tab under **General**.

### Sequence bar

Shown above the **Graph** and **Data** tabs. One cell per residue in file order, with its one-letter
code.

| Mark | Meaning |
|---|---|
| Letter | One-letter code; a modified amino acid shows the code of its parent amino acid |
| Grey "?" | Not an amino acid: ligand, ion or water |
| Small gap between cells | A new chain starts |
| Number above a cell | Chain identifier and number at the start of each chain (for example "A 1"), then every residue number divisible by 10 (dropped where labels would overlap) |
| Red vertical line between two cells | Chain break: no peptide bond between these residues. Hover for the residue names and the C–N distance. |
| Dotted vertical line between two cells | cis peptide bond. Hover for the residues and ω. |
| Green dot, top right | Modified amino acid |
| Dot in sulfur colour, top left | Cysteine in a disulfide bridge |
| Purple bar at the bottom | The alternate conformations of this position are different amino acids (microheterogeneity) |
| Light highlight | Residues currently visible in the graph or in the data table |
| Strong highlight | The selected residue |

Hovering over a cell shows residue name, number and chain, the parent amino acid of a modified
residue ("Modified amino acid, derived from …"), the amino acid per conformation, and the disulfide
partners.

Clicking a cell **selects** the residue: it is highlighted in the sequence bar, in the graph (its
label), in the data table (its rows) and in the Ramachandran plot (its point), and the visible view
scrolls to it. The other view catches up when you switch to it. The bar scrolls along with the
graph or table.

### Graph tab

The graph shows every residue as a horizontal band, stacked from top to bottom in file order.
Within an amino acid the backbone atoms N, CA and C form a vertical column; O (and OXT at the
C-terminus) sit to the left; the side chain grows to the right from CA as a branching tree, with
rings (His, Phe, Tyr, Trp, and Pro in the proportional layout) drawn as closed polygons. The left
margin of each band shows the one-letter code, the three-letter code and the residue number
(for example "#52"). Peptide bonds connect C of one band with N of the next.

**Ligands** get a band of their own, laid out from their 3D coordinates (projected onto their two
main directions, then untangled so atoms do not overlap). **Disulfide bridges** connect distant
residues; they are drawn as brackets in a separate column on the right, one lane per bridge.
**Water** appears as a single oxygen atom. Hydrogen atoms, single-atom ions and ligand atoms
without any bond are not drawn (an ion only gets an empty band with its label), but they appear in
the [data table](#data-tab) and in the statistics.

Every atom is a circle labelled with its atom name and coloured by element. Every bond shows its
length in Å next to it (hidden where the bond is drawn very short).

**Line styles** (see also the legend):

| Line | Bond type |
|---|---|
| Solid | Backbone and peptide bond |
| Dashed | Side chain (including CA–CB) |
| Dash-dot | Bond within a ligand |
| Dotted bracket | Disulfide bridge |

**Toolbar**

| Button | Function |
|---|---|
| **By length** | Bond colour from red (shortest bond shown) to green (longest bond shown). Disulfide bridges are excluded from the range. |
| **Uniform** | All bonds in one neutral colour |
| **Reference** | Bond colour by the deviation from the reference value, in standard deviations (see [Z-score](#z-score-rmsz-warnings-and-outliers)): green for \|Z\| ≤ 2, orange for 2 < \|Z\| ≤ 4, red for \|Z\| > 4 (outlier), grey for bonds without reference value |
| **B-factor** | Atoms and bonds coloured by B-factor, from blue (low) via violet to red (high); a bond gets the mean B-factor of its two atoms (see [colour scales](#colour-scales)) |
| **Schematic** | Fixed grid: every residue band has the same height and all bonds have the same drawn length, independent of the real bond lengths |
| **Proportional** | Drawn bond lengths are proportional to the real bond lengths (55 pixels per Å). This holds exactly for the backbone, for unbranched side-chain segments and for rings; at branch points the branches are spread vertically, so these bonds are drawn somewhat longer. |
| Circular arrow (**Reset zoom**) | Returns to the original zoom level and position |

Switching style or layout animates the change. Bondamino remembers your style and layout in this
browser.

**Navigating the graph**

| Action | Effect |
|---|---|
| Mouse wheel | Scrolls up and down |
| Shift + mouse wheel | Scrolls left and right |
| Ctrl + mouse wheel, or pinching on a touchpad | Zooms (from 0.1× to 4×) |
| Dragging with the mouse | Moves the view |

**Tooltips** (hold the mouse pointer still over an element):

- **Atom**: atom name, residue, number and chain; the B-factor (if present); every bond angle with
  this atom at its centre, with its value and, if a reference exists, its Z-score.
- **Bond**: the two atoms, the length in Å, the bond type, and the reference value with standard
  deviation and Z-score ("Reference 1.459 ± 0.020 Å · Z = −0.1"), or "No reference value".

**Legend** (right side):

- The line styles **— Backbone**, **- - Side chain**, and, if present, **··· Disulfide bridge** and
  **-·- Ligand**.
- The element colours of the atoms in the model (C, N, O, S, H, P first, then alphabetically).
  Hidden in **B-factor** mode.
- In **Reference** mode: **Deviation |Z|** with ≤ 2, 2–4, **> 4 (outlier)** and **no reference**.
- In **B-factor** mode: **B-factor (Å²)** with a colour bar and the values at both ends (a "≤" or
  "≥" in front means that more extreme values exist and are shown in the end colour), and
  **without B-factor**. If no atom has a B-factor, it says **no B-factors**.

There is no image export of the graph itself.

### Ramachandran tab

A Ramachandran plot: one point per amino acid that has both φ and ψ (see
[dihedral angles](#dihedral-angles)), with φ on the horizontal and ψ on the vertical axis, each from
−180° to 180°. Point colours: blue for **Other amino acids**, orange for **Glycine**, green for
**Proline** (modified amino acids count as their parent amino acid). Outliers have a red ring and
are drawn on top.

Hovering over a point shows residue, number, chain, φ, ψ, the category and the region (for example
"General · Favoured"). Clicking a point selects the residue (see [sequence bar](#sequence-bar));
the graph and data table jump to it when you switch to them.

The side panel contains:

| Element | Function |
|---|---|
| **Category** | **All** shows every point. **General**, **Glycine**, **cis-proline**, **trans-proline**, **Pre-proline** or **Ile and Val** shows only the amino acids of that category, together with that category's **Favoured** (darker) and **Allowed** (lighter) regions as shaded areas (see [Ramachandran evaluation](#ramachandran-evaluation)). |
| **Legend** | Point colours with counts, **Outliers** with count, and the region shading when a category is selected |
| **Download**: **PNG**, **JPEG**, **SVG** | Saves the plot as an image, see [Exports](#exports) |

The numerical Ramachandran evaluation is in the **Statistics** tab under **Angles**.

### Data tab

A table with one row per atom of the displayed model and conformation, in file order. A line
separates residues; HETATM rows are shown in grey; the rows of the selected residue are highlighted.
Hovering over a column heading explains it.

| Column | Content |
|---|---|
| **Record** | ATOM or HETATM |
| **No.** | Atom serial number |
| **Atom** | Atom name |
| **Alt** | Alternate location identifier |
| **Residue** | Residue name |
| **Chain** | Chain identifier |
| **Pos.** | Residue sequence number |
| **iCode** | Insertion code |
| **x**, **y**, **z** | Coordinates in Å (3 decimals) |
| **Occupancy** | Occupancy (2 decimals) |
| **B-factor** | Temperature factor in Å² (2 decimals) |
| **Element** | Element symbol |
| **Charge** | Formal charge |

For the **Average** model, coordinates, occupancy and B-factor are the averages over all models.

### Statistics tab

All numbers refer to the displayed model and conformation, except the rows marked **All models**
and the **Average (Å)** and **Median (Å)** columns of **Shortest and longest bond**, which cover all
models.
Tables with a small arrow in a row can be expanded by clicking the row; they stay expanded when you
change the model.

**General**

- **Elements**: count and share of each element; **Total atoms**.
- Key figures: **Models** (number of models in the file), and, if given in the file header,
  **Method** (for example X-ray crystallography, Solution NMR), **Resolution** (Å), R-work
  and R-free. R-work is taken from "R VALUE (WORKING SET)" or, failing that,
  "R VALUE (WORKING + TEST SET)" in REMARK 3.
- **Amino acids**: count and share of each amino acid (modified ones as "… (modified)");
  **Total amino acids**. Ligands, ions and water are not counted.

**Elements**

- **Elements per chain**: count of each element in each chain, and its share within that chain.
- **Atoms per chain**: number of atoms per chain and share of all atoms.

**Amino acids**

- **Composition per chain**: count of each amino acid per chain and its share within that chain.
- **Amino acids per chain**: number of amino acids per chain and share of all.

**Bonds**

- **Deviation from reference values**:
  - Key figures **Assessed bonds** (bonds with a reference value, "… of … (…%)"), **RMSZ**, and
    **Outliers (|Z| > 4)**.
  - A table with one row for the **Backbone** and one row per amino acid for its side chain, each
    with **Count**, **RMSZ** and **Outliers**. Click a row to show its groups: **Bond** (atom pair),
    **Class** (see [reference values](#reference-values)), **Reference (Å)** (mean ± standard
    deviation), **Count**, **Average (Å)** of the bonds in the structure, **Dev. (Å)** (average
    minus reference), **RMSZ**, **Outliers**. RMSZ values above 2 are shown in orange, above 4 in
    red; non-zero outlier counts in red.
  - The sources of the reference values.
  - **Outliers (|Z| > 4)**: every outlier bond, sorted by decreasing \|Z\|, with **Atoms**,
    **Chain**, **Pos.**, **Amino acid**, **Length (Å)**, **Reference (Å)** and **Z**.
- **Peptide bonds**: **Existent** and **Non-existent** (**Chain breaks**) relative to
  **Theoretically possible** (pairs of consecutive amino acids of the same chain in the file).
- **Disulfide bridges** (only if present): **Chain**, **Pos.**, **Length (Å)** and **Source**:
  **SSBOND record** or **Distance** (see [disulfide bridges](#disulfide-bridges)).
- **Shortest and longest bond**: for **All models** (with the model it occurs in) and for the
  **Selected model**: **Model**, **Length (Å)**, **Atoms**, **Chain**, **Pos.**, **Atom no.**,
  **Amino acid**, and the **Average (Å)** and **Median (Å)** of that same bond over all models,
  each followed by the difference of this bond from it in parentheses.
- **Distribution of bond lengths**: a histogram of all bond lengths of the model in bins of 0.01 Å;
  the vertical axis is the share of all bonds. Hover over a column for the bin range, the share and
  the number of bonds.
- **Bonds by atoms involved**: one row per atom pair (for example N–CA, CB–CG), with **Chains**,
  **Count**, **Average (Å)**, **Median (Å)**, **Min (Å)**, **Max (Å)** and **Δ Max − Min (Å)**;
  **Total bonds** at the bottom. Click a row to list every single bond with **Length (Å)**,
  **Dev. from avg. (Å)**, **Chain**, **Pos.**, **Atom no.** and **Amino acid**; the shortest and
  longest are marked **Min** and **Max**.

**Angles**

- **Ramachandran evaluation**: **Assessed amino acids**, **Favoured**, **Allowed** and **Outliers**
  with counts and shares; a table per **Category**; the source of the regions; and a list of
  **Outliers** with **Chain**, **Pos.**, **Amino acid**, **Category**, φ and ψ.
- **cis and trans peptide bonds**: **Peptide bonds**, **trans**, **cis** (with the number before a
  proline) and **twisted** (see [peptide bond conformation](#peptide-bond-conformation)). If any
  peptide bond is not trans, a table lists them with **Chain**, **Pos.**, **Amino acids**, ω and
  **Type** (**cis**, **cis (before proline)** or **twisted**). A cis bond that is not before a proline
  is highlighted in red, because it is unusual.
- **Deviation from reference values**: the same structure as for bonds, for bond angles in degrees
  (**Assessed angles**, **Angle**, **Reference (°)**, **Mean (°)**, **Dev. (°)**).
- **Bond angles**: one row per angle type (for example N–CA–C) with **Count**, **Mean (°)**, **SD (°)**,
  **Min (°)** and **Max (°)**; **Total bond angles** at the bottom. Click a row to list the individual
  angles with **Angle (°)**, **Dev. from mean (°)**, **Chain**, **Pos.**, **Atom no.**, **Amino acid**
  and **Min**/**Max** marks.

The **Download** panel on the right contains the exports described next.

### Exports

All exports are created in your browser from the displayed model and conformation. File names have
the form `<file name>_<content>_<model>`, for example `1LYZ_statistics_model-1.pdf` or
`1EJG_bonds_model-1-conformation-a.csv`.

**Statistics tab, Download panel**

| Button | Content |
|---|---|
| **PDF** | A printable A4 report of the statistics (see below). Option **Individual bonds**: also list every single bond per atom pair (this can run to many pages for large proteins). Bondamino remembers this option in your browser. |
| **JSON** | All values of the statistics and angle sections in a machine-readable text format. Always contains every individual bond. |
| CSV **Atoms** | The data table: same columns as the [Data tab](#data-tab) |
| CSV **Bonds** | One row per bond: chain, position, amino acid, atom name and atom number of both atoms, **Type**, **Length (Å)**, **Reference (Å)**, **σ (Å)**, **Z** |
| CSV **Bond angles** | One row per bond angle: the three atoms (as above), **Type** (for example N–CA–C), **Angle (°)**, **Reference (°)**, **σ (°)**, **Z** |
| CSV **Dihedral angles** | One row per amino acid: **Chain**, **Pos.**, **Amino acid**, φ, ψ, ω (the peptide bond before this amino acid), **Category** and **Region** of the Ramachandran evaluation |
| **Histogram**: **PNG**, **JPEG**, **SVG** | The bond length histogram as an image (1200 × 480, raster images at double resolution), always with light colours and white background |

CSV files use UTF-8 encoding. Their number format follows the interface language so that a
spreadsheet program in the same language opens them directly: in English, columns are separated
by commas and decimals use a point; in German, columns are separated by semicolons and decimals use
a comma. Lengths have 3 decimals, angles 1 decimal, Z-scores 2 decimals.

**Ramachandran tab, Download panel**: **PNG**, **JPEG**, **SVG** save the plot (800 × 800, raster
images at double resolution) with light colours and white background. The selected **Category**
applies: with a category selected, only its points and its regions are exported.

**PDF report contents**, in this order: title with file name, model and creation date; **General**
(models, amino acids, atoms, chains, bonds; method, resolution, R-work, R-free if known; elements;
amino acids); elements and amino acids per chain; **Bonds** (deviation from reference values with
outliers, peptide bonds, disulfide bridges, shortest and longest bond, histogram, bonds by atoms
involved, optionally the individual bonds); **Angles** (Ramachandran plot of all categories with
legend, Ramachandran evaluation and outliers, cis and trans peptide bonds, bond angles, deviation
from reference values). Each page has a footer with file, model and page number. Because the
built-in PDF fonts lack Greek letters, the PDF writes "phi", "psi", "omega" and "Delta".

### Notes tab

Notes are short texts on a file, for example observations or to-dos.

**Writing a note**: type into the text field ("Write a note …", up to 10,000 characters), choose
**Visible to**: **Only me** or **Working groups**. With **Working groups**, tick one or more
**Working group …** boxes. Then click **Save note**. The new note appears in the list without
reloading the page.

**List of notes**, newest first, with a filter: **All visible**, **Only me** (notes not shared with
any group, that is your personal notes) or a working group (notes readable by that group). Each note
shows:

- the author's name ("–" if the author's account was deleted), date and time, and **edited** if the
  note was changed (hover for the time of the last change);
- one blue badge per working group that may read it, or a grey **personal** badge;
- for leaders of a group: a small "×" on that group's badge (**Remove from group**), with
  confirmation; the note remains available to its author;
- for the author: a pencil button (**Edit note**) and a bin button (**Delete note**, with the
  confirmation "Permanently delete the note?").

**Edit note** turns the text into an editable field with the visibility options; **Save** stores the
change, **Cancel** discards it.

Notes are deleted together with their file. Only notes you are allowed to read are counted in the
**Latest note** column of the start page.

### Account page

Opened via **Account** in the account menu.

- **User account**: **Username**, **Name**, **Email** and your **Working groups**. These are
  maintained by the administrator and cannot be changed here.
- **Change password**: **Old password**, **New password**, **New password confirmation**, button
  **Change password**. The same password rules apply as when
  [setting a password](#invitation-and-setting-a-password). You stay signed in.
- **Settings**: **Language** and **Display** (**Light** or **Dark**), button **Save**. These are the
  same settings as in the account menu.

### Administration

The administration area at `/admin/` is available to administrators only. Its sign-in is limited
in the same way as the normal sign-in. Installation, the guided setup and server settings are
described in [INSTALL.md](INSTALL.md).

**Users**

- **Add user**: username, email address (required), first and last name, and working group
  memberships. No password is entered: after saving, Bondamino emails an invitation with a link to
  set the password (valid 24 hours, usable once). If sending fails, the user is still created and an
  error message names the problem.
- The user list has a column **Password set**, which shows whether the user has accepted the
  invitation.
- The action **Resend invitation** (select users, choose the action, click "Go") sends a new
  invitation to each selected user who has no password yet, has an email address and is active.
  Other selected users are skipped and counted in a warning.
- The page of an existing user also shows its **Working groups** and its **Settings**
  (**Language**, **Display**).
- Deleting a user triggers the clean-up described in
  [What happens when something is deleted](#what-happens-when-something-is-deleted).
- The Django "Groups" in the same section are permission groups for administrators. They are not
  working groups.

**PDB Grapher** section

| Entry | Function |
|---|---|
| **Languages** | The languages users can choose. Each has a **Language** code and a **Display name** written in the language itself (for example "English"). Translations exist for German (`de`) and British English (`en-gb`). Only languages entered here can be selected. |
| **Working groups** | **Name**, **Description**, **Group leader** and **Members**. The list shows the number of members and can be searched by group name, leader or member. The leader is added to the members automatically. |
| **PDB files** | All files, independent of visibility: **File**, **File name**, **Uploaded by**, **Uploaded on**, **Personal**, **Working groups**. Filter by working group and personal; search by file name or uploader. |
| **Notes** | All notes: **PDB file**, **Author**, **Text**, **Readable by** (working groups), **Created on**, **Edited on**. Filter by working group; search by text, file name or author. |

---

## Limits and behaviour worth knowing

| Topic | Behaviour |
|---|---|
| File format | Only the classic PDB format. The file name must end in `.pdb` (not case-sensitive). mmCIF/PDBx files (`.cif`) are not supported. |
| File size | At most 20 MB (20 × 1024 × 1024 bytes), for uploads and for RCSB downloads. |
| PDB IDs | Four characters: a digit 1–9 followed by three letters or digits. Very large structures that RCSB provides only as mmCIF cannot be loaded. The download times out after 20 seconds. |
| Unreadable files | A new upload that cannot be read is discarded. If an existing file cannot be read, the detail page returns to the start page with an error message. Lines that are not understood are skipped. |
| Processing time | Files are read anew whenever the detail page opens; nothing derived is stored. On opening, the statistics over all models are computed, which takes longer for files with many models. |
| Sign-in limit | After 5 failed attempts for the same username from the same network address, or 20 failed attempts in total from one network address, sign-in from that address is locked for 15 minutes. A successful sign-in resets only the counter for that username. The same applies to the administration sign-in. |
| Forgotten password limit | At most 3 requests per email address and 10 per network address within one hour. Further requests show the usual confirmation, but no email is sent. |
| Password links | Invitation and forgotten-password links are valid for 24 hours and work only once. |
| Notes | Up to 10,000 characters. |
| Hydrogen atoms | Read and listed in the data table and the element statistics, but not bonded and not drawn in the graph. |
| Browser memory | Graph style and layout and the PDF option **Individual bonds** are remembered per browser. |

---

## Scientific background

All lengths are in Å, all angles in degrees.

### Reading the PDB file

- Bondamino reads the ATOM and HETATM records (coordinates, occupancy, B-factor, element, charge)
  and the records MODEL/ENDMDL, MODRES, SSBOND and CONECT. From the header it reads HEADER (PDB ID),
  EXPDTA (method), REMARK 2 (resolution) and REMARK 3 (R-work, R-free).
- A residue is identified by chain identifier, residue number and insertion code, so 52 and 52A are
  different residues. Residues are kept in file order.
- If the element column (columns 77–78) is empty, the element is taken from the first letter of the
  atom name.
- If two different residue names share the same chain, number and insertion code without alternate
  location identifiers (for example a ligand numbered like a protein residue), each name is treated
  as a residue of its own.
- **Standard amino acids** are the 20 proteinogenic ones. A HETATM residue that has the atoms N, CA
  and C is a **modified amino acid**; its parent amino acid comes from the MODRES record or, if there
  is none, from a built-in list (MSE, FME → Met; SEP → Ser; TPO → Thr; PTR → Tyr; CSO, CME, OCS,
  CSD → Cys; HYP → Pro; MLY, M3L, ALY, KCX, LLP → Lys; CGU → Glu; NEP → His). Every other residue,
  except water (HOH, DOD, WAT, H2O), is a **ligand**; a single-atom ligand is an ion.

### Models and the averaged model

- Files with several models (typically NMR ensembles) show the first model initially. Files without
  MODEL records are treated as one model.
- The pseudo-model **Average**: atoms and bonds are taken from the first model. Each bond length is
  the mean of the distances between the two atoms over all models in which both atoms exist; in a
  model where a peptide bond is a chain break (C–N > 2.0 Å) that model is left out for this bond.
  Distances are averaged unrounded and then rounded to 3 decimals. Coordinates, occupancy and
  B-factor of each atom are averaged over the models in which the atom occurs.
- In the averaged model, φ, ψ, ω and the bond angles are calculated from the **averaged
  coordinates**. They are therefore not the mean of the per-model angles. Likewise, the averaged
  bond lengths are generally not equal to the distances between averaged coordinates. Reference
  values and residue information are those of the first model.

### Alternate conformations

- Exactly one conformation is used per residue. With **Conformation: automatic**, it is the
  conformation with the highest mean occupancy of its atoms; on a tie, the first one in the file.
  Atoms without an alternate location identifier always belong to the residue.
- With a specific conformation selected (for example **Conformation B**), that conformation is used
  for every residue that has it; residues without it use the automatic choice.
- If an atom name occurs more than once in the selected set, the first occurrence is used.
- If the conformations of a position are different amino acids (microheterogeneity), the residue
  name of the chosen conformation is used.

### How bonds are determined

For standard amino acids, bonds are not read from the file; they follow from fixed topologies and
distance rules. Hydrogen atoms are ignored.

**Backbone** (every amino acid): N–CA, CA–C, C–O and, at the C-terminus, C–OXT.

**Side chains** (heavy atoms; CA–CB counts as side chain):

| Amino acid | Side-chain bonds |
|---|---|
| Ala | CA–CB |
| Arg | CA–CB, CB–CG, CG–CD, CD–NE, NE–CZ, CZ–NH1, CZ–NH2 |
| Asn | CA–CB, CB–CG, CG–OD1, CG–ND2 |
| Asp | CA–CB, CB–CG, CG–OD1, CG–OD2 |
| Cys | CA–CB, CB–SG |
| Gln | CA–CB, CB–CG, CG–CD, CD–OE1, CD–NE2 |
| Glu | CA–CB, CB–CG, CG–CD, CD–OE1, CD–OE2 |
| Gly | none |
| His | CA–CB, CB–CG, CG–ND1, CG–CD2, ND1–CE1, CD2–NE2, CE1–NE2 |
| Ile | CA–CB, CB–CG1, CB–CG2, CG1–CD1 |
| Leu | CA–CB, CB–CG, CG–CD1, CG–CD2 |
| Lys | CA–CB, CB–CG, CG–CD, CD–CE, CE–NZ |
| Met | CA–CB, CB–CG, CG–SD, SD–CE |
| Phe | CA–CB, CB–CG, CG–CD1, CG–CD2, CD1–CE1, CD2–CE2, CE1–CZ, CE2–CZ |
| Pro | CA–CB, CB–CG, CG–CD, CD–N (ring closure) |
| Ser | CA–CB, CB–OG |
| Thr | CA–CB, CB–OG1, CB–CG2 |
| Trp | CA–CB, CB–CG, CG–CD1, CG–CD2, CD1–NE1, CD2–CE2, CD2–CE3, NE1–CE2, CE2–CZ2, CE3–CZ3, CZ2–CH2, CZ3–CH2 |
| Tyr | CA–CB, CB–CG, CG–CD1, CG–CD2, CD1–CE1, CD2–CE2, CE1–CZ, CE2–CZ, CZ–OH |
| Val | CA–CB, CB–CG1, CB–CG2 |

A bond is drawn only if both atoms are present in the file.

**Peptide bonds and chain breaks**: between two consecutive residues of the same chain (in file
order), the bond C(i)–N(i+1) is formed if its length is at most **2.0 Å** (ideal value about
1.33 Å). A longer distance is a **chain break**, typically caused by residues missing from the model.

#### Disulfide bridges

Candidates are cysteines and cysteine-derived modified amino acids that have an SG atom. Two of them
are bonded (SG–SG):

1. if an SSBOND record names the pair, regardless of the distance (records that refer to a
   symmetry copy, symmetry operator other than 1555, are ignored, because those atoms are not in
   the file); shown as **SSBOND record**;
2. otherwise, if the SG–SG distance is at most **2.5 Å** (typical value about 2.0 Å); shown as
   **Distance**.

#### Ligands and modified amino acids

For residues that are neither standard amino acids nor water, bonds between their own heavy atoms
are taken from the **CONECT** records (only bonds within the residue; all residues with the same
name share the bonds found). If there are no CONECT records for that residue name, two atoms are
bonded if their distance $d$ satisfies

```math
d \le r_1 + r_2 + 0.45\ \text{Å}
```

where $r_1$ and $r_2$ are covalent radii (Cordero et al. 2008):

| Element | C | N | O | S | Se | P | F | Cl | Br | I | B |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Radius (Å) | 0.76 | 0.71 | 0.66 | 1.05 | 1.20 | 1.07 | 0.57 | 1.02 | 1.20 | 1.39 | 0.84 |

| Element | Fe | Zn | Mg | Cu | Mn | Co | Ni | Ca | Mo | W | V |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Radius (Å) | 1.32 | 1.22 | 1.41 | 1.32 | 1.39 | 1.26 | 1.24 | 1.76 | 1.54 | 1.62 | 1.53 |

Other elements use 0.80 Å. Bonds between different residues (for example to metal ions) are not
included. For a modified amino acid, the backbone bonds are the usual ones and all its other bonds
form its side chain.

### Bond angles

Every pair of bonded neighbours of an atom forms a bond angle with that atom at the centre:

```math
\theta = \arccos\frac{\vec{u}\cdot\vec{v}}{|\vec{u}|\,|\vec{v}|}
```

where $\vec{u}$ and $\vec{v}$ point from the centre atom to the two neighbours. This includes angles
across the peptide bond (CA–C–N, O–C–N, C–N–CA) and at disulfide bridges (CB–SG–SG); the latter have
no reference value.

### Reference values

Each bond and bond angle that has a reference value gets a mean $\mu$ and a standard deviation
$\sigma$.

**Backbone** (atoms N, CA, C, O, OXT only, and the peptide bond): Engh & Huber (2001), with the
assignment rules of the program Dangle (Richardson Lab). The **class** depends on the residue itself
(Gly, Pro), on the residue after the peptide bond (X–Gly, X–Pro, X–Pro (cis) for a cis peptide
bond before proline) or on the residue before it (Gly–X); **other** applies otherwise. Where several
classes apply, the first one listed wins.

| Bond | Values: class μ ± σ (Å) |
|---|---|
| N–CA | Gly 1.456 ± 0.015; Pro 1.468 ± 0.017; other 1.459 ± 0.020 |
| CA–C | Gly 1.514 ± 0.016; Pro 1.524 ± 0.020; other 1.525 ± 0.026 |
| C–O | Gly 1.232 ± 0.016; Pro 1.228 ± 0.020; other 1.229 ± 0.019 |
| C–N (peptide) | X–Gly 1.326 ± 0.018; X–Pro 1.338 ± 0.019; other 1.336 ± 0.023 |

| Angle | Values: class μ ± σ (°) |
|---|---|
| N–CA–C | Gly 113.1 ± 2.5; Pro 112.1 ± 2.6; other 111.0 ± 2.7 |
| CA–C–O | Gly 120.6 ± 1.8; Pro 120.2 ± 2.4; other 120.1 ± 2.1 |
| CA–C–N | Gly–X 116.2 ± 2.0; X–Pro 117.1 ± 2.8; other 117.2 ± 2.2 |
| O–C–N | X–Pro 121.1 ± 1.9; Gly–X 123.2 ± 1.7; other 122.7 ± 1.6 |
| C–N–CA | X–Pro (cis) 127.0 ± 2.4; X–Pro 119.3 ± 1.5; X–Gly 122.3 ± 2.1; other 121.7 ± 2.5 |

For the class X–Pro (cis), a peptide bond counts as cis if $|\omega| < 30°$. Modified amino acids use
the class of their parent amino acid. Peptide-bond references require both residues to be (possibly
modified) amino acids with known parent.

**Side chains**: the CCP4 Monomer Library. Used are all bond lengths and bond angles without hydrogen
that involve at least one atom outside the backbone (so also CA–CB and N–CA–CB). Values are looked
up by the residue's own name. Entries exist for the 20 standard amino acids and for the modified
amino acids ALY, CGU, CME, CSD, CSO, FME, HYP, KCX, LLP, M3L, MLY, MSE, NEP, OCS, PTR, SEP and TPO.

**No reference** exists for ligand bonds, disulfide bridges, angles involving them, and side chains
of other modified amino acids.

### Z-score, RMSZ, warnings and outliers

For a bond length or bond angle $x$ with reference $\mu \pm \sigma$:

```math
Z = \frac{x - \mu}{\sigma}
```

A positive Z means longer (or wider) than the reference. For a set of $n$ assessed values:

```math
\mathrm{RMSZ} = \sqrt{\frac{1}{n}\sum_{i=1}^{n} Z_i^2}
```

An RMSZ near 1 means the deviations are as large as the reference standard deviations predict.

| \|Z\| | Rating | Colour |
|---|---|---|
| ≤ 2 | within expectation | green |
| > 2 and ≤ 4 | warning | orange |
| > 4 | outlier | red |
| no reference | not assessed | grey |

Z-scores of bond lengths are calculated from the lengths as displayed (rounded to 0.001 Å); those
of bond angles from the unrounded angles.

### Dihedral angles

The dihedral (torsion) angle of four atoms $p_0$–$p_1$–$p_2$–$p_3$ is the angle between the plane
of $p_0, p_1, p_2$ and the plane of $p_1, p_2, p_3$. With
$\vec{b}_1 = p_1 - p_0$, $\vec{b}_2 = p_2 - p_1$, $\vec{b}_3 = p_3 - p_2$:

```math
\theta = \operatorname{atan2}\!\left(|\vec{b}_2|\;\vec{b}_1\cdot(\vec{b}_2\times\vec{b}_3),\;(\vec{b}_1\times\vec{b}_2)\cdot(\vec{b}_2\times\vec{b}_3)\right)
```

The range is −180° to 180°. The sign follows the IUPAC-IUB convention: looking along the central
bond from $p_1$ to $p_2$, the angle is positive if the front bond must be turned clockwise to cover
the rear bond.

| Angle | Atoms |
|---|---|
| φ (phi) | C(i−1)–N(i)–CA(i)–C(i) |
| ψ (psi) | N(i)–CA(i)–C(i)–N(i+1) |
| ω (omega) | CA(i)–C(i)–N(i+1)–CA(i+1) |

φ and ψ are calculated only for residues that are joined by peptide bonds to both neighbours, so
not at chain ends, next to chain breaks, or when atoms are missing. ω is calculated for every
peptide bond. In the **Dihedral angles** CSV, ω is listed with the residue after the peptide bond.

### Peptide bond conformation

| \|ω\| | Classification |
|---|---|
| < 30° | **cis** |
| 30° to < 150° | **twisted** |
| ≥ 150° | **trans** |

cis peptide bonds before a proline are common; elsewhere they are rare and are highlighted. The same
30° limit decides the cis class of the backbone reference values. For the Ramachandran category of
proline, a wider limit applies (see below).

### Ramachandran evaluation

Bondamino scores φ/ψ as MolProbity does, using the Top8000 reference distributions of the Richardson
Lab (Williams et al. 2018).

**Category** of each amino acid (the first that applies; modified amino acids count as their parent):

1. **Glycine**: Gly.
2. **cis-proline** or **trans-proline**: Pro; cis if the preceding peptide bond has $|\omega| < 90°$,
   otherwise trans.
3. **Pre-proline**: the next residue is a proline.
4. **Ile and Val**: Ile or Val.
5. **General**: all others.

**Value**: each category has a grid of 180 × 180 cells of 2° × 2° over φ and ψ. The value at (φ, ψ)
is interpolated linearly between the four nearest cell centres; the grid wraps around at ±180°.

**Region** from the value $v$:

| Region | Condition |
|---|---|
| **Favoured** | $v \ge 0.02$ (2 %) |
| **Allowed** | below favoured, and $v \ge$ the category limit: 0.0005 (0.05 %) for General; 0.002 (0.2 %) for cis-proline; 0.001 (0.1 %) for Glycine, trans-proline, Pre-proline, Ile and Val |
| **Outliers** | below the category limit |

The shaded regions in the plot show the same favoured and allowed areas of the selected category.
The point colours of the plot (other, glycine, proline) are a simpler grouping than the six
categories.

### Statistics over all models

For every bond, Bondamino collects its length in every model in which it is formed (with the
selected conformation). From these unrounded lengths it computes the mean and the median, rounded to
0.001 Å; these appear as **Average (Å)** and **Median (Å)** in **Shortest and longest bond**. It also
finds the single shortest and longest bond over all models and the model in which each occurs
(**All models**). Files without MODEL records count as one model. The mean over all models equals the
bond length of the **Average** model.

### Colour scales

| Mode | Scale |
|---|---|
| **By length** | Hue from 0° (red, shortest displayed bond) to 120° (green, longest), linear in bond length, in HSL colour space. Disulfide bridges are not used to set the range. If all bonds have the same length, they are yellow. |
| **Reference** | By \|Z\| class, see the table above. |
| **B-factor** | Hue from 230° (blue, low) via violet to 360° (red, high), linear in B. The ends are the 2nd and 98th percentiles of the B-factors of the drawn atoms, so single extreme values do not compress the scale; values beyond the ends get the end colours. A bond gets the colour of the mean B-factor of its two atoms. Atoms or bonds without B-factor are grey. |

### References

Brief list; full references are in [CITATION.md](CITATION.md).

- Engh, R. A. & Huber, R. (2001). Structure quality and target parameters. *International Tables for
  Crystallography*, Vol. F, ch. 18.3, 382–392. Backbone reference values.
- Dangle, Richardson Lab, Duke University. Assignment rules of the backbone reference classes.
- CCP4 Monomer Library. Side-chain reference values.
- Williams, C. J. et al. (2018). MolProbity: More and better reference data for improved all-atom
  structure validation. *Protein Science* 27, 293–315. Top8000 Ramachandran distributions.
- Cordero, B. et al. (2008). Covalent radii revisited. *Dalton Transactions*, 2832–2838. Covalent radii.
- IUPAC-IUB Commission on Biochemical Nomenclature (1970). Abbreviations and symbols for the
  description of the conformation of polypeptide chains. Sign convention of dihedral angles.

---

## Glossary

| Term | Meaning |
|---|---|
| Alternate conformation, altLoc | Two or more positions of the same atoms in a crystal structure, marked by a letter (A, B, …) in column 17 of the PDB file |
| Atom serial number | The running number of an atom in the PDB file (**No.**, **Atom no.**) |
| B-factor | Temperature factor (Å²): describes how much an atom's position is smeared out by motion or disorder; higher means less well defined |
| Backbone | The repeating chain N–CA–C(=O) of a protein |
| Chain | One polypeptide (or other molecule) in a structure, identified by a letter |
| Chain break | Two consecutive residues of a chain that are not joined by a peptide bond, usually because residues are missing from the model |
| CONECT record | Line in a PDB file that lists bonds of an atom; used for ligands |
| CSV | A plain-text table format that spreadsheet programs open directly |
| Dihedral angle, torsion angle | The angle between two planes defined by four consecutive atoms |
| Disulfide bridge | A covalent S–S bond between two cysteines |
| HETATM | PDB record for atoms that are not part of a standard residue: ligands, ions, water, modified amino acids |
| Insertion code | A letter after the residue number (for example 52A) that distinguishes inserted residues |
| JSON | A plain-text format for structured data, readable by programs |
| Ligand | A small molecule bound to the protein, for example ATP or haem |
| Microheterogeneity | Different amino acids at the same position in alternate conformations |
| mmCIF (PDBx/mmCIF) | The current standard file format of the Protein Data Bank; not supported by Bondamino |
| Model | One complete set of coordinates; NMR files usually contain many models |
| Modified amino acid | A chemically altered amino acid, for example selenomethionine (MSE) or phosphoserine (SEP) |
| MODRES record | Line in a PDB file that names the parent amino acid of a modified residue |
| Occupancy | Fraction (0 to 1) of molecules in the crystal in which an atom is at the given position |
| PDB | Protein Data Bank; also its classic fixed-column file format |
| PDB ID | The four-character code of an entry in the Protein Data Bank, for example 1LYZ |
| φ, ψ, ω (phi, psi, omega) | Backbone dihedral angles, see [dihedral angles](#dihedral-angles) |
| Ramachandran plot | Plot of ψ against φ for every amino acid; shows which backbone conformations occur |
| RCSB | Research Collaboratory for Structural Bioinformatics, which runs the US Protein Data Bank site from which Bondamino downloads structures |
| Reference value | Expected mean and standard deviation of a bond length or angle, from high-quality structures or small-molecule data |
| Residue | One amino acid (or ligand, ion, water molecule) in a structure |
| RMSZ | Root mean square of Z-scores, see [Z-score](#z-score-rmsz-warnings-and-outliers) |
| R-work, R-free | Agreement between a crystal structure model and the measured data; R-free uses data left out of refinement |
| SSBOND record | Line in a PDB file that lists a disulfide bridge |
| Working group | A set of Bondamino users who share files and notes |
| Z-score | Deviation from a reference value in units of its standard deviation |
