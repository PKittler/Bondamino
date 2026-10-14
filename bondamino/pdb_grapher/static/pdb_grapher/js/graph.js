/**
 * @fileoverview Front-end of the Bondamino detail page (template
 * `pdb_grapher/templates/pdb_grapher/graph.html`).
 *
 * The page shows one uploaded PDB structure in five Bootstrap tabs; this file
 * drives the first four (the "notes" tab is handled by `notes.js`):
 *
 *   - `graph`        Interactive d3 graph: one horizontal band per residue,
 *                    backbone N–CA–C vertically on the left, side chains as a
 *                    tree (rings as polygons) to the right, ligands laid out from
 *                    their 3D coordinates. Bonds are coloured by length, uniformly,
 *                    by Z-score against reference values or by B-factor. Two
 *                    layouts: `schematic` (fixed grid) and `proportional` (edge
 *                    length proportional to the bond length in Å).
 *   - `ramachandran` Ramachandran plot (φ/ψ), ω (cis/trans) and bond angles,
 *                    each with reference statistics.
 *   - `data`         Atom table (one row per ATOM/HETATM record).
 *   - `stats`        Composition, bond statistics (reference Z-scores, peptide
 *                    bonds, disulfides, extremes over all models, length
 *                    histogram, bonds grouped by atom pair) and exports.
 *
 * A sequence bar (one-letter codes) above the graph and data tabs is linked to
 * both views: clicking a cell jumps to the residue, scrolling highlights the
 * residues currently in view.
 *
 * Data flow:
 *   - `#graph-data` (JSON, built by `utils.build_graph`): `{residues, nodes,
 *     edges}` of the initially selected model. Atom ids have the form
 *     `<chain>_<pos>_<atom>`; edges carry `type` (backbone, sidechain, peptide,
 *     disulfide, ligand), `length` (Å) and `ref` ([mean, sigma, group] or null);
 *     residues may carry `rama` (φ/ψ and Top8000 evaluation from the server).
 *   - `#bond-ensemble` (JSON, `utils.build_bond_ensemble`): per-bond mean/median
 *     over all models plus the shortest and longest bond of all models.
 *   - `#structure-info` (JSON): header metadata (method, resolution, R values)
 *     and the names of the reference sources; read lazily via `structureInfo()`.
 *   - `GET /pdb/<pk>/models/?model=&altloc=&ensemble=1` (URL in
 *     `#structure-selects[data-model-url]`): returns `{graph, ensemble?}` when the
 *     user picks another model or alternate conformation; everything is then
 *     re-rendered from the new `graphData` (see `loadStructure` in the init block).
 *   - Ramachandran regions are a static JSON file fetched on demand
 *     (`#angles-export[data-rama-regions]`, see `loadRamaRegions`).
 *
 * Exports (stats and Ramachandran tabs): JSON and PDF are both generated from
 * `buildStatsExport()`; the PDF uses jsPDF + jspdf-autotable, loaded lazily on
 * the first PDF export (`loadScript`). CSV tables come from `csvTable()`, the
 * histogram and Ramachandran plot can be saved as SVG/PNG/JPEG.
 *
 * Key state:
 *   - `canvas` (from `setupCanvas`): d3 selections, the zoom behaviour and
 *     `canvas.state` (active tab, selected residue, band bounds, pending jumps).
 *   - `graphData`, `bondEnsemble`, `activeStyle`, `activeLayout` live in the
 *     closure of the init IIFE at the end of the file; style and layout are
 *     persisted in localStorage.
 *   - Module-level: `ramaView`, `histogram` and the `expanded*` sets that keep
 *     collapsible table rows open across re-renders (e.g. after a model change).
 *
 * All user-visible strings go through Django's JavaScript catalogue (`gettext`,
 * `ngettext`, `interpolate`); the German source strings are the message ids and
 * must not be changed here without updating the translations.
 *
 * File layout (sections marked with `// ===== ... =====`): constants, helpers,
 * ring geometry, graph layout, rendering, B-factor colouring, sequence bar, data
 * table, statistics, legend, reference values, angles/Ramachandran, bond
 * statistics, histogram, JSON/PDF export, CSV export, initialisation.
 */
'use strict';

// ===== Constants =====
/** @const {number} Height of one residue band in the schematic layout (px). */
const BAND_HEIGHT = 200;
/** @const {number} Width reserved for the residue labels on the left (px). */
const LABEL_W = 90;
/** @const {number} x position of the backbone column (px). */
const BACKBONE_X = LABEL_W + 80;
/** @const {number} Horizontal step per side-chain tree level, schematic layout (px). */
const SC_COL_W = 82;
/** @const {number} Radius of an atom circle (px). */
const NODE_R = 16;
/** @const {number} Vertical space per leaf of the side-chain tree, schematic layout (px). */
const LEAF_H = 46;
const PROP_SCALE = 55; // px per Å in proportional layout
const SC_RING_SIDE = 60; // Edge length of rings (and ligand bonds) in the schematic layout (px)
/** @const {number} Length of the separator line between two residue bands (px). */
const RESIDUE_SEPARATOR_W = 100;

/** @const {!Object<string, number>} y offsets of backbone atoms within a band (schematic layout). */
const BB_Y = { N: 30, CA: 90, C: 150, O: 160 };
/** @const {!Set<string>} Atom names that belong to the backbone. */
const BACKBONE_ATOMS = new Set(['N', 'CA', 'C', 'O', 'OXT']);

/**
 * Ring-shaped side chains: atom names in ring order (the last one closes back to
 * the first). TRP has two fused rings sharing the edge CD2–CE2. The PRO ring
 * closes back into the backbone atoms N/CA.
 * @const {!Object<string, !Array<!Array<string>>>}
 */
const SIDECHAIN_RINGS = {
    HIS: [['CG', 'ND1', 'CE1', 'NE2', 'CD2']],
    PHE: [['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']],
    TYR: [['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']],
    TRP: [['CG', 'CD1', 'NE1', 'CE2', 'CD2'], ['CD2', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3']],
    PRO: [['N', 'CA', 'CB', 'CG', 'CD']],
};

// ===== Helpers =====
/**
 * Page language (`<html lang>`) used for number and date formats; texts come from
 * Django's JavaScript catalogue (gettext, ngettext, interpolate).
 * @const {string}
 */
const LOCALE = document.documentElement.lang || 'de';

/**
 * Reads a CSS custom property from the document root (follows the light/dark theme).
 * @param {string} name Property name including the leading `--`.
 * @return {string} The trimmed value, or '' if undefined.
 */
function cssVar(name) {
    return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

/**
 * Colour of a chemical element from the CSS variable `--atom-<ELEMENT>`; elements
 * without their own colour fall back to `--atom-default`.
 * @param {string} element Element symbol (any case).
 * @return {string} A CSS colour value.
 */
function getElementColor(element) {
    return cssVar(`--atom-${(element || '').toUpperCase()}`) || cssVar('--atom-default');
}

// ===== Ring geometry =====
// Side-chain rings are drawn with their real bond lengths in the proportional
// layout and as regular polygons in the schematic layout (`bondPx` then returns a
// constant edge length).

/**
 * Computes the vertices of a closed polygon with the given edge lengths whose
 * vertices all lie on a common circle (cyclic polygon). This way every edge,
 * including the ring-closing one, gets exactly its real proportional length
 * without stretching any edge to close the ring. The circumradius is found by
 * bisection so that the central angles sum up to 2π.
 * @param {!Array<number>} sides Edge lengths in ring order (px).
 * @return {!Array<{x: number, y: number}>} Vertices around the origin.
 */
function layoutCyclicPolygon(sides) {
    const n = sides.length;
    const maxSide = Math.max(...sides);
    const angleSum = R => sides.reduce((s, len) => s + 2 * Math.asin(Math.min(1, len / (2 * R))), 0);

    let lo = maxSide / 2 + 1e-6;
    let hi = lo * 2;
    while (angleSum(hi) > 2 * Math.PI) hi *= 2;
    for (let i = 0; i < 60; i++) {
        const mid = (lo + hi) / 2;
        if (angleSum(mid) > 2 * Math.PI) lo = mid; else hi = mid;
    }
    const R = (lo + hi) / 2;

    const pts = [];
    let cum = 0;
    for (let i = 0; i < n; i++) {
        pts.push({ x: R * Math.cos(cum), y: R * Math.sin(cum) });
        cum += 2 * Math.asin(Math.min(1, sides[i] / (2 * R)));
    }
    return pts;
}

/**
 * Rotates a point around the origin.
 * @param {{x: number, y: number}} p The point.
 * @param {number} cosA Cosine of the rotation angle.
 * @param {number} sinA Sine of the rotation angle.
 * @return {{x: number, y: number}} The rotated point.
 */
function rotatePoint(p, cosA, sinA) {
    return { x: p.x * cosA - p.y * sinA, y: p.x * sinA + p.y * cosA };
}

/**
 * Attaches a single ring vertex to an already fixed world position and orients the
 * ring so that it "opens" from there in direction `outDir` (away from the backbone).
 * @param {!Array<{x: number, y: number}>} localPts Ring vertices from layoutCyclicPolygon.
 * @param {number} vertexIdx Index of the vertex to pin.
 * @param {{x: number, y: number}} worldPos Target position of that vertex.
 * @param {{x: number, y: number}} outDir Direction from the vertex to the ring centre.
 * @return {!Array<{x: number, y: number}>} Ring vertices in world coordinates.
 */
function anchorRingSingle(localPts, vertexIdx, worldPos, outDir) {
    const p = localPts[vertexIdx];
    const localAngleToCenter = Math.atan2(-p.y, -p.x);
    const worldAngle = Math.atan2(outDir.y, outDir.x);
    const rotation = worldAngle - localAngleToCenter;
    const cosA = Math.cos(rotation), sinA = Math.sin(rotation);
    const pRot = rotatePoint(p, cosA, sinA);
    const tx = worldPos.x - pRot.x, ty = worldPos.y - pRot.y;
    return localPts.map(pt => {
        const r = rotatePoint(pt, cosA, sinA);
        return { x: r.x + tx, y: r.y + ty };
    });
}

/**
 * Attaches two adjacent ring vertices to two already fixed world positions (e.g.
 * backbone N/CA for PRO, or the shared edge of two fused rings for TRP). Since two
 * points fix the ring only up to a reflection, `chooseMirror` picks the side.
 * @param {!Array<{x: number, y: number}>} localPts Ring vertices from layoutCyclicPolygon.
 * @param {number} idxA Index of the first pinned vertex.
 * @param {number} idxB Index of the second pinned vertex.
 * @param {{x: number, y: number}} worldA Target position of vertex `idxA`.
 * @param {{x: number, y: number}} worldB Target position of vertex `idxB`.
 * @param {function(!Array, !Array): !Array} chooseMirror Receives the unmirrored
 *     and the mirrored placement and returns the one to use.
 * @return {!Array<{x: number, y: number}>} Ring vertices in world coordinates.
 */
function anchorRingDouble(localPts, idxA, idxB, worldA, worldB, chooseMirror) {
    function build(mirror) {
        const pts = mirror ? localPts.map(p => ({ x: p.x, y: -p.y })) : localPts;
        const a = pts[idxA], b = pts[idxB];
        const localAngle = Math.atan2(b.y - a.y, b.x - a.x);
        const worldAngle = Math.atan2(worldB.y - worldA.y, worldB.x - worldA.x);
        const rotation = worldAngle - localAngle;
        const cosA = Math.cos(rotation), sinA = Math.sin(rotation);
        const aRot = rotatePoint(a, cosA, sinA);
        const tx = worldA.x - aRot.x, ty = worldA.y - aRot.y;
        return pts.map(pt => {
            const r = rotatePoint(pt, cosA, sinA);
            return { x: r.x + tx, y: r.y + ty };
        });
    }
    return chooseMirror(build(false), build(true));
}

/**
 * Places all rings of one side chain, in order, so that a second fused ring can use
 * the already placed atoms of the first one as anchors (see TRP). A ring with two
 * already positioned atoms is anchored on that edge; otherwise its first atom is
 * placed to the right of its (tree-placed) neighbour outside the ring.
 * @param {!Array<!Array<string>>} ringDefs Rings as atom names (see SIDECHAIN_RINGS).
 * @param {!Object<string, {x: number, y: number}>} positions Atom id -> position;
 *     updated in place.
 * @param {!Object<string, !Array<string>>} scAdj Side-chain adjacency by atom id.
 * @param {function(string, string): number} bondPx Edge length between two atom ids (px).
 * @param {function(string): string} pfx Maps an atom name to the residue's atom id.
 */
function layoutSidechainRings(ringDefs, positions, scAdj, bondPx, pfx) {
    const priorRingPositions = []; // Atoms of previous rings of this side chain (e.g. TRP ring 1)

    ringDefs.forEach(ringNames => {
        const ids = ringNames.map(pfx);
        const n = ids.length;
        const sides = ids.map((id, i) => bondPx(id, ids[(i + 1) % n]));
        const localPts = layoutCyclicPolygon(sides);

        const fixedIdx = [];
        ids.forEach((id, i) => { if (positions[id]) fixedIdx.push(i); });
        const others = priorRingPositions.filter(p => !ids.includes(p.id));

        let placed;
        if (fixedIdx.length >= 2) {
            const idxA = fixedIdx[0], idxB = fixedIdx[1];
            const worldA = positions[ids[idxA]], worldB = positions[ids[idxB]];
            placed = anchorRingDouble(localPts, idxA, idxB, worldA, worldB, (optA, optB) => {
                // For a fused ring (e.g. TRP) pick the side pointing away from the
                // already placed first ring instead of overlapping it; otherwise
                // (e.g. PRO) simply the side with larger x (further from the backbone).
                const score = pts => others.length > 0
                    ? pts.reduce((s, p) => s + others.reduce((s2, o) => s2 + Math.hypot(p.x - o.pos.x, p.y - o.pos.y), 0), 0)
                    : pts.reduce((s, p) => s + p.x, 0);
                return score(optA) >= score(optB) ? optA : optB;
            });
        } else {
            const anchorAtomId = (scAdj[ids[0]] || []).find(nid => positions[nid] && !ids.includes(nid));
            if (!anchorAtomId) return;
            const anchorPos = positions[anchorAtomId];
            const vertex0Pos = {
                x: anchorPos.x + bondPx(anchorAtomId, ids[0]),
                y: anchorPos.y,
            };
            placed = anchorRingSingle(localPts, 0, vertex0Pos, { x: 1, y: 0 });
        }

        ids.forEach((id, i) => { if (!positions[id]) positions[id] = placed[i]; });
        ids.forEach(id => priorRingPositions.push({ id, pos: positions[id] }));
    });
}

/**
 * Rotates the ring system of a side chain around its first atom (e.g. CG) so that
 * its centroid lies at the height of that atom. Otherwise fused rings (TRP) bend
 * diagonally downwards and reach far into the next residue's band.
 * @param {!Array<!Array<string>>} ringDefs Rings as atom names.
 * @param {!Object<string, {x: number, y: number}>} positions Updated in place.
 * @param {function(string): string} pfx Maps an atom name to the residue's atom id.
 */
function balanceRingSystem(ringDefs, positions, pfx) {
    const ids = [...new Set(ringDefs.flat().map(pfx))].filter(id => positions[id]);
    const pivotId = pfx(ringDefs[0][0]);
    const pivot = positions[pivotId];
    if (!pivot || ids.length < 3) return;
    const cx = d3.mean(ids, id => positions[id].x);
    const cy = d3.mean(ids, id => positions[id].y);
    const angle = -Math.atan2(cy - pivot.y, cx - pivot.x);
    const cosA = Math.cos(angle), sinA = Math.sin(angle);
    ids.forEach(id => {
        if (id === pivotId) return;
        const r = rotatePoint({ x: positions[id].x - pivot.x, y: positions[id].y - pivot.y }, cosA, sinA);
        positions[id] = { x: pivot.x + r.x, y: pivot.y + r.y };
    });
}

/**
 * Places atoms branching off a ring (e.g. OH of TYR) radially away from the
 * (approximate) ring centre, at distance `bondPx`.
 * @param {!Array<!Array<string>>} ringDefs Rings as atom names.
 * @param {!Object<string, {x: number, y: number}>} positions Updated in place.
 * @param {!Object<string, !Object>} atomById Graph nodes by atom id.
 * @param {!Object<string, !Array<string>>} scAdj Side-chain adjacency by atom id.
 * @param {function(string, string): number} bondPx Edge length between two atom ids (px).
 * @param {function(string): string} pfx Maps an atom name to the residue's atom id.
 * @param {string} resPrefix Atom id prefix of the residue (`<chain>_<pos>_`).
 */
function placeRingBranches(ringDefs, positions, atomById, scAdj, bondPx, pfx, resPrefix) {
    ringDefs.forEach(ringNames => {
        const ids = ringNames.map(pfx);
        const centroid = ids.reduce((c, id) => ({ x: c.x + positions[id].x / ids.length, y: c.y + positions[id].y / ids.length }), { x: 0, y: 0 });

        ids.forEach(ringId => {
            (scAdj[ringId] || []).forEach(nid => {
                if (positions[nid] || !nid.startsWith(resPrefix)) return;
                const a = atomById[nid];
                if (!a || BACKBONE_ATOMS.has(a.label) || ids.includes(nid)) return;

                const ringPos = positions[ringId];
                let dirX = ringPos.x - centroid.x, dirY = ringPos.y - centroid.y;
                const len = Math.hypot(dirX, dirY) || 1;
                dirX /= len; dirY /= len;
                const dist = bondPx(ringId, nid);
                positions[nid] = { x: ringPos.x + dirX * dist, y: ringPos.y + dirY * dist };
            });
        });
    });
}

const LIGAND_PAD = 28; // Padding between ligand and the edge of its band (px)

/**
 * 2D layout of a ligand derived from its 3D coordinates: projection onto the two
 * principal axes (the longer one horizontal), then a short force simulation
 * untangles overlapping atoms.
 * @param {!Array<!Object>} atoms Graph nodes of the ligand (with x/y/z in Å).
 * @param {!Array<!Object>} bonds Ligand edges (`source`, `target`, `length`).
 * @param {function({length: number}): number} linkPx Target length of a bond in px.
 * @return {{positions: !Object<string, {x: number, y: number}>, width: number,
 *     height: number}} Positions relative to the top-left corner, and the extent.
 */
function layoutLigand(atoms, bonds, linkPx) {
    const n = atoms.length;
    const mean = ['x', 'y', 'z'].map(k => d3.mean(atoms, a => a[k]));
    const centered = atoms.map(a => [a.x - mean[0], a.y - mean[1], a.z - mean[2]]);
    const cov = [0, 1, 2].map(i => [0, 1, 2].map(j => d3.sum(centered, p => p[i] * p[j])));
    const mul = v => cov.map(row => row[0] * v[0] + row[1] * v[1] + row[2] * v[2]);
    const norm = v => { const l = Math.hypot(...v) || 1; return v.map(c => c / l); };
    // Principal axes of the covariance matrix via power iteration (deterministic,
    // no random start); the second axis is kept orthogonal to the first.
    const axis = orthogonalTo => {
        let v = orthogonalTo ? norm([orthogonalTo[1], -orthogonalTo[0], 0.3]) : [1, 0.6, 0.3];
        for (let i = 0; i < 50; i++) {
            v = mul(v);
            if (orthogonalTo) {
                const d = v[0] * orthogonalTo[0] + v[1] * orthogonalTo[1] + v[2] * orthogonalTo[2];
                v = v.map((c, k) => c - d * orthogonalTo[k]);
            }
            v = norm(v);
        }
        return v;
    };
    const e1 = axis(null), e2 = axis(e1);
    const meanBond = d3.mean(bonds, b => b.length) || 1.5;
    const nodes = atoms.map((a, i) => {
        const p = centered[i];
        const scale = linkPx({ length: meanBond }) / meanBond;
        return { id: a.id, x: (p[0] * e1[0] + p[1] * e1[1] + p[2] * e1[2]) * scale,
                 y: (p[0] * e2[0] + p[1] * e2[1] + p[2] * e2[2]) * scale };
    });
    if (n > 1) {
        d3.forceSimulation(nodes)
            .force('link', d3.forceLink(bonds.map(b => ({ ...b }))).id(d => d.id).distance(linkPx).strength(1))
            .force('repel', d3.forceManyBody().strength(-40).distanceMax(linkPx({ length: meanBond }) * 2.5))
            .force('collide', d3.forceCollide(NODE_R + 6))
            .stop()
            .tick(300);
    }
    const [x0, x1] = d3.extent(nodes, d => d.x), [y0, y1] = d3.extent(nodes, d => d.y);
    return {
        positions: Object.fromEntries(nodes.map(d => [d.id, { x: d.x - x0, y: d.y - y0 }])),
        width: x1 - x0, height: y1 - y0,
    };
}

// ===== Layout =====
/**
 * Computes the 2D position of every atom and the vertical extent of every residue
 * band. Residues are stacked top to bottom in file order. Per amino acid: backbone
 * N, CA, C in a vertical column at BACKBONE_X, O/OXT to the left, the side chain as
 * a tree growing to the right from CA, rings placed as polygons. Ligands (residues
 * with ligand edges) get their own band laid out by layoutLigand.
 *
 * In the `schematic` layout bands have a fixed height (BAND_HEIGHT) and fixed
 * spacings; in the `proportional` layout all distances are bond lengths ×
 * PROP_SCALE and band heights follow from the peptide bonds.
 * @param {!Object} graphData Graph data (`residues`, `nodes`, `edges`).
 * @param {string} layout `'schematic'` or `'proportional'`.
 * @return {{positions: !Object<string, {x: number, y: number}>,
 *     bandBounds: !Array<{yStart: number, yEnd: number}>}} Atom positions by id
 *     and band bounds indexed by residue index.
 */
function computePositions(graphData, layout) {
    const isProp = layout === 'proportional';
    const positions = {};
    const atomById = {};
    graphData.nodes.forEach(n => { atomById[n.id] = n; });

    // Sidechain adjacency
    const scAdj = {};
    graphData.nodes.forEach(n => { scAdj[n.id] = []; });
    graphData.edges.filter(e => e.type === 'sidechain').forEach(e => {
        scAdj[e.source].push(e.target);
        scAdj[e.target].push(e.source);
    });

    // Bond length lookup (all edge types)
    const bl = {};
    graphData.edges.forEach(e => {
        if (!bl[e.source]) bl[e.source] = {};
        if (!bl[e.target]) bl[e.target] = {};
        bl[e.source][e.target] = e.length;
        bl[e.target][e.source] = e.length;
    });
    // Falls back to 1.5 Å (typical C–C distance) for unbonded or missing atom pairs.
    function getBL(a, b) { return (bl[a] && bl[a][b]) ? bl[a][b] : 1.5; }

    const leafH = isProp ? 1.5 * PROP_SCALE : LEAF_H;
    // y of the next N atom in the proportional layout.
    let yPropCursor = 20;
    const bandBounds = [];

    // Bonds within ligands, grouped by residue (atom id without the atom name).
    const ligandBonds = d3.group(graphData.edges.filter(e => e.type === 'ligand'),
        e => e.source.slice(0, e.source.lastIndexOf('_')));

    graphData.residues.forEach(res => {
        const ri = res.index;
        const pfx = s => `${res.chain}_${res.pos}_${s}`;

        // Bands are stacked without gaps; ligand bands are as tall as needed.
        const yStart = ri === 0 ? 0 : bandBounds[ri - 1].yEnd;

        // --- Ligand: own band; in the schematic layout at least as tall as an amino acid's ---
        const ligBonds = ligandBonds.get(`${res.chain}_${res.pos}`);
        if (ligBonds) {
            const ids = [...new Set(ligBonds.flatMap(e => [e.source, e.target]))];
            const linkPx = isProp ? b => b.length * PROP_SCALE : () => SC_RING_SIDE;
            const lig = layoutLigand(ids.map(id => atomById[id]), ligBonds, linkPx);
            const bandH = isProp ? lig.height + 2 * LIGAND_PAD : Math.max(BAND_HEIGHT, lig.height + 2 * LIGAND_PAD);
            const top = yStart + (bandH - lig.height) / 2;
            Object.entries(lig.positions).forEach(([id, p]) => {
                positions[id] = { x: BACKBONE_X + p.x, y: top + p.y };
            });
            const yEnd = yStart + bandH;
            bandBounds.push({ yStart, yEnd });
            if (isProp) yPropCursor = yEnd + PROP_SCALE * 0.66;
            return;
        }

        const nId = pfx('N'), caId = pfx('CA'), cId = pfx('C');
        const oId = pfx('O'), oxtId = pfx('OXT');
        const bandY = yStart;

        // --- Backbone Y positions ---
        let yN, yCA, yC;
        if (isProp) {
            yN  = yPropCursor;
            yCA = yN  + (atomById[nId]  && atomById[caId] ? getBL(nId,  caId) : 1.46) * PROP_SCALE;
            yC  = yCA + (atomById[caId] && atomById[cId]  ? getBL(caId, cId)  : 1.52) * PROP_SCALE;
        } else {
            yN  = bandY + BB_Y.N;
            yCA = bandY + BB_Y.CA;
            yC  = bandY + BB_Y.C;
        }

        if (atomById[nId])   positions[nId]  = { x: BACKBONE_X, y: yN  };
        if (atomById[caId])  positions[caId] = { x: BACKBONE_X, y: yCA };
        if (atomById[cId])   positions[cId]  = { x: BACKBONE_X, y: yC  };

        if (atomById[oId]) {
            const dx = isProp ? getBL(cId, oId) * PROP_SCALE : 58;
            positions[oId] = { x: BACKBONE_X - dx, y: isProp ? yC : bandY + BB_Y.O };
        }
        if (atomById[oxtId]) {
            // Offset from O (otherwise the two carboxylate oxygens almost coincide).
            // In the proportional layout via an angle instead of a pixel offset so
            // that the edge length still equals the real bond length.
            if (isProp) {
                const dist = getBL(cId, oxtId) * PROP_SCALE;
                const angle = Math.PI * (35 / 180);
                positions[oxtId] = {
                    x: BACKBONE_X - dist * Math.cos(angle),
                    y: yC + dist * Math.sin(angle),
                };
            } else {
                positions[oxtId] = { x: BACKBONE_X - 58 - 40, y: bandY + BB_Y.C + 30 };
            }
        }

        // --- Sidechain tree ---
        const caPos = positions[caId];
        if (caPos) {
            const resPrefix = `${res.chain}_${res.pos}_`;
            const cbList = (scAdj[caId] || []).filter(id => {
                const a = atomById[id];
                return a && id.startsWith(resPrefix) && !BACKBONE_ATOMS.has(a.label);
            });

            if (cbList.length > 0) {
                const scVisited = new Set([...BACKBONE_ATOMS].map(a => pfx(a)));

                // Ring atoms are placed separately with ring geometry and therefore
                // excluded from the regular tree traversal; drawn as a tree, the
                // ring-closing bonds would cross each other (e.g. for TRP).
                // In the schematic layout PRO stays a tree: its ring runs through the
                // backbone atoms N/CA, which are fixed to the grid.
                // Modified amino acids (e.g. PTR, HYP) inherit the rings of their parent
                // amino acid as long as the atom names match.
                const ringDefs = (SIDECHAIN_RINGS[res.name] || SIDECHAIN_RINGS[res.parent] || [])
                    .filter(ring => ring.every(name => atomById[pfx(name)]))
                    .filter(ring => isProp || !ring.some(name => BACKBONE_ATOMS.has(name)));
                const bondPx = isProp ? (a, b) => getBL(a, b) * PROP_SCALE : () => SC_RING_SIDE;
                ringDefs.forEach(ringNames => ringNames.forEach(name => scVisited.add(pfx(name))));

                // Depth-first spanning tree of the side chain within this residue.
                function buildTree(id, parentId) {
                    scVisited.add(id);
                    const children = (scAdj[id] || []).filter(nid => {
                        if (scVisited.has(nid)) return false;
                        if (!nid.startsWith(resPrefix)) return false;
                        const a = atomById[nid];
                        return a && !BACKBONE_ATOMS.has(a.label);
                    });
                    return { id, parentId, children: children.map(c => buildTree(c, id)) };
                }

                function countLeaves(node) {
                    return node.children.length === 0 ? 1
                        : node.children.reduce((s, c) => s + countLeaves(c), 0);
                }

                // Each tree level moves one step right; every subtree gets vertical
                // space proportional to its number of leaves and is centred within it.
                function assignPos(node, xParent, yTop) {
                    const dx = isProp ? getBL(node.parentId, node.id) * PROP_SCALE : SC_COL_W;
                    const x = xParent + dx;
                    const totalH = countLeaves(node) * leafH;
                    const y = yTop + totalH / 2;
                    positions[node.id] = { x, y };
                    let cy = yTop;
                    node.children.forEach(child => {
                        const ch = countLeaves(child) * leafH;
                        assignPos(child, x, cy);
                        cy += ch;
                    });
                }

                const trees = cbList
                    .filter(cbId => !scVisited.has(cbId))
                    .map(cbId => buildTree(cbId, caId));
                const totalLeaves = trees.reduce((s, t) => s + countLeaves(t), 0);
                let yTop = caPos.y - (totalLeaves * leafH) / 2;
                trees.forEach(tree => {
                    const lH = countLeaves(tree) * leafH;
                    assignPos(tree, caPos.x, yTop);
                    yTop += lH;
                });

                // Place ring atoms only after their anchor (e.g. CB) got a position above.
                if (ringDefs.length > 0) {
                    layoutSidechainRings(ringDefs, positions, scAdj, bondPx, pfx);
                    // Do not rotate the PRO ring: it hangs on the fixed backbone atoms N/CA.
                    const anchoredToBackbone = ringDefs.some(ring => ring.some(name => BACKBONE_ATOMS.has(name)));
                    if (!anchoredToBackbone) balanceRingSystem(ringDefs, positions, pfx);
                    placeRingBranches(ringDefs, positions, atomById, scAdj, bondPx, pfx, resPrefix);
                }
            }
        }

        // --- Band bounds & cursor advance ---
        // Proportional layout: the next N follows at the distance of the peptide
        // bond C–N; the band boundary lies in the middle of that bond.
        if (isProp) {
            const nextRes = graphData.residues[ri + 1];
            const peptideAdv = nextRes
                ? getBL(cId, `${nextRes.chain}_${nextRes.pos}_N`) * PROP_SCALE
                : 30;
            const yBandStart = ri === 0 ? 0 : bandBounds[ri - 1].yEnd;
            const yBandEnd   = yC + peptideAdv / 2;
            bandBounds.push({ yStart: yBandStart, yEnd: yBandEnd });
            yPropCursor = yC + peptideAdv;
        } else {
            bandBounds.push({ yStart: bandY, yEnd: bandY + BAND_HEIGHT });
        }
    });

    return { positions, bandBounds };
}

// ===== Render (D3) =====
/** @const {number} Duration of layout/style transitions (ms). */
const TRANSITION_MS = 600;
/** @const {number} Extra pan space around the content (px). */
const VIEW_PAD = 40;

/**
 * Stable key of a residue (d3 data joins, selection state, sequence bar).
 * @param {!Object} r Residue from `graphData.residues`.
 * @return {string} `<chain>_<pos>_<name>`.
 */
const residueKey = r => `${r.chain}_${r.pos}_${r.name}`;
/**
 * Stable key of an edge; also matches the keys of `bondEnsemble.bonds`.
 * @param {!Object} e Edge from `graphData.edges`.
 * @return {string} `<source>|<target>`.
 */
const edgeKey = e => `${e.source}|${e.target}`;

const d3DefaultConstrain = d3.zoom().constrain();
/**
 * Zoom constraint: d3's default constraint, plus the view may never start right of
 * the left content edge (translate-x ≤ 0), not even when zooming out, where d3
 * would otherwise centre the content.
 * @param {!d3.ZoomTransform} transform Proposed transform.
 * @param {!Array<!Array<number>>} extent Viewport extent.
 * @param {!Array<!Array<number>>} translateExtent Allowed content extent.
 * @return {!d3.ZoomTransform} The constrained transform.
 */
function constrainLeftEdge(transform, extent, translateExtent) {
    const t = d3DefaultConstrain(transform, extent, translateExtent);
    return t.x > 0 ? d3.zoomIdentity.translate(0, t.y).scale(t.k) : t;
}

/**
 * One-time setup of the graph SVG: a viewport group (transformed by d3.zoom) with
 * the layers bands, labels, edges and nodes (drawn in this order). The mouse wheel
 * pans vertically (Shift + wheel horizontally), Ctrl + wheel (or trackpad pinch)
 * zooms, dragging pans the view.
 * @return {!Object} The canvas object shared by all render functions: `svg`,
 *     `viewport`, `zoom` and `state`. Render functions later attach further
 *     selections to it (`seqCells`, `tableRows`, `tableScroll`,
 *     `firstRowByResidue`, `ramaPoints`, `ramaAngles`).
 */
function setupCanvas() {
    const svg = d3.select('#protein-svg');
    const viewport = svg.append('g').attr('class', 'viewport');
    ['bands', 'labels', 'edges', 'nodes'].forEach(name => {
        viewport.append('g').attr('class', `layer-${name}`);
    });

    let canvas;
    const zoom = d3.zoom()
        .scaleExtent([0.1, 4])
        .extent(() => [[0, 0], [svg.node().clientWidth, svg.node().clientHeight]])
        // The zoom behaviour only handles Ctrl + wheel; plain wheel events are
        // handled by the 'wheel.pan' listener below.
        .filter(event => (event.type === 'wheel' ? event.ctrlKey : !event.ctrlKey) && !event.button)
        // d3 multiplies the delta by ten while Ctrl is held (for pinch gestures); use
        // a uniform factor instead so that one mouse-wheel notch does not zoom too far.
        .wheelDelta(event => -event.deltaY * (event.deltaMode === 1 ? 0.05 : event.deltaMode ? 1 : 0.005))
        .constrain(constrainLeftEdge)
        .on('zoom', event => {
            viewport.attr('transform', event.transform);
            if (canvas && canvas.state.activeTab === 'graph') updateSequenceView(canvas, event.transform);
        });

    svg.call(zoom)
        // Plain wheel pans the view; deltaMode 1/2 = lines/pages, converted to px.
        .on('wheel.pan', event => {
            if (event.ctrlKey) return;
            event.preventDefault();
            const unit = event.deltaMode === 1 ? 16 : event.deltaMode === 2 ? svg.node().clientHeight : 1;
            let dx = event.deltaX * unit;
            let dy = event.deltaY * unit;
            if (event.shiftKey && !dx) { dx = dy; dy = 0; }
            const k = d3.zoomTransform(svg.node()).k;
            svg.call(zoom.translateBy, -dx / k, -dy / k);
        }, { passive: false });

    canvas = {
        svg, viewport, zoom,
        // residues/bandBounds: last rendered layout; selectedKey/selectedResidue: residue
        // selected via sequence bar or Ramachandran plot; activeTab: data-tab of the
        // visible tab; pendingJump: views that still have to jump to the selection
        // the next time they are shown.
        state: { residues: [], bandBounds: [], selectedKey: null, activeTab: 'graph', pendingJump: new Set() },
    };
    return canvas;
}

/**
 * Draws (or updates) the graph tab with d3 data joins, so that switching layout or
 * bond style animates atoms and bonds to their new positions/colours. Also re-renders
 * the sequence bar and the B-factor legend.
 * @param {!Object} canvas Canvas object from setupCanvas.
 * @param {!Object} graphData Graph data (`residues`, `nodes`, `edges`).
 * @param {string} bondStyle `'gradient'` (length, red → green), `'mono'`,
 *     `'reference'` (Z-score class) or `'bfactor'`.
 * @param {string} layout `'schematic'` or `'proportional'`.
 * @param {boolean} animate Whether to use transitions.
 */
function render(canvas, graphData, bondStyle, layout, animate) {
    const { positions, bandBounds } = computePositions(graphData, layout);

    const allPos = Object.values(positions);
    if (allPos.length === 0) return;

    const { svg, viewport, zoom } = canvas;
    const viewW = svg.node().clientWidth;
    const viewH = svg.node().clientHeight;
    // Gutter for disulfide bridges right of all atoms, one lane per bridge
    const disulfideCount = graphData.edges.filter(e => e.type === 'disulfide').length;
    const gutterX = d3.max(allPos, p => p.x) + NODE_R + SS_GUTTER_GAP;
    const maxX = d3.max(allPos, p => p.x) + NODE_R + 60 + (disulfideCount ? SS_GUTTER_GAP + disulfideCount * SS_LANE_W : 0);
    const maxY = d3.max(allPos, p => p.y) + NODE_R + 20;
    const contentW = Math.max(maxX, 500, viewW);
    const contentH = layout === 'proportional'
        ? maxY + 20
        : (bandBounds.length ? bandBounds[bandBounds.length - 1].yEnd : 0) + 10;

    canvas.state.residues = graphData.residues;
    canvas.state.bandBounds = bandBounds;
    renderSequence(canvas, graphData);

    zoom.translateExtent([[-VIEW_PAD, -VIEW_PAD], [contentW + VIEW_PAD, Math.max(contentH, viewH) + VIEW_PAD]]);
    // Clamp the current transform to the new bounds (e.g. after a layout switch).
    svg.call(zoom.translateBy, 0, 0);

    const t = d3.transition().duration(animate ? TRANSITION_MS : 0);
    const fadeIn  = sel => sel.attr('opacity', 0).call(s => s.transition(t).attr('opacity', 1));
    const fadeOut = sel => sel.transition(t).attr('opacity', 0).remove();

    const gBands  = viewport.select('.layer-bands');
    const gLabels = viewport.select('.layer-labels');
    const gEdges  = viewport.select('.layer-edges');
    const gNodes  = viewport.select('.layer-nodes');

    // --- Bands & labels ---
    const bandResidues = graphData.residues.filter(r => bandBounds[r.index]);
    const bandY = r => bandBounds[r.index].yStart;

    // Short separator line on the left at every boundary between two residues
    // (instead of coloured bands, which side chains of neighbouring residues may
    // reach into).
    gBands.selectAll('line.residue-separator')
        .data(bandResidues.filter(r => r.index > 0), residueKey)
        .join(
            enter => enter.append('line')
                .attr('class', 'residue-separator')
                .attr('x1', 0).attr('x2', RESIDUE_SEPARATOR_W)
                .attr('y1', bandY).attr('y2', bandY)
                .call(fadeIn),
            update => update,
            fadeOut,
        )
        .transition(t)
        .attr('y1', bandY)
        .attr('y2', bandY);

    const pos = (r, atom) => positions[`${r.chain}_${r.pos}_${atom}`];
    // Faint vertical guide along the backbone from N to C of each amino acid.
    gBands.selectAll('line.backbone-guide')
        .data(bandResidues.filter(r => pos(r, 'N') && pos(r, 'C')), residueKey)
        .join(
            enter => enter.append('line')
                .attr('class', 'backbone-guide')
                .attr('x1', BACKBONE_X).attr('x2', BACKBONE_X)
                .attr('y1', r => pos(r, 'N').y - 5).attr('y2', r => pos(r, 'C').y + 5)
                .attr('stroke-width', 0.5)
                .call(fadeIn),
            update => update,
            fadeOut,
        )
        .attr('stroke', cssVar('--bs-border-color'))
        .attr('stroke-opacity', 0.3)
        .transition(t)
        .attr('y1', r => pos(r, 'N').y - 5)
        .attr('y2', r => pos(r, 'C').y + 5);

    // Residue label (one-letter code, three-letter code, position), vertically centred in the band.
    const labelTransform = r => `translate(8,${(bandBounds[r.index].yStart + bandBounds[r.index].yEnd) / 2})`;
    gLabels.selectAll('g.residue-label')
        .data(bandResidues, residueKey)
        .join(
            enter => {
                const g = enter.append('g').attr('class', 'residue-label').attr('transform', labelTransform);
                g.append('text').attr('class', 'residue-label-text label-one').attr('y', -10).attr('font-weight', 700);
                g.append('text').attr('class', 'residue-label-text label-three').attr('y', 6);
                g.append('text').attr('class', 'residue-label-text label-seq').attr('y', 20)
                    .attr('font-size', '10px').attr('fill-opacity', 0.7);
                return g.call(fadeIn);
            },
            update => update,
            fadeOut,
        )
        .classed('selected', r => residueKey(r) === canvas.state.selectedKey)
        .call(g => g.select('.label-one').text(r => r.one_letter))
        .call(g => g.select('.label-three').text(r => r.name))
        .call(g => g.select('.label-seq').text(r => `#${r.pos}`))
        .transition(t)
        .attr('transform', labelTransform);

    // --- Bond color ---
    const visibleEdges = graphData.edges.filter(e => positions[e.source] && positions[e.target]);
    // Disulfide bridges connect distant residues and are drawn as arcs; the colour
    // scale is based on the other bonds only (SG–SG is considerably longer).
    const straightEdges = visibleEdges.filter(e => e.type !== 'disulfide');
    const arcEdges = visibleEdges.filter(e => e.type === 'disulfide');
    const [minLen, maxLen] = d3.extent(straightEdges, e => e.length);
    const hue = d3.scaleLinear().domain([minLen, maxLen]).range([0, 120]).clamp(true);
    const monoColor = cssVar('--bond-mono');

    const zColors = Object.fromEntries(['ok', 'warn', 'outlier', 'none'].map(c => [c, cssVar(`--z-${c}`)]));
    const nodeById = new Map(graphData.nodes.map(n => [n.id, n]));
    const bScale = bFactorScale(graphData.nodes.filter(n => positions[n.id]));
    renderBFactorLegend(bScale);
    const bColor = b => b === null || b === undefined || !bScale ? zColors.none : bScale.color(b);
    // Bond colour for the selected style; 'gradient' maps the length range of the
    // shown bonds to hue 0 (red, shortest) … 120 (green, longest). In 'bfactor'
    // mode a bond gets the mean B-factor of its two atoms.
    function bondColor(e) {
        if (bondStyle === 'mono') return monoColor;
        if (bondStyle === 'reference') return zColors[zClass(zScore(e))];
        if (bondStyle === 'bfactor') {
            const bs = [nodeById.get(e.source), nodeById.get(e.target)].map(n => n?.temp_factor);
            return bColor(bs.some(b => b === null || b === undefined) ? null : d3.mean(bs));
        }
        return d3.hsl(minLen === maxLen ? 60 : hue(e.length), 0.8, 0.42).toString();
    }
    const atomColor = n => bondStyle === 'bfactor' ? bColor(n.temp_factor) : getElementColor(n.element);

    // --- Edges ---
    // Straight bonds; side-chain bonds dashed, ligand bonds dash-dotted. The
    // tooltip shows length, bond type and reference value.
    const atomLabel = id => id.slice(id.lastIndexOf('_') + 1);
    gEdges.selectAll('line.bond')
        .data(straightEdges, edgeKey)
        .join(
            enter => enter.append('line')
                .attr('class', 'bond')
                .attr('stroke-linecap', 'round')
                .attr('x1', e => positions[e.source].x).attr('y1', e => positions[e.source].y)
                .attr('x2', e => positions[e.target].x).attr('y2', e => positions[e.target].y)
                .attr('stroke', e => bondColor(e))
                .call(sel => sel.append('title'))
                .call(fadeIn),
            update => update,
            fadeOut,
        )
        .attr('stroke-width', e => e.type === 'sidechain' || e.type === 'ligand' ? 1.8 : 2.5)
        .attr('stroke-dasharray', e => e.type === 'sidechain' ? '5,3' : e.type === 'ligand' ? '8,3,2,3' : null)
        .call(sel => sel.select('title').text(e =>
            `${atomLabel(e.source)}–${atomLabel(e.target)}: ${e.length} Å (${BOND_TYPE_LABELS[e.type] || e.type})\n${referenceText(e)}`))
        .transition(t)
        .attr('stroke', e => bondColor(e))
        .attr('x1', e => positions[e.source].x).attr('y1', e => positions[e.source].y)
        .attr('x2', e => positions[e.target].x).attr('y2', e => positions[e.target].y);

    // Disulfide bridge as a bracket: from SG horizontally right into its own lane of
    // the gutter, there vertically to the partner and horizontally back (rounded
    // corners). This way it crosses no side chains of other residues. `mid` is the
    // middle of the vertical segment.
    const lanes = new Map(arcEdges.map((e, i) => [edgeKey(e), i]));
    function arcGeometry(e) {
        const src = positions[e.source], tgt = positions[e.target];
        const x = gutterX + lanes.get(edgeKey(e)) * SS_LANE_W;
        const dir = tgt.y >= src.y ? 1 : -1;
        const r = Math.min(10, Math.abs(tgt.y - src.y) / 2, (x - Math.max(src.x, tgt.x)) / 2);
        return {
            d: `M${src.x},${src.y} H${x - r} Q${x},${src.y} ${x},${src.y + dir * r} ` +
               `V${tgt.y - dir * r} Q${x},${tgt.y} ${x - r},${tgt.y} H${tgt.x}`,
            mid: { x, y: (src.y + tgt.y) / 2 },
        };
    }
    gEdges.selectAll('path.disulfide-bond')
        .data(arcEdges, edgeKey)
        .join(
            enter => enter.append('path')
                .attr('class', 'disulfide-bond')
                .attr('d', e => arcGeometry(e).d)
                .attr('stroke', e => bondColor(e))
                .call(sel => sel.append('title'))
                .call(fadeIn),
            update => update,
            fadeOut,
        )
        .call(sel => sel.select('title').text(e =>
            `${atomLabel(e.source)}–${atomLabel(e.target)}: ${e.length} Å (${BOND_TYPE_LABELS.disulfide})\n${referenceText(e)}`))
        .transition(t)
        .attr('stroke', e => bondColor(e))
        .attr('d', e => arcGeometry(e).d);
    d3.select('#legend-disulfide').property('hidden', !graphData.edges.some(e => e.type === 'disulfide'));
    d3.select('#legend-ligand').property('hidden', !graphData.edges.some(e => e.type === 'ligand'));

    // Length label next to the edge midpoint (offset perpendicular to the edge);
    // hidden for edges shorter than 20 px.
    function bondLabel(e) {
        if (e.type === 'disulfide') {
            const { mid } = arcGeometry(e);
            return { x: mid.x + 14, y: mid.y, visible: true };
        }
        const src = positions[e.source], tgt = positions[e.target];
        const dx = tgt.x - src.x, dy = tgt.y - src.y;
        const pxLen = Math.sqrt(dx * dx + dy * dy) || 1;
        return {
            x: (src.x + tgt.x) / 2 - dy / pxLen * 9,
            y: (src.y + tgt.y) / 2 + dx / pxLen * 9,
            visible: pxLen > 20,
        };
    }

    gEdges.selectAll('text.bond-length-text')
        .data(visibleEdges, edgeKey)
        .join(
            enter => enter.append('text')
                .attr('class', 'bond-length-text')
                .attr('x', e => bondLabel(e).x).attr('y', e => bondLabel(e).y)
                .attr('opacity', 0),
            update => update,
            fadeOut,
        )
        .text(e => e.length)
        .transition(t)
        .attr('x', e => bondLabel(e).x)
        .attr('y', e => bondLabel(e).y)
        .attr('opacity', e => bondLabel(e).visible ? 1 : 0);

    // --- Nodes ---
    // Bond angles with this atom at the centre (for the tooltip)
    const anglesByCenter = analyzeAngles(graphData).byCenter;
    const nodeTransform = n => `translate(${positions[n.id].x},${positions[n.id].y})`;
    gNodes.selectAll('g.atom-group')
        .data(graphData.nodes.filter(n => positions[n.id]), n => n.id)
        .join(
            enter => {
                const g = enter.append('g').attr('class', 'atom-group').attr('transform', nodeTransform);
                g.append('title');
                g.append('circle').attr('class', 'atom-circle').attr('r', NODE_R)
                    .attr('fill', atomColor);
                g.append('text').attr('class', 'atom-text');
                return g.call(fadeIn);
            },
            update => update,
            fadeOut,
        )
        .call(g => g.select('title').text(n => [
            `${n.label} · ${n.residue} ${n.residue_seq}${n.i_code || ''} (${interpolate(gettext('Kette %s'), [n.chain])})`,
            ...(n.temp_factor === null || n.temp_factor === undefined ? []
                : [interpolate(gettext('B-Faktor: %s Å²'), [formatBFactor(n.temp_factor)])]),
            ...(anglesByCenter.get(n.id) || []).map(a =>
                `${a.type}: ${formatAngle(a.value)}${a.z === null ? '' : ` · Z = ${formatZ(a.z)}`}`),
        ].join('\n')))
        .call(g => g.select('text').text(n => n.label))
        .call(g => g.select('circle').transition(t).attr('fill', atomColor))
        .transition(t)
        .attr('transform', nodeTransform);
}

// ===== Colouring by B-factor =====
// Blue (low) via violet to red (high). The scale spans the 2nd to the 98th
// percentile of the displayed atoms so that single extreme values do not compress
// the colours.
const B_HUE_LOW = 230;
const B_HUE_HIGH = 360;
const formatBFactor = v => v.toLocaleString(LOCALE, { minimumFractionDigits: 2, maximumFractionDigits: 2 });

/**
 * Builds the B-factor colour scale for the given atoms.
 * @param {!Array<!Object>} nodes Graph nodes (`temp_factor` may be null).
 * @return {?{low: number, high: number, clippedLow: boolean, clippedHigh: boolean,
 *     color: function(number): string}} The scale, or null if no atom has a B-factor.
 */
function bFactorScale(nodes) {
    const values = nodes.map(n => n.temp_factor).filter(b => b !== null && b !== undefined).sort(d3.ascending);
    if (!values.length) return null;
    const low = d3.quantileSorted(values, 0.02);
    const high = d3.quantileSorted(values, 0.98);
    const hue = d3.scaleLinear().domain([low, high]).range([B_HUE_LOW, B_HUE_HIGH]).clamp(true);
    return {
        low, high,
        clippedLow: values[0] < low, clippedHigh: values[values.length - 1] > high,
        color: b => d3.hsl(low === high ? (B_HUE_LOW + B_HUE_HIGH) / 2 : hue(b), 0.75, 0.45).toString(),
    };
}

/**
 * Updates the B-factor legend (`#legend-bfactor`): gradient bar and end values,
 * prefixed with ≤/≥ when values outside the percentile range were clipped.
 * @param {?Object} scale Result of bFactorScale.
 */
function renderBFactorLegend(scale) {
    const legend = d3.select('#legend-bfactor');
    if (legend.empty()) return;
    legend.select('.bfactor-empty').property('hidden', Boolean(scale));
    legend.select('.bfactor-scale').property('hidden', !scale);
    if (!scale) return;
    const stops = d3.range(0, 1.0001, 0.25).map(f => scale.color(scale.low + f * (scale.high - scale.low)));
    legend.select('.bfactor-bar').style('background', `linear-gradient(to right, ${stops.join(', ')})`);
    legend.select('.bfactor-low').text(`${scale.clippedLow ? '≤ ' : ''}${formatBFactor(scale.low)}`);
    legend.select('.bfactor-high').text(`${scale.clippedHigh ? '≥ ' : ''}${formatBFactor(scale.high)}`);
}

// ===== Sequence bar =====
// Disulfide gutter in the graph (used by render): gap to the atoms and width per lane
const SS_GUTTER_GAP = 40;
const SS_LANE_W = 16;

// Sequence bar geometry (px): cell size, height of the number row above the cells,
// extra gap between chains and at chain breaks / cis peptide bonds.
const SEQ_CELL_W = 16;
const SEQ_CELL_H = 22;
const SEQ_TICK_H = 14;
const SEQ_CHAIN_GAP = 10;
const SEQ_BREAK_GAP = 4;

/**
 * Analyses peptide bonds between consecutive amino acids of the same chain.
 * `possible` counts all such pairs, `breaks` the chain breaks among them, i.e.
 * pairs without a peptide bond because the server rejected the C–N distance as
 * too large (e.g. residues missing from the model). Ligands (without N/CA/C) do
 * not count.
 * @param {!Object} graphData Graph data.
 * @return {{possible: number, formed: number,
 *     breaks: !Map<number, {prev: !Object, distance: number}>}} `breaks` is keyed
 *     by the index of the residue after the gap and holds the residue before it
 *     and the C–N distance in Å.
 */
function analyzePeptideBonds(graphData) {
    const nodes = new Map(graphData.nodes.map(n => [n.id, n]));
    const peptides = new Set(graphData.edges.filter(e => e.type === 'peptide').map(edgeKey));
    const atom = (res, name) => nodes.get(`${res.chain}_${res.pos}_${name}`);
    const isAminoAcid = res => ['N', 'CA', 'C'].every(name => atom(res, name));
    const breaks = new Map();
    let possible = 0;
    graphData.residues.forEach((res, i) => {
        const prev = graphData.residues[i - 1];
        if (!prev || prev.chain !== res.chain || !isAminoAcid(prev) || !isAminoAcid(res)) return;
        possible++;
        const c = atom(prev, 'C'), n = atom(res, 'N');
        if (peptides.has(`${c.id}|${n.id}`)) return;
        breaks.set(res.index, { prev, distance: Math.hypot(c.x - n.x, c.y - n.y, c.z - n.z) });
    });
    return { possible, formed: possible - breaks.size, breaks };
}

// Width of one character of the number labels (9px Roboto Mono) plus minimum gap
const SEQ_TICK_CHAR_W = 5.5;
const SEQ_TICK_GAP = 6;

/**
 * Numbers above the sequence bar: chain start ("A 1") and every multiple of ten.
 * If a label would overlap the previous one, the multiple-of-ten label is dropped;
 * chain starts are always kept.
 * @param {!Array<{res: !Object, x: number, chainStart: boolean}>} cells Sequence cells.
 * @return {!Array<{cell: !Object, label: (string|number), chainStart: boolean}>} Ticks.
 */
function sequenceTicks(cells) {
    const end = t => t.cell.x + t.label.length * SEQ_TICK_CHAR_W + SEQ_TICK_GAP;
    const ticks = [];
    cells.forEach(cell => {
        const res = cell.res;
        if (!cell.chainStart && (res.seq % 10 !== 0 || res.i_code)) return;
        const tick = { cell, label: cell.chainStart ? `${res.chain} ${res.pos}` : res.pos, chainStart: cell.chainStart };
        const prev = ticks[ticks.length - 1];
        if (prev && end(prev) > cell.x) {
            if (!tick.chainStart) return;
            if (!prev.chainStart) ticks.pop();
        }
        ticks.push(tick);
    });
    return ticks;
}

/**
 * Microheterogeneity: the alternate conformations (altLoc) of a residue are
 * different amino acids.
 * @param {!Object} res Residue with `conformations` (altLoc -> residue name).
 * @return {?string} E.g. "PRO (A) · SER (B, C)", or null if all conformations
 *     are the same amino acid.
 */
function residueVariants(res) {
    const byName = d3.groups(Object.entries(res.conformations || {}), ([, name]) => name);
    if (byName.length < 2) return null;
    return byName.map(([name, items]) => `${name} (${items.map(([altloc]) => altloc).join(', ')})`).join(' · ');
}

/**
 * Renders the one-letter sequence as a row of cells (`#sequence-svg`). Chain
 * changes get a gap, chain breaks within a chain a red line, cis peptide bonds
 * their own marker. Cells are marked for modified residues, disulfide cysteines
 * and microheterogeneity. Clicking a cell selects the residue and jumps to it in
 * the visible view (see jumpToResidue). Stores the cell selection in
 * `canvas.seqCells`.
 * @param {!Object} canvas Canvas object.
 * @param {!Object} graphData Graph data.
 */
function renderSequence(canvas, graphData) {
    const svg = d3.select('#sequence-svg');
    if (svg.empty()) return;

    const residues = graphData.residues;
    const { breaks } = analyzePeptideBonds(graphData);
    const bridges = disulfidePartners(graphData);
    // cis peptide bonds, keyed by the index of the residue after the bond
    const cisBonds = new Map(analyzeAngles(graphData).omegas.filter(o => o.kind === 'cis').map(o => [o.res.index, o]));
    const cisMarks = [];
    let x = 0;
    const gaps = [];
    const cells = residues.map((res, i) => {
        const chainStart = i === 0 || res.chain !== residues[i - 1].chain;
        if (chainStart && i > 0) x += SEQ_CHAIN_GAP;
        const chainBreak = breaks.get(res.index);
        if (chainBreak) {
            x += SEQ_BREAK_GAP;
            // Middle of the space between the two cells (cells are SEQ_CELL_W - 2 wide)
            gaps.push({ ...chainBreak, res, x: x - (SEQ_BREAK_GAP + 2) / 2 });
        }
        const cis = cisBonds.get(res.index);
        if (cis) {
            x += SEQ_BREAK_GAP;
            cisMarks.push({ ...cis, x: x - (SEQ_BREAK_GAP + 2) / 2 });
        }
        const cell = { res, x, chainStart };
        x += SEQ_CELL_W;
        return cell;
    });
    const cellKey = c => residueKey(c.res);

    // Extra width so that a tick label at the last cell is not cut off.
    svg.attr('width', x + 40).attr('height', SEQ_TICK_H + SEQ_CELL_H);

    svg.selectAll('text.seq-tick')
        .data(sequenceTicks(cells), t => cellKey(t.cell))
        .join('text')
        .attr('class', 'seq-tick')
        .attr('x', t => t.cell.x + 1)
        .attr('y', SEQ_TICK_H - 4)
        .text(t => t.label);

    canvas.seqCells = svg.selectAll('g.seq-cell')
        .data(cells, cellKey)
        .join(enter => {
            const g = enter.append('g');
            g.append('title');
            g.append('rect').attr('width', SEQ_CELL_W - 2).attr('height', SEQ_CELL_H).attr('rx', 3);
            // Dot top right: modified amino acid; top left: cysteine with disulfide bridge
            g.append('circle').attr('class', 'seq-mod-mark').attr('cx', SEQ_CELL_W - 5).attr('cy', 3.5).attr('r', 2);
            g.append('circle').attr('class', 'seq-ss-mark').attr('cx', 3).attr('cy', 3.5).attr('r', 2);
            // Mark at the bottom edge: conformations are different amino acids
            g.append('rect').attr('class', 'seq-conf-mark')
                .attr('x', 2).attr('y', SEQ_CELL_H - 3).attr('width', SEQ_CELL_W - 6).attr('height', 2.5).attr('rx', 1.25);
            g.append('text').attr('x', (SEQ_CELL_W - 2) / 2).attr('y', SEQ_CELL_H / 2);
            return g;
        })
        .attr('class', c => `seq-cell${c.res.one_letter === '?' ? ' het' : ''}${c.res.parent ? ' modified' : ''}${residueVariants(c.res) ? ' variant' : ''}${bridges.has(residueKey(c.res)) ? ' disulfide' : ''}`)
        .attr('transform', c => `translate(${c.x},${SEQ_TICK_H})`)
        .on('click', (event, c) => jumpToResidue(canvas, c.res))
        .call(g => g.select('title').text(c => {
            const lines = [`${c.res.name} ${c.res.pos} (${interpolate(gettext('Kette %s'), [c.res.chain])})`];
            if (c.res.parent) lines.push(interpolate(gettext('Modifizierte Aminosäure, abgeleitet von %s'), [aminoAcidLabel(c.res.parent)]));
            const variants = residueVariants(c.res);
            if (variants) lines.push(interpolate(gettext('Aminosäure je Konformation: %s'), [variants]));
            (bridges.get(residueKey(c.res)) || []).forEach(partner =>
                lines.push(interpolate(gettext('Disulfidbrücke zu %s'), [partner])));
            return lines.join('\n');
        }))
        .call(g => g.select('text').text(c => c.res.one_letter));

    svg.selectAll('g.seq-break')
        .data(gaps, g => residueKey(g.res))
        .join(enter => {
            const g = enter.append('g').attr('class', 'seq-break');
            g.append('title');
            // Wider invisible hit area so that the tooltip is easy to trigger
            g.append('rect').attr('class', 'seq-break-hit').attr('x', -4).attr('width', 8)
                .attr('y', 0).attr('height', SEQ_TICK_H + SEQ_CELL_H);
            g.append('line').attr('y1', 1).attr('y2', SEQ_TICK_H + SEQ_CELL_H);
            return g;
        })
        .attr('transform', g => `translate(${g.x},0)`)
        .call(g => g.select('title').text(g => interpolate(
            gettext('Kettenlücke: keine Peptidbindung zwischen %(from)s und %(to)s (Abstand %(distance)s Å)'),
            {
                from: `${g.prev.name} ${g.prev.pos}`,
                to: `${g.res.name} ${g.res.pos}`,
                distance: g.distance.toLocaleString(LOCALE, { maximumFractionDigits: 1 }),
            }, true)));

    svg.selectAll('g.seq-cis')
        .data(cisMarks, m => residueKey(m.res))
        .join(enter => {
            const g = enter.append('g').attr('class', 'seq-cis');
            g.append('title');
            g.append('rect').attr('class', 'seq-break-hit').attr('x', -4).attr('width', 8)
                .attr('y', 0).attr('height', SEQ_TICK_H + SEQ_CELL_H);
            g.append('line').attr('y1', 1).attr('y2', SEQ_TICK_H + SEQ_CELL_H);
            return g;
        })
        .attr('transform', m => `translate(${m.x},0)`)
        .call(g => g.select('title').text(m => interpolate(
            gettext('cis-Peptidbindung zwischen %(from)s und %(to)s (ω = %(omega)s)'),
            { from: `${m.prev.name} ${m.prev.pos}`, to: `${m.res.name} ${m.res.pos}`, omega: formatAngle(m.omega) }, true)));

    updateSelection(canvas);
}

/**
 * Highlights the sequence cells whose band currently lies in the visible part of
 * the graph (called on every zoom/pan while the graph tab is active).
 * @param {!Object} canvas Canvas object.
 * @param {!d3.ZoomTransform} transform Current zoom transform.
 */
function updateSequenceView(canvas, transform) {
    if (!canvas || !canvas.seqCells) return;
    const { bandBounds } = canvas.state;
    const viewH = canvas.svg.node().clientHeight;
    const yTop = -transform.y / transform.k;
    const yBottom = (viewH - transform.y) / transform.k;
    markSequenceInView(canvas, c => {
        const b = bandBounds[c.res.index];
        return b && b.yEnd > yTop && b.yStart < yBottom;
    });
}

/**
 * Highlights the cells for which `inView` is true and scrolls the sequence bar
 * along if that range leaves its visible part (centred if it fits, otherwise
 * aligned to the first cell).
 * @param {!Object} canvas Canvas object.
 * @param {function(!Object): boolean} inView Predicate on a sequence cell.
 */
function markSequenceInView(canvas, inView) {
    if (!canvas.seqCells) return;
    canvas.seqCells.classed('in-view', inView);

    const scroller = document.getElementById('sequence-scroll');
    const visible = canvas.seqCells.filter(inView).data();
    if (!scroller || visible.length === 0) return;
    const first = visible[0].x;
    const last = visible[visible.length - 1].x + SEQ_CELL_W;
    const left = scroller.scrollLeft, width = scroller.clientWidth;
    if (first >= left && last <= left + width) return;
    scroller.scrollLeft = last - first > width
        ? first - 8
        : (first + last) / 2 - width / 2;
}

/**
 * Applies the `selected` class for the selected residue in all views: graph
 * labels, sequence cells, data table rows and Ramachandran points (the selected
 * point is raised to the top).
 * @param {!Object} canvas Canvas object.
 */
function updateSelection(canvas) {
    const key = canvas.state.selectedKey;
    canvas.viewport.selectAll('g.residue-label').classed('selected', r => residueKey(r) === key);
    if (canvas.seqCells) canvas.seqCells.classed('selected', c => residueKey(c.res) === key);
    if (canvas.tableRows) canvas.tableRows.classed('selected', n => nodeResidueKey(n) === key);
    if (canvas.ramaPoints) canvas.ramaPoints.classed('selected', p => residueKey(p.res) === key).filter('.selected').raise();
}

/**
 * Selects a residue and jumps to it in the visible view (graph or data table);
 * the other of the two views catches up on the jump as soon as its tab is shown
 * (via `state.pendingJump`). When called from another tab (Ramachandran plot),
 * both views jump later.
 * @param {!Object} canvas Canvas object.
 * @param {!Object} res Residue from `graphData.residues`.
 */
function jumpToResidue(canvas, res) {
    canvas.state.selectedKey = residueKey(res);
    canvas.state.selectedResidue = res;
    updateSelection(canvas);

    const { activeTab, pendingJump } = canvas.state;
    pendingJump.clear();
    ['graph', 'data'].forEach(view => {
        if (view !== activeTab) pendingJump.add(view);
    });
    if (activeTab === 'graph') jumpGraphToResidue(canvas, res, 500);
    if (activeTab === 'data') scrollTableToResidue(canvas, res, true);
}

/**
 * Scrolls the graph vertically so that the residue's band is centred.
 * @param {!Object} canvas Canvas object.
 * @param {!Object} res Residue.
 * @param {number} duration Transition duration (ms); 0 jumps immediately.
 */
function jumpGraphToResidue(canvas, res, duration) {
    const b = canvas.state.bandBounds[res.index];
    if (!b) return;
    // Keep horizontal position and zoom level, only jump vertically to the band.
    const node = canvas.svg.node();
    const t = d3.zoomTransform(node);
    const mid = (b.yStart + b.yEnd) / 2;
    const dy = (node.clientHeight / 2 - (mid * t.k + t.y)) / t.k;
    canvas.svg.transition().duration(duration)
        .call(canvas.zoom.translateBy, 0, dy);
}

// ===== Data table =====
const fixed = digits => v => (v === null || v === undefined ? '' : v.toFixed(digits));
/**
 * Columns of the atom table (data tab); also the column set of the "atoms" CSV
 * export. `value` extracts the cell from a graph node, `num` right-aligns it.
 * @const {!Array<{label: string, title: string, value: function(!Object): *,
 *     num: (boolean|undefined)}>}
 */
const DATA_COLUMNS = [
    { label: gettext('Record'),    title: gettext('Record-Typ (ATOM/HETATM)'), value: n => n.record },
    { label: gettext('Nr.'),       title: gettext('Atom-Seriennummer'),        value: n => n.serial, num: true },
    { label: gettext('Atom'),      title: gettext('Atomname'),                 value: n => n.label },
    { label: gettext('Alt'),       title: gettext('Alternate Location'),       value: n => n.alt_loc },
    { label: gettext('Residue'),   title: gettext('Residuename'),              value: n => n.residue },
    { label: gettext('Chain'),     title: gettext('Chain-ID'),                 value: n => n.chain },
    { label: gettext('Pos.'),      title: gettext('Residue-Sequenznummer'),    value: n => n.residue_seq, num: true },
    { label: gettext('iCode'),     title: gettext('Insertion Code'),           value: n => n.i_code },
    { label: 'x',         title: gettext('x-Koordinate (Å)'),         value: n => fixed(3)(n.x), num: true },
    { label: 'y',         title: gettext('y-Koordinate (Å)'),         value: n => fixed(3)(n.y), num: true },
    { label: 'z',         title: gettext('z-Koordinate (Å)'),         value: n => fixed(3)(n.z), num: true },
    { label: gettext('Besetzung'), title: gettext('Occupancy'),                value: n => fixed(2)(n.occupancy), num: true },
    { label: gettext('B-Faktor'),  title: gettext('Temperaturfaktor (Å²)'),    value: n => fixed(2)(n.temp_factor), num: true },
    { label: gettext('Element'),   title: gettext('Element'),                  value: n => n.element },
    { label: gettext('Ladung'),    title: gettext('Formalladung'),             value: n => n.charge },
];

/**
 * residueKey() computed from a graph node instead of a residue.
 * @param {!Object} n Graph node.
 * @return {string} `<chain>_<pos>_<name>`.
 */
const nodeResidueKey = n => `${n.chain}_${n.residue_seq}${n.i_code || ''}_${n.residue}`;

/**
 * Renders the atom table: one row per atom in PDB file order; the first row of
 * every residue gets a separator line. Stores `tableRows`, `tableScroll` and
 * `firstRowByResidue` (residue index -> first row element) on the canvas.
 * @param {!Object} canvas Canvas object.
 * @param {!Array<!Object>} nodes Graph nodes.
 */
function renderDataTable(canvas, nodes) {
    const table = d3.select('#data-table-scroll table');
    if (table.empty()) return;
    canvas.tableScroll = document.getElementById('data-table-scroll');

    table.select('thead tr').selectAll('th')
        .data(DATA_COLUMNS)
        .join('th')
        .attr('scope', 'col')
        .attr('class', c => c.num ? 'text-end' : null)
        .attr('title', c => c.title)
        .text(c => c.label);

    const residueStarts = new Set(nodes
        .filter((n, i) => i === 0 || nodes[i - 1].residue_index !== n.residue_index)
        .map(n => n.id));

    canvas.tableRows = table.select('tbody').selectAll('tr')
        .data(nodes, n => n.id)
        .join('tr')
        .classed('residue-start', n => residueStarts.has(n.id))
        .classed('hetatm', n => n.hetatm);

    canvas.tableRows.selectAll('td')
        .data(n => DATA_COLUMNS.map(c => c.value(n)))
        .join('td')
        .attr('class', (v, i) => DATA_COLUMNS[i].num ? 'text-end' : null)
        .text(v => v ?? '');

    canvas.firstRowByResidue = new Map();
    canvas.tableRows.each(function (n) {
        if (!canvas.firstRowByResidue.has(n.residue_index)) canvas.firstRowByResidue.set(n.residue_index, this);
    });

    updateSelection(canvas);
}

/**
 * Scrolls the atom table so that the first row of the residue sits just below the
 * sticky header.
 * @param {!Object} canvas Canvas object.
 * @param {!Object} res Residue.
 * @param {boolean} smooth Whether to scroll smoothly.
 */
function scrollTableToResidue(canvas, res, smooth) {
    const row = canvas.firstRowByResidue && canvas.firstRowByResidue.get(res.index);
    if (!row) return;
    const headH = canvas.tableScroll.querySelector('thead').offsetHeight;
    canvas.tableScroll.scrollTo({ top: row.offsetTop - headH, behavior: smooth ? 'smooth' : 'auto' });
}

/**
 * Highlights in the sequence bar the residues whose table rows are currently
 * visible (binary search over the row offsets).
 * @param {!Object} canvas Canvas object.
 */
function updateTableView(canvas) {
    const rows = canvas.tableRows ? canvas.tableRows.nodes() : [];
    if (rows.length === 0) return;
    const scroller = canvas.tableScroll;
    const top = scroller.scrollTop + scroller.querySelector('thead').offsetHeight;
    const bottom = scroller.scrollTop + scroller.clientHeight;

    // Rows are sorted by offsetTop: first row whose bottom edge lies below `y`
    // (1 px tolerance so that a row ending exactly at the edge does not count).
    const firstBelow = y => {
        let lo = 0, hi = rows.length;
        while (lo < hi) {
            const mid = (lo + hi) >> 1;
            if (rows[mid].offsetTop + rows[mid].offsetHeight > y + 1) hi = mid; else lo = mid + 1;
        }
        return lo;
    };
    const firstIdx = Math.min(firstBelow(top), rows.length - 1);
    const lastIdx = Math.max(Math.min(firstBelow(bottom), rows.length - 1), firstIdx);
    const first = d3.select(rows[firstIdx]).datum().residue_index;
    const last = d3.select(rows[lastIdx]).datum().residue_index;
    markSequenceInView(canvas, c => c.res.index >= first && c.res.index <= last);
}

// ===== Statistics (composition) =====
const formatInt = v => v.toLocaleString(LOCALE);
// Percentages (0–100) in the page language's format, e.g. "12,5 %" or "12.5%".
const formatPercent = v => (v / 100).toLocaleString(LOCALE, { style: 'percent', minimumFractionDigits: 1, maximumFractionDigits: 1 });
const formatPercentTick = v => (v / 100).toLocaleString(LOCALE, { style: 'percent', maximumFractionDigits: 2 });
// Element symbol in the usual spelling (CA → Ca); the colour comes from the atom
// CSS variables (as a var() reference, so it follows the theme without re-rendering).
const elementSymbol = el => el ? el[0].toUpperCase() + el.slice(1).toLowerCase() : '?';
const elementColor = el => `var(--atom-${(el || '').toUpperCase()}, var(--atom-default))`;

/** @const {!Object<string, string>} Translated names of the 20 standard amino acids. */
const AMINO_ACID_NAMES = {
    ALA: gettext('Alanin'), ARG: gettext('Arginin'), ASN: gettext('Asparagin'), ASP: gettext('Asparaginsäure'),
    CYS: gettext('Cystein'), GLN: gettext('Glutamin'), GLU: gettext('Glutaminsäure'), GLY: gettext('Glycin'),
    HIS: gettext('Histidin'), ILE: gettext('Isoleucin'), LEU: gettext('Leucin'), LYS: gettext('Lysin'),
    MET: gettext('Methionin'), PHE: gettext('Phenylalanin'), PRO: gettext('Prolin'), SER: gettext('Serin'),
    THR: gettext('Threonin'), TRP: gettext('Tryptophan'), TYR: gettext('Tyrosin'), VAL: gettext('Valin'),
};

/**
 * Display name of an amino acid; modified ones (e.g. MSE) are named after their
 * parent amino acid if known.
 * @param {string} name Three-letter residue name.
 * @param {?string=} parent Parent amino acid of a modified residue.
 * @return {string} Translated name, "<parent> (modified)" or the raw name.
 */
function aminoAcidLabel(name, parent) {
    if (AMINO_ACID_NAMES[name]) return AMINO_ACID_NAMES[name];
    if (parent && AMINO_ACID_NAMES[parent]) return interpolate(gettext('%s (modifiziert)'), [AMINO_ACID_NAMES[parent]]);
    return name;
}
// Parent amino acid per residue name (modified residues only)
const parentsByName = residues => new Map(residues.filter(r => r.parent).map(r => [r.name, r.parent]));

/** @const {!Object<string, string>} Translated labels of the edge types. */
const BOND_TYPE_LABELS = {
    backbone: gettext('Backbone'), sidechain: gettext('Seitenkette'), peptide: gettext('Peptidbindung'),
    disulfide: gettext('Disulfidbrücke'), ligand: gettext('Ligand'),
};

/**
 * Renders a "count and share" table with a bar per row; sorted by descending
 * count by default.
 * @param {string} tbodySelector Selector of the target `<tbody>`.
 * @param {!Array<{key: string, label: string, sub: (string|undefined),
 *     count: number, color: string, dot: (boolean|undefined)}>} rows Rows;
 *     `dot` shows a colour dot before the label, `sub` a secondary text.
 * @param {number} total Reference value for the share column.
 * @param {{sort: (boolean|undefined)}=} options `sort: false` keeps the row order.
 */
function renderShareTable(tbodySelector, rows, total, { sort = true } = {}) {
    if (sort) rows.sort((a, b) => d3.descending(a.count, b.count) || d3.ascending(a.label, b.label));
    const maxCount = d3.max(rows, r => r.count) || 1;

    d3.select(tbodySelector).selectAll('tr')
        .data(rows, r => r.key)
        .join(enter => {
            const tr = enter.append('tr');
            tr.append('td')
                .call(td => td.append('span').attr('class', 'legend-dot me-2 align-middle'))
                .call(td => td.append('span').attr('class', 'share-label'))
                .call(td => td.append('span').attr('class', 'share-sub text-body-secondary small ms-2'));
            tr.append('td').attr('class', 'share-count text-end');
            tr.append('td').attr('class', 'share-pct text-end');
            tr.append('td').attr('class', 'bar-cell').append('div').attr('class', 'stats-bar');
            return tr;
        })
        .order()
        .call(tr => tr.select('.legend-dot')
            .style('display', r => r.dot ? null : 'none')
            .style('background', r => r.color))
        .call(tr => tr.select('.share-label').text(r => r.label))
        .call(tr => tr.select('.share-sub').text(r => r.sub || ''))
        .call(tr => tr.select('.share-count').text(r => formatInt(r.count)))
        .call(tr => tr.select('.share-pct').text(r => formatPercent(100 * r.count / total)))
        .call(tr => tr.select('.stats-bar')
            .style('width', r => `${100 * r.count / maxCount}%`)
            .style('background', r => r.color));
}

/**
 * Renders the composition part of the statistics tab: counts, elements, amino
 * acids and the per-chain breakdowns.
 * @param {!Object} graphData Graph data.
 */
function renderStats(graphData) {
    if (d3.select('#stat-atoms').empty()) return;
    const { residues, nodes } = graphData;

    // Count amino acids only (including modified ones with a known one-letter
    // code), no ligands, ions or water.
    const aminoAcids = residues.filter(r => r.one_letter !== '?');
    d3.select('#stat-amino-acids').text(formatInt(aminoAcids.length));
    d3.select('#stat-atoms').text(formatInt(nodes.length));
    d3.selectAll('.share-total-pct').text(formatPercent(100));

    renderShareTable('#stat-elements',
        d3.rollups(nodes, v => v.length, n => n.element || '').map(([el, count]) => ({
            key: el, label: elementSymbol(el), count, color: elementColor(el), dot: true,
        })),
        nodes.length);

    const parents = parentsByName(aminoAcids);
    renderShareTable('#stat-amino-acid-counts',
        d3.rollups(aminoAcids, v => v.length, r => r.name).map(([name, count]) => ({
            key: name,
            label: aminoAcidLabel(name, parents.get(name)),
            sub: `${name} · ${aminoAcids.find(r => r.name === name).one_letter}`,
            count,
            color: 'var(--bs-primary)',
        })),
        aminoAcids.length);

    renderChainStats(nodes, aminoAcids);
}

const chainLabel = c => interpolate(gettext('Kette %s'), [c.trim() || '–']);

/**
 * "Per chain" section: a matrix (rows × chains) with shares relative to each
 * chain, plus the distribution of the counted objects over the chains (absolute
 * and relative).
 * @param {{matrix: string, share: string, shareTotal: string, items: !Array<!Object>,
 *     keyOf: function(!Object): string, rowLabel: string,
 *     renderLabel: function(!d3.Selection, string)}} config `matrix`/`share`/
 *     `shareTotal`: selectors of the target elements; `items`: counted objects
 *     (atoms or amino acids) with `chain`; `keyOf`: row key; `rowLabel`: header
 *     of the first column; `renderLabel(th, key)`: renders a row label.
 */
function renderChainBreakdown({ matrix, share, shareTotal, items, keyOf, rowLabel, renderLabel }) {
    const table = d3.select(matrix);
    if (table.empty()) return;

    const chains = [...new Set(items.map(n => n.chain))]; // Order as in the file
    const chainTotals = d3.rollup(items, v => v.length, n => n.chain);
    const byChain = d3.rollup(items, v => v.length, n => n.chain, keyOf);
    const keys = d3.rollups(items, v => v.length, keyOf)
        .sort((a, b) => d3.descending(a[1], b[1]) || d3.ascending(a[0], b[0]))
        .map(([key]) => key);

    // Table 1: shares are relative to the respective chain.
    table.selectAll('*').remove();
    const head = table.append('thead');
    const groupRow = head.append('tr');
    groupRow.append('th').attr('scope', 'col').attr('rowspan', 2).attr('class', 'align-bottom')
        .text(rowLabel);
    chains.forEach(c => groupRow.append('th')
        .attr('scope', 'colgroup').attr('colspan', 2).attr('class', 'text-center chain-start')
        .text(chainLabel(c)));
    const subRow = head.append('tr');
    chains.forEach(() => {
        subRow.append('th').attr('scope', 'col').attr('class', 'text-end chain-start').text(gettext('Anzahl'));
        subRow.append('th').attr('scope', 'col').attr('class', 'text-end').text(gettext('Anteil'));
    });

    const body = table.append('tbody');
    keys.forEach(key => {
        const tr = body.append('tr');
        renderLabel(tr.append('th').attr('scope', 'row'), key);
        chains.forEach(c => {
            const count = byChain.get(c).get(key) || 0;
            tr.append('td').attr('class', 'text-end chain-start').text(formatInt(count));
            tr.append('td').attr('class', 'text-end').text(formatPercent(100 * count / chainTotals.get(c)));
        });
    });

    const total = table.append('tfoot').append('tr').attr('class', 'share-total');
    total.append('th').attr('scope', 'row').text(gettext('Gesamt'));
    chains.forEach(c => {
        const sum = d3.sum(keys, key => byChain.get(c).get(key) || 0);
        total.append('td').attr('class', 'text-end chain-start').text(formatInt(sum));
        total.append('td').attr('class', 'text-end').text(formatPercent(100 * sum / chainTotals.get(c)));
    });

    // Table 2: distribution over the chains.
    renderShareTable(share,
        chains.map(c => ({ key: c, label: chainLabel(c), count: chainTotals.get(c), color: 'var(--bs-primary)' })),
        items.length,
        { sort: false });
    d3.select(shareTotal).text(formatInt(items.length));
}

/**
 * Renders the per-chain breakdowns of elements and amino acids.
 * @param {!Array<!Object>} nodes Graph nodes (atoms).
 * @param {!Array<!Object>} aminoAcids Amino-acid residues.
 */
function renderChainStats(nodes, aminoAcids) {
    renderChainBreakdown({
        matrix: '#stat-elements-by-chain', share: '#stat-chain-atoms', shareTotal: '#stat-chain-atoms-total',
        items: nodes,
        keyOf: n => n.element || '',
        rowLabel: gettext('Element'),
        renderLabel: (th, el) => {
            th.append('span').attr('class', 'legend-dot me-2 align-middle').style('background', elementColor(el));
            th.append('span').text(elementSymbol(el));
        },
    });
    // Only chains with amino acids appear (pure ligand chains would have no shares).
    renderChainBreakdown({
        matrix: '#stat-aa-by-chain', share: '#stat-chain-aa', shareTotal: '#stat-chain-aa-total',
        items: aminoAcids,
        keyOf: r => r.name,
        rowLabel: gettext('Aminosäure'),
        renderLabel: (th, name) => {
            th.append('span').text(aminoAcidLabel(name, parentsByName(aminoAcids).get(name)));
            th.append('span').attr('class', 'text-body-secondary small ms-2')
                .text(`${name} · ${aminoAcids.find(r => r.name === name).one_letter}`);
        },
    });
}

// ===== Legend =====
// Elements of the displayed model: the main protein elements first, then alphabetically.
const LEGEND_ORDER = ['C', 'N', 'O', 'S', 'H', 'P'];

/**
 * Renders the element colour legend of the graph tab (`#legend-elements`).
 * @param {!Array<!Object>} nodes Graph nodes.
 */
function renderLegend(nodes) {
    const rank = el => { const i = LEGEND_ORDER.indexOf(el.toUpperCase()); return i < 0 ? LEGEND_ORDER.length : i; };
    const elements = [...new Set(nodes.map(n => n.element || ''))]
        .sort((a, b) => d3.ascending(rank(a), rank(b)) || d3.ascending(elementSymbol(a), elementSymbol(b)));

    d3.select('#legend-elements').selectAll('.legend-item')
        .data(elements, el => el)
        .join(enter => {
            const item = enter.append('span').attr('class', 'legend-item d-flex align-items-center gap-1');
            item.append('span').attr('class', 'legend-dot');
            item.append('span').attr('class', 'legend-label');
            return item;
        })
        .order()
        .call(item => item.select('.legend-dot').style('background', elementColor))
        .call(item => item.select('.legend-label').text(elementSymbol));
}

// ===== Reference values (Engh & Huber 2001 for the backbone, CCP4 Monomer Library for side chains) =====
// ref = [mean, standard deviation, group] from the server, or null (no reference).
// Bonds carry it as `edge.ref`, bond angles get it from `node.angle_refs` of the
// centre atom (see analyzeAngles). Group: "backbone:N-CA:GLY" (section, atoms and
// class) or "LYS:CB-CG" (side chain of an amino acid, no class).
/** @const {number} |Z| above which a value is flagged as a warning. */
const Z_WARN = 2;
/** @const {number} |Z| above which a value counts as an outlier. */
const Z_OUTLIER = 4;
const zOf = (value, ref) => ref ? (value - ref[0]) / ref[1] : null;
const zScore = e => zOf(e.length, e.ref);
// CSS class suffix (z-ok, z-warn, z-outlier, z-none) and colour key of a Z-score.
const zClass = z => z === null ? 'none' : Math.abs(z) > Z_OUTLIER ? 'outlier' : Math.abs(z) > Z_WARN ? 'warn' : 'ok';
const formatZ = z => z.toLocaleString(LOCALE, { minimumFractionDigits: 1, maximumFractionDigits: 1, signDisplay: 'exceptZero' });
const formatRmsz = v => v === null ? '–' : v.toLocaleString(LOCALE, { minimumFractionDigits: 2, maximumFractionDigits: 2 });

// Residue classes of backbone references; classes across the peptide bond use X
// for the respective other amino acid.
const REFERENCE_CLASS_LABELS = {
    GLY: 'Gly', PRO: 'Pro', 'X-GLY': 'X–Gly', 'X-PRO': 'X–Pro', 'GLY-X': 'Gly–X', 'cis-X-PRO': 'X–Pro (cis)',
};
/**
 * Splits a reference group into display parts.
 * @param {string} group E.g. "backbone:C-N:X-PRO".
 * @return {!Array<string>} [atoms, class label, section], e.g.
 *     ["C–N", "X–Pro", "backbone"]; the class label is '' if there is none.
 */
function referenceGroupLabel(group) {
    const [section, atoms, cls] = group.split(':');
    const label = !cls ? '' : cls === 'other' ? gettext('übrige') : REFERENCE_CLASS_LABELS[cls] || cls;
    return [atoms.replaceAll('-', '–'), label, section];
}
const capitalize = s => s.charAt(0).toUpperCase() + s.slice(1);
const referenceSectionLabel = (section, parents) =>
    section === 'backbone' ? capitalize(gettext('Backbone')) : aminoAcidLabel(section, parents.get(section));
// Source note below reference tables; `sources` = structure-info's reference_sources.
const referenceSourceText = sources => sources
    ? interpolate(gettext('Referenzwerte: %(backbone)s (Backbone), %(sidechain)s (Seitenketten)'), sources, true) : '';

/**
 * Tooltip line describing a bond's reference value and Z-score.
 * @param {!Object} e Edge.
 * @return {string} Translated text.
 */
function referenceText(e) {
    if (!e.ref) return gettext('Ohne Referenzwert');
    return interpolate(gettext('Referenz %(ref)s ± %(sigma)s Å · Z = %(z)s'), {
        ref: formatLength(e.ref[0]), sigma: formatLength(e.ref[1]), z: formatZ(zScore(e)),
    }, true);
}

/**
 * Rates bonds or angles against their reference values: overall, per group, per
 * section (backbone first, then the side chains per amino acid) and the outliers
 * sorted by descending |Z|. RMSZ = root mean square of the Z-scores.
 * @param {!Array<{value: number, ref: ?Array, z: ?number, atoms: !Array}>} items
 *     Bonds or angles.
 * @return {{total: number, rated: number, rmsz: ?number, groups: !Array<!Object>,
 *     sections: !Array<!Object>, outliers: !Array<!Object>}} The rating.
 */
function rateReferences(items) {
    const rated = items.filter(i => i.z !== null);
    const rmsz = list => list.length ? Math.sqrt(d3.mean(list, i => i.z * i.z)) : null;
    const outliersOf = list => list.filter(i => Math.abs(i.z) > Z_OUTLIER).length;
    const groups = d3.groups(rated, i => i.ref[2])
        .map(([group, list]) => ({
            group, section: group.split(':')[0], ref: list[0].ref[0], sigma: list[0].ref[1], count: list.length,
            mean: d3.mean(list, i => i.value), rmsz: rmsz(list), outliers: outliersOf(list),
        }))
        .sort((a, b) => d3.descending(a.count, b.count) || d3.ascending(a.group, b.group));
    // Within a section, classes of the same atoms stay together (atoms ordered by total count)
    const atomsOf = g => g.group.split(':')[1];
    const sectionGroups = section => {
        const list = groups.filter(g => g.section === section);
        const perAtoms = d3.rollup(list, v => d3.sum(v, g => g.count), atomsOf);
        return list.sort((a, b) => d3.descending(perAtoms.get(atomsOf(a)), perAtoms.get(atomsOf(b)))
            || d3.ascending(atomsOf(a), atomsOf(b)) || d3.descending(a.count, b.count));
    };
    const sections = d3.groups(rated, i => i.ref[2].split(':')[0])
        .map(([section, list]) => ({
            section, count: list.length, rmsz: rmsz(list), outliers: outliersOf(list),
            groups: sectionGroups(section),
        }))
        .sort((a, b) => d3.descending(a.section === 'backbone', b.section === 'backbone')
            || d3.descending(a.count, b.count) || d3.ascending(a.section, b.section));
    return {
        total: items.length, rated: rated.length, rmsz: rmsz(rated), groups, sections,
        outliers: rated.filter(i => Math.abs(i.z) > Z_OUTLIER).sort((a, b) => d3.descending(Math.abs(a.z), Math.abs(b.z))),
    };
}

/**
 * Reference rating of all bonds of the graph.
 * @param {!Object} graphData Graph data.
 * @return {!Object} Result of rateReferences.
 */
function referenceStats(graphData) {
    const nodeById = new Map(graphData.nodes.map(n => [n.id, n]));
    return rateReferences(graphData.edges
        .filter(e => nodeById.has(e.source) && nodeById.has(e.target))
        .map(e => ({ ...e, value: e.length, atoms: [nodeById.get(e.source), nodeById.get(e.target)], z: zScore(e) })));
}

// Labels and formatters for bond lengths ('bond') or bond angles ('angle'); lets
// the reference tables (tab and PDF) share one implementation.
const referenceKind = kind => kind === 'angle' ? {
    item: gettext('Winkel'), rated: gettext('Bewertete Winkel'), value: gettext('Winkel (°)'),
    reference: gettext('Referenz (°)'), mean: gettext('Durchschnitt (°)'), delta: gettext('Abw. (°)'),
    format: formatDegrees, formatDelta: formatDegreesDelta,
    describe: atoms => ({ ...describeAngle(atoms), atoms: atoms.map(a => a.label).join('–') }),
} : {
    item: gettext('Bindung'), rated: gettext('Bewertete Bindungen'), value: gettext('Länge (Å)'),
    reference: gettext('Referenz (Å)'), mean: gettext('Durchschnitt (Å)'), delta: gettext('Abw. (Å)'),
    format: formatLength, formatDelta, describe: describeBond,
};

// Expanded sections ("<prefix>|<section>"); kept across re-renders.
const expandedReferenceSections = new Set();

/**
 * Renders a reference rating: key figures, a table per section (expandable to
 * the groups) and the outliers, into the elements #<prefix>-facts, -groups,
 * -outliers, -outliers-block and -outliers-title.
 * @param {string} prefix Element id prefix ('stat-ref' or 'angle-ref').
 * @param {!Object} s Result of rateReferences.
 * @param {string} kind 'bond' or 'angle' (see referenceKind).
 * @param {!Map<string, string>} parents Parent amino acid by residue name.
 */
function renderReferenceStats(prefix, s, kind, parents) {
    const facts = d3.select(`#${prefix}-facts`);
    if (facts.empty()) return;
    const k = referenceKind(kind);

    facts.selectAll('div')
        .data([
            [k.rated, interpolate(gettext('%(rated)s von %(total)s (%(share)s)'), {
                rated: formatInt(s.rated), total: formatInt(s.total), share: formatPercent(s.total ? 100 * s.rated / s.total : 0),
            }, true)],
            [gettext('RMSZ'), formatRmsz(s.rmsz)],
            [interpolate(gettext('Ausreißer (|Z| > %s)'), [Z_OUTLIER]), formatInt(s.outliers.length)],
        ])
        .join(enter => {
            const div = enter.append('div');
            div.append('dt');
            div.append('dd');
            return div;
        })
        .call(div => div.select('dt').text(d => d[0]))
        .call(div => div.select('dd').text(d => d[1]));

    const table = d3.select(`#${prefix}-groups`);
    table.selectAll('*').remove();
    const columns = [
        [k.item], [gettext('Klasse')], [k.reference, true], [gettext('Anzahl'), true],
        [k.mean, true], [k.delta, true], [gettext('RMSZ'), true], [gettext('Ausreißer'), true],
    ];
    table.append('thead').append('tr').selectAll('th').data(columns).join('th')
        .attr('scope', 'col').attr('class', c => c[1] ? 'text-end' : null).text(c => c[0]);
    const rmszCell = (row, value) => row.append('td').attr('class', `text-end z-${zClass(value)}`).text(formatRmsz(value));
    const outlierCell = (row, count) => row.append('td').attr('class', `text-end${count ? ' z-outlier' : ''}`).text(formatInt(count));
    s.sections.forEach(section => {
        const key = `${prefix}|${section.section}`;
        const body = table.append('tbody').attr('class', 'bond-group');
        const summary = body.append('tr').attr('class', 'bond-summary');
        const toggle = summary.append('th').attr('scope', 'row').attr('colspan', 3)
            .append('button').attr('type', 'button').attr('class', 'bond-toggle').html(CHEVRON_SVG);
        toggle.append('span').text(referenceSectionLabel(section.section, parents));
        summary.append('td').attr('class', 'text-end').text(formatInt(section.count));
        summary.append('td').attr('colspan', 2);
        rmszCell(summary, section.rmsz);
        outlierCell(summary, section.outliers);

        section.groups.forEach(g => {
            const [atoms, cls] = referenceGroupLabel(g.group);
            const tr = body.append('tr').attr('class', 'ref-detail');
            tr.append('th').attr('scope', 'row').text(atoms);
            tr.append('td').text(cls);
            tr.append('td').attr('class', 'text-end').text(`${k.format(g.ref)} ± ${k.format(g.sigma)}`);
            tr.append('td').attr('class', 'text-end').text(formatInt(g.count));
            tr.append('td').attr('class', 'text-end').text(k.format(g.mean));
            tr.append('td').attr('class', 'text-end').text(k.formatDelta(g.mean - g.ref));
            rmszCell(tr, g.rmsz);
            outlierCell(tr, g.outliers);
        });
        const details = body.selectAll('tr.ref-detail');
        const setExpanded = expanded => {
            toggle.attr('aria-expanded', expanded);
            details.property('hidden', !expanded);
            if (expanded) expandedReferenceSections.add(key);
            else expandedReferenceSections.delete(key);
        };
        summary.on('click', () => setExpanded(!expandedReferenceSections.has(key)));
        setExpanded(expandedReferenceSections.has(key));
    });
    const total = table.append('tfoot').append('tr').attr('class', 'share-total');
    total.append('th').attr('scope', 'row').attr('colspan', 3).text(gettext('Gesamt'));
    total.append('td').attr('class', 'text-end').text(formatInt(s.rated));
    total.append('td').attr('colspan', 2);
    rmszCell(total, s.rmsz);
    outlierCell(total, s.outliers.length);

    // Individual outliers, by descending |Z|
    d3.select(`#${prefix}-outliers-title`).text(interpolate(gettext('Ausreißer (|Z| > %s)'), [Z_OUTLIER]));
    const list = d3.select(`#${prefix}-outliers`);
    list.selectAll('*').remove();
    d3.select(`#${prefix}-outliers-block`).property('hidden', !s.outliers.length);
    if (!s.outliers.length) return;
    list.append('thead').append('tr').selectAll('th')
        .data([[gettext('Atome')], [gettext('Kette')], [gettext('Pos.')], [gettext('Aminosäure')],
               [k.value, true], [k.reference, true], ['Z', true]])
        .join('th').attr('scope', 'col').attr('class', c => c[1] ? 'text-end' : null).text(c => c[0]);
    const body = list.append('tbody');
    s.outliers.forEach(item => {
        const d = k.describe(item.atoms);
        const tr = body.append('tr');
        [d.atoms, d.chain, d.pos, d.residues].forEach(v => tr.append('td').text(v));
        tr.append('td').attr('class', 'text-end').text(k.format(item.value));
        tr.append('td').attr('class', 'text-end').text(`${k.format(item.ref[0])} ± ${k.format(item.ref[1])}`);
        tr.append('td').attr('class', 'text-end z-outlier').text(formatZ(item.z));
    });
}

// Reference rating for the JSON/PDF export; same structure for bonds and angles.
// `digits` = decimals of the values (3 for Å, 1 for °).
const referenceExport = (s, unit, digits) => ({
    unit,
    sources: structureInfo()?.reference_sources ?? null,
    outlier_threshold: Z_OUTLIER,
    rated: s.rated,
    rmsz: roundTo(s.rmsz, 3),
    outlier_count: s.outliers.length,
    // section: "backbone" or the residue name of the side chain
    sections: s.sections.map(sec => ({
        section: sec.section, count: sec.count, rmsz: roundTo(sec.rmsz, 3), outliers: sec.outliers,
        groups: sec.groups.map(g => ({
            atoms: referenceGroupLabel(g.group)[0], class: referenceGroupLabel(g.group)[1] || null,
            reference: g.ref, sigma: g.sigma, count: g.count, mean: roundTo(g.mean, digits),
            rmsz: roundTo(g.rmsz, 3), outliers: g.outliers,
        })),
    })),
    outliers: s.outliers.map(i => ({
        value: roundTo(i.value, digits), reference: i.ref[0], sigma: i.ref[1], z: roundTo(i.z, 2),
        atoms: i.atoms.map(exportAtom),
    })),
});

// ===== Angles: bond angles, torsion angles φ/ψ and ω =====
// Bond angles and ω are computed in the browser from the atom coordinates of the
// displayed model (or the model average); φ/ψ and their Ramachandran evaluation
// come precomputed from the server (`residue.rama`).
const CIS_MAX_OMEGA = 30;     // |ω| < 30°: cis peptide bond
const TRANS_MIN_OMEGA = 150;  // |ω| ≥ 150°: trans; in between: twisted
const BACKBONE_ORDER = ['N', 'CA', 'C', 'O', 'OXT'];
const DEG = 180 / Math.PI;

// Minimal 3D vector helpers on atoms ({x, y, z}) and arrays ([x, y, z]).
const vec = (a, b) => [b.x - a.x, b.y - a.y, b.z - a.z];
const dot = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
const cross = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];

/**
 * Bond angle a–b–c.
 * @param {{x: number, y: number, z: number}} a First atom.
 * @param {{x: number, y: number, z: number}} b Centre atom.
 * @param {{x: number, y: number, z: number}} c Third atom.
 * @return {number} Angle in degrees (0 … 180).
 */
function bondAngle(a, b, c) {
    const u = vec(b, a), v = vec(b, c);
    const cos = dot(u, v) / Math.sqrt(dot(u, u) * dot(v, v));
    return Math.acos(Math.max(-1, Math.min(1, cos))) * DEG;
}

/**
 * Torsion (dihedral) angle p0–p1–p2–p3 around the axis p1–p2.
 * @param {{x: number, y: number, z: number}} p0 First atom.
 * @param {{x: number, y: number, z: number}} p1 Second atom.
 * @param {{x: number, y: number, z: number}} p2 Third atom.
 * @param {{x: number, y: number, z: number}} p3 Fourth atom.
 * @return {number} Angle in degrees (−180 … 180, IUPAC convention).
 */
function dihedral(p0, p1, p2, p3) {
    const b0 = vec(p1, p0), b2 = vec(p2, p3);
    let b1 = vec(p1, p2);
    const len = Math.sqrt(dot(b1, b1));
    b1 = b1.map(v => v / len);
    const v = b0.map((c, i) => c - dot(b0, b1) * b1[i]);
    const w = b2.map((c, i) => c - dot(b2, b1) * b1[i]);
    return Math.atan2(dot(cross(b1, v), w), dot(v, w)) * DEG;
}

// Rounded to 0.1°; `|| 0` turns −0 into 0 (otherwise "−0.0°")
const roundAngle = v => Math.round(v * 10) / 10 || 0;
const formatAngle = v => `${roundAngle(v).toLocaleString(LOCALE, { minimumFractionDigits: 1, maximumFractionDigits: 1 })}°`;

/**
 * Colour class in the Ramachandran plot; modified amino acids count like their
 * parent amino acid. (Not to be confused with the Top8000 category `rama.type`.)
 * @param {!Object} res Residue.
 * @return {string} 'gly', 'pro' or 'general'.
 */
function ramaClass(res) {
    const name = res.parent || res.name;
    return name === 'GLY' ? 'gly' : name === 'PRO' ? 'pro' : 'general';
}
const RAMA_CLASSES = ['general', 'gly', 'pro'];
const ramaClassLabel = cls => ({ general: gettext('Übrige Aminosäuren'), gly: gettext('Glycin'), pro: gettext('Prolin') })[cls];

// Evaluation from the server (Top8000 as in MolProbity): category and region per amino acid
const RAMA_TYPES = ['general', 'glycine', 'cis-proline', 'trans-proline', 'pre-proline', 'ile-val'];
const ramaTypeLabel = type => ({
    general: gettext('Allgemein'), glycine: gettext('Glycin'), 'cis-proline': gettext('cis-Prolin'),
    'trans-proline': gettext('trans-Prolin'), 'pre-proline': gettext('Vor Prolin'), 'ile-val': gettext('Ile und Val'),
})[type];
const RAMA_REGIONS = ['favored', 'allowed', 'outlier'];
const ramaRegionLabel = region => ({
    favored: gettext('Bevorzugt'), allowed: gettext('Erlaubt'), outlier: gettext('Ausreißer'),
})[region];
/**
 * Ramachandran view state: category shown in the plot ('' = all) and the region
 * grids per category (loaded on demand by loadRamaRegions).
 */
const ramaView = { type: '', regions: null };

/**
 * Collects φ/ψ per amino acid (from the server), ω per peptide bond and all bond
 * angles. Torsions across a peptide bond exist only where the server detected a
 * bond (no chain break).
 * @param {!Object} graphData Graph data.
 * @return {{ramachandran: !Array<!Object>, omegas: !Array<!Object>,
 *     bondAngles: !Array<!Object>, angleTypes: !Array<!Object>,
 *     byCenter: !Map<string, !Array<!Object>>}} `ramachandran`: points
 *     {res, phi, psi, cls, type, region}; `omegas`: {prev, res, omega, kind} with
 *     kind 'cis'/'twisted'/'trans'; `bondAngles`: {atoms, type, value, ref, z};
 *     `angleTypes`: bond angles grouped by type with statistics; `byCenter`:
 *     bond angles by centre atom id.
 */
function analyzeAngles(graphData) {
    const nodes = new Map(graphData.nodes.map(n => [n.id, n]));
    const peptides = new Set(graphData.edges.filter(e => e.type === 'peptide').map(edgeKey));
    const atom = (res, name) => nodes.get(`${res.chain}_${res.pos}_${name}`);
    const linked = (a, b) => {
        const c = atom(a, 'C'), n = atom(b, 'N');
        return Boolean(c && n && peptides.has(`${c.id}|${n.id}`));
    };

    const residues = graphData.residues;
    // φ/ψ, category and region are computed by the server (utils.evaluate_ramachandran)
    const ramachandran = residues.filter(res => res.rama).map(res => ({
        res, phi: res.rama.phi, psi: res.rama.psi, cls: ramaClass(res), type: res.rama.type, region: res.rama.region,
    }));
    const omegas = [];
    residues.forEach((res, i) => {
        const next = residues[i + 1];
        const ca = atom(res, 'CA'), c = atom(res, 'C');
        if (!ca || !c) return;
        if (next && linked(res, next) && atom(next, 'CA')) {
            const omega = dihedral(ca, c, atom(next, 'N'), atom(next, 'CA'));
            const abs = Math.abs(omega);
            const kind = abs < CIS_MAX_OMEGA ? 'cis' : abs < TRANS_MIN_OMEGA ? 'twisted' : 'trans';
            omegas.push({ prev: res, res: next, omega, kind });
        }
    });

    // Bond angles: per atom all pairs of its bonding partners. Order of the end
    // atoms: earlier residue first (C–N–CA across the peptide bond), within a
    // residue by backbone order (N, CA, C, O, OXT), then by name; the partner of a
    // disulfide bridge always comes last (CB–SG–SG). The reference is looked up in
    // the centre atom's `angle_refs` by the sorted end-atom names.
    const neighbors = new Map();
    const disulfide = new Set();
    graphData.edges.forEach(e => {
        if (!nodes.has(e.source) || !nodes.has(e.target)) return;
        [[e.source, e.target], [e.target, e.source]].forEach(([a, b]) => {
            neighbors.set(a, [...(neighbors.get(a) || []), b]);
            if (e.type === 'disulfide') disulfide.add(`${a}|${b}`);
        });
    });
    const priority = n => { const i = BACKBONE_ORDER.indexOf(n.label); return i < 0 ? BACKBONE_ORDER.length : i; };
    const rank = (center, n) => [disulfide.has(`${center.id}|${n.id}`) ? 1 : 0, n.residue_index, priority(n), n.label];
    const compare = (center, a, b) => {
        const ra = rank(center, a), rb = rank(center, b);
        for (let i = 0; i < ra.length; i++) if (ra[i] !== rb[i]) return ra[i] < rb[i] ? -1 : 1;
        return 0;
    };
    const bondAngles = [];
    const byCenter = new Map();
    neighbors.forEach((ids, centerId) => {
        const center = nodes.get(centerId);
        for (let i = 0; i < ids.length; i++) {
            for (let j = i + 1; j < ids.length; j++) {
                const [a, c] = [nodes.get(ids[i]), nodes.get(ids[j])].sort((x, y) => compare(center, x, y));
                const value = bondAngle(a, center, c);
                const ref = center.angle_refs?.[[a.label, c.label].sort().join('|')] ?? null;
                const angle = { atoms: [a, center, c], type: `${a.label}–${center.label}–${c.label}`, value, ref, z: zOf(value, ref) };
                bondAngles.push(angle);
                byCenter.set(centerId, [...(byCenter.get(centerId) || []), angle]);
            }
        }
    });
    const angleTypes = d3.groups(bondAngles, a => a.type)
        .map(([type, items]) => ({
            type, items, count: items.length, mean: d3.mean(items, a => a.value),
            sd: items.length > 1 ? d3.deviation(items, a => a.value) : null,
            min: d3.min(items, a => a.value), max: d3.max(items, a => a.value),
        }))
        .sort((a, b) => d3.descending(a.count, b.count) || d3.ascending(a.type, b.type));

    return { ramachandran, omegas, bondAngles, angleTypes, byCenter };
}

const RAMA_SIZE = 480;
const RAMA_MAJOR_STEP = 60; // Major grid and axis labels (°)
const RAMA_MINOR_STEP = 20; // Minor grid (°)
const RAMA_MARGIN = { top: 12, right: 14, bottom: 46, left: 56 };

/**
 * Draws a Ramachandran plot into `svg` (used by the tab and the image export).
 * Colours come from CSS classes (.rama-general, .rama-gly, .rama-pro, ...). If a
 * category is selected and its regions are loaded, the favoured/allowed regions
 * are drawn as contours and only points of that category are shown.
 * @param {!d3.Selection} svg Target SVG selection.
 * @param {!Array<!Object>} points Ramachandran points from analyzeAngles.
 * @param {number} size Width and height of the plot (px).
 * @param {{type: string, regions: ?Object}=} view View state (default ramaView).
 * @return {!d3.Selection} Selection of the point circles.
 */
function drawRamachandran(svg, points, size, view = ramaView) {
    const m = RAMA_MARGIN;
    const x = d3.scaleLinear().domain([-180, 180]).range([m.left, size - m.right]);
    const y = d3.scaleLinear().domain([-180, 180]).range([size - m.bottom, m.top]);
    const ticks = d3.range(-180, 181, RAMA_MAJOR_STEP);
    const degrees = d => `${d}°`;

    // Allowed and favoured region of the selected category: a string with one digit
    // per grid cell (`bins` × `bins`, currently 180 × 180 at 2°, index φ-bin * bins +
    // ψ-bin; 0 = outlier, 1 = allowed, 2 = favoured), transposed into the row-major
    // layout d3.contours expects and projected into plot coordinates.
    const regions = view.type && view.regions?.regions[view.type];
    if (regions) {
        const n = view.regions.bins;
        const values = new Array(n * n);
        for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) values[j * n + i] = +regions[i * n + j];
        const step = 360 / n;
        const path = d3.geoPath(d3.geoTransform({
            point(px, py) { this.stream.point(x(-180 + px * step), y(-180 + py * step)); },
        }));
        svg.append('g').selectAll('path')
            .data(d3.contours().size([n, n]).thresholds([1, 2])(values))
            .join('path')
            .attr('class', c => c.value === 2 ? 'rama-region-favored' : 'rama-region-allowed')
            .attr('d', path);
        points = points.filter(p => p.type === view.type);
    }

    // Minor grid (finer than the major grid), below everything else
    const minor = d3.range(-180, 181, RAMA_MINOR_STEP).filter(v => v % RAMA_MAJOR_STEP !== 0);
    const minorGrid = svg.append('g').attr('class', 'rama-minor');
    minorGrid.selectAll('line.v').data(minor).join('line').attr('class', 'v')
        .attr('x1', x).attr('x2', x).attr('y1', y(-180)).attr('y2', y(180));
    minorGrid.selectAll('line.h').data(minor).join('line').attr('class', 'h')
        .attr('x1', x(-180)).attr('x2', x(180)).attr('y1', y).attr('y2', y);

    const xAxis = svg.append('g').attr('class', 'hist-axis').attr('transform', `translate(0,${y(-180)})`)
        .call(d3.axisBottom(x).tickValues(ticks).tickSize(-(y(-180) - y(180))).tickPadding(8).tickFormat(degrees));
    const yAxis = svg.append('g').attr('class', 'hist-axis').attr('transform', `translate(${x(-180)},0)`)
        .call(d3.axisLeft(y).tickValues(ticks).tickSize(-(x(180) - x(-180))).tickPadding(8).tickFormat(degrees));
    [xAxis, yAxis].forEach(axis => axis.selectAll('.tick line').attr('class', 'hist-grid'));
    svg.append('rect').attr('class', 'rama-frame').attr('fill', 'none')
        .attr('x', x(-180)).attr('y', y(180)).attr('width', x(180) - x(-180)).attr('height', y(-180) - y(180));
    svg.append('line').attr('class', 'rama-zero').attr('x1', x(0)).attr('x2', x(0)).attr('y1', y(-180)).attr('y2', y(180));
    svg.append('line').attr('class', 'rama-zero').attr('x1', x(-180)).attr('x2', x(180)).attr('y1', y(0)).attr('y2', y(0));

    svg.append('text').attr('class', 'hist-label').attr('text-anchor', 'middle')
        .attr('x', (x(-180) + x(180)) / 2).attr('y', size - 8).text('φ');
    svg.append('text').attr('class', 'hist-label').attr('text-anchor', 'middle')
        .attr('transform', `translate(16,${(y(-180) + y(180)) / 2}) rotate(-90)`).text('ψ');

    // Draw glycine and proline above the others and outliers on top so that they stay visible
    const order = p => RAMA_CLASSES.indexOf(p.cls) + (p.region === 'outlier' ? RAMA_CLASSES.length : 0);
    return svg.append('g').selectAll('circle')
        .data([...points].sort((a, b) => order(a) - order(b)))
        .join('circle')
        .attr('class', p => `rama-point rama-${p.cls}${p.region === 'outlier' ? ' rama-outlier' : ''}`)
        .attr('cx', p => x(p.phi)).attr('cy', p => y(p.psi)).attr('r', p => p.region === 'outlier' ? 4.2 : 3.2);
}

/**
 * Loads the region grids per category once (static file
 * `data/rama8000_regions.json`, URL in `#angles-export[data-rama-regions]`).
 * @return {!Promise<!Object>} `{source, license, bins, regions: {type: string}}`.
 */
async function loadRamaRegions() {
    if (ramaView.regions) return ramaView.regions;
    const src = document.getElementById('angles-export')?.dataset.ramaRegions;
    ramaView.regions = await (await fetch(src)).json();
    return ramaView.regions;
}

/**
 * Renders plot and legend in the Ramachandran tab. Clicking a point selects the
 * residue (jumpToResidue). Stores `ramaAngles` and `ramaPoints` on the canvas so
 * that the plot can be redrawn when the category filter changes.
 * @param {!Object} canvas Canvas object.
 * @param {!Object} angles Result of analyzeAngles.
 */
function renderRamachandranPlot(canvas, angles) {
    const container = d3.select('#angle-rama');
    if (container.empty()) return;
    canvas.ramaAngles = angles;
    const shown = ramaView.type && ramaView.regions
        ? angles.ramachandran.filter(p => p.type === ramaView.type) : angles.ramachandran;
    container.selectAll('*').remove();
    const svg = container.append('svg')
        .attr('viewBox', `0 0 ${RAMA_SIZE} ${RAMA_SIZE}`)
        .attr('role', 'img')
        .attr('aria-label', interpolate(gettext('Ramachandran-Plot mit %s Aminosäuren'), [formatInt(shown.length)]));
    canvas.ramaPoints = drawRamachandran(svg, angles.ramachandran, RAMA_SIZE)
        .on('click', (event, p) => jumpToResidue(canvas, p.res));
    canvas.ramaPoints.append('title').text(p =>
        `${p.res.name} ${p.res.pos} (${chainLabel(p.res.chain)})\nφ = ${formatAngle(p.phi)}, ψ = ${formatAngle(p.psi)}`
        + `\n${ramaTypeLabel(p.type)} · ${ramaRegionLabel(p.region)}`);

    // Legend: classes (point colours), outliers and, if a category is selected, its regions
    const counts = d3.rollup(shown, v => v.length, p => p.cls);
    const entries = [
        ...RAMA_CLASSES.filter(cls => !ramaView.type || counts.get(cls)).map(cls => ({
            swatch: 'legend-dot', style: `background: var(--rama-${cls})`,
            label: `${ramaClassLabel(cls)} (${formatInt(counts.get(cls) || 0)})`,
        })),
        { swatch: 'legend-dot rama-legend-outlier',
          label: `${ramaRegionLabel('outlier')} (${formatInt(shown.filter(p => p.region === 'outlier').length)})` },
        ...(ramaView.type && ramaView.regions ? ['favored', 'allowed'].map(region => ({
            swatch: `rama-legend-region rama-legend-${region}`, label: ramaRegionLabel(region),
        })) : []),
    ];
    d3.select('#angle-rama-legend').selectAll('div')
        .data(entries)
        .join(enter => {
            // Align the symbol with the first text line (long labels wrap in the sidebar)
            const div = enter.append('div').attr('class', 'd-flex align-items-start gap-2');
            div.append('span');
            div.append('span').attr('class', 'rama-legend-label');
            return div;
        })
        .call(div => div.select('span').attr('class', e => `${e.swatch} rama-legend-dot`).attr('style', e => e.style || null))
        .call(div => div.select('.rama-legend-label').text(e => e.label));
    updateSelection(canvas);
}

/**
 * Ramachandran evaluation per category, totals and outliers (for the tab and the export).
 * @param {!Array<!Object>} points Ramachandran points.
 * @return {{total: number, favored: number, allowed: number, outlier: number,
 *     byType: !Array<!Object>, outliers: !Array<!Object>}} Counts per region,
 *     per category (only non-empty ones) and the outlier points.
 */
function ramachandranSummary(points) {
    const countRegions = list => Object.fromEntries(RAMA_REGIONS.map(r => [r, list.filter(p => p.region === r).length]));
    return {
        total: points.length,
        ...countRegions(points),
        byType: RAMA_TYPES.map(type => {
            const list = points.filter(p => p.type === type);
            return { type, count: list.length, ...countRegions(list) };
        }).filter(t => t.count),
        outliers: points.filter(p => p.region === 'outlier'),
    };
}

/**
 * Renders the Ramachandran key figures, the table per category and the outliers.
 * @param {!Object} angles Result of analyzeAngles.
 */
function renderRamachandranStats(angles) {
    const facts = d3.select('#rama-facts');
    if (facts.empty()) return;
    const s = ramachandranSummary(angles.ramachandran);
    const share = count => formatPercent(s.total ? 100 * count / s.total : 0);
    facts.selectAll('div')
        .data([
            [gettext('Bewertete Aminosäuren'), formatInt(s.total)],
            ...RAMA_REGIONS.map(r => [ramaRegionLabel(r), `${formatInt(s[r])} (${share(s[r])})`]),
        ])
        .join(enter => {
            const div = enter.append('div');
            div.append('dt');
            div.append('dd');
            return div;
        })
        .call(div => div.select('dt').text(d => d[0]))
        .call(div => div.select('dd').text(d => d[1]));

    const table = d3.select('#rama-types');
    table.selectAll('*').remove();
    table.append('thead').append('tr').selectAll('th')
        .data([[gettext('Kategorie')], [gettext('Anzahl'), true], ...RAMA_REGIONS.map(r => [ramaRegionLabel(r), true])])
        .join('th').attr('scope', 'col').attr('class', c => c[1] ? 'text-end' : null).text(c => c[0]);
    const row = (tr, label, t) => {
        tr.append('th').attr('scope', 'row').text(label);
        tr.append('td').attr('class', 'text-end').text(formatInt(t.count ?? t.total));
        RAMA_REGIONS.forEach(r => tr.append('td')
            .attr('class', `text-end${r === 'outlier' && t[r] ? ' z-outlier' : ''}`).text(formatInt(t[r])));
    };
    const body = table.append('tbody');
    s.byType.forEach(t => row(body.append('tr'), ramaTypeLabel(t.type), t));
    row(table.append('tfoot').append('tr').attr('class', 'share-total'), gettext('Gesamt'), s);

    const list = d3.select('#rama-outliers');
    list.selectAll('*').remove();
    d3.select('#rama-outliers-block').property('hidden', !s.outliers.length);
    if (!s.outliers.length) return;
    list.append('thead').append('tr').selectAll('th')
        .data([[gettext('Kette')], [gettext('Pos.')], [gettext('Aminosäure')], [gettext('Kategorie')], ['φ', true], ['ψ', true]])
        .join('th').attr('scope', 'col').attr('class', c => c[1] ? 'text-end' : null).text(c => c[0]);
    const rows = list.append('tbody').selectAll('tr').data(s.outliers).join('tr');
    rows.append('td').text(p => chainId(p.res.chain));
    rows.append('td').text(p => p.res.pos);
    rows.append('td').text(p => p.res.name);
    rows.append('td').text(p => ramaTypeLabel(p.type));
    rows.append('td').attr('class', 'text-end').text(p => formatAngle(p.phi));
    rows.append('td').attr('class', 'text-end').text(p => formatAngle(p.psi));
}

/**
 * Renders the whole Ramachandran tab: plot, Ramachandran statistics, ω
 * (cis/trans/twisted), bond-angle references and bond-angle groups.
 * @param {!Object} canvas Canvas object.
 * @param {!Object} graphData Graph data.
 */
function renderAngles(canvas, graphData) {
    if (d3.select('#angle-rama').empty()) return;
    const angles = analyzeAngles(graphData);
    renderRamachandranPlot(canvas, angles);
    renderRamachandranStats(angles);

    // ω
    const kinds = d3.rollup(angles.omegas, v => v.length, o => o.kind);
    const cis = angles.omegas.filter(o => o.kind === 'cis');
    const cisPro = cis.filter(o => ramaClass(o.res) === 'pro').length;
    d3.select('#angle-omega-facts').selectAll('div')
        .data([
            [gettext('Peptidbindungen'), formatInt(angles.omegas.length)],
            ['trans', formatInt(kinds.get('trans') || 0)],
            ['cis', interpolate(gettext('%(count)s (davon vor Prolin: %(pro)s)'), { count: formatInt(cis.length), pro: formatInt(cisPro) }, true)],
            [gettext('verdreht'), formatInt(kinds.get('twisted') || 0)],
        ])
        .join(enter => {
            const div = enter.append('div');
            div.append('dt');
            div.append('dd');
            return div;
        })
        .call(div => div.select('dt').text(d => d[0]))
        .call(div => div.select('dd').text(d => d[1]));

    const omegaTable = d3.select('#angle-omega');
    omegaTable.selectAll('*').remove();
    // List only if there are peptide bonds that are not trans; cis bonds not
    // preceding a proline are highlighted as unusual.
    const unusual = angles.omegas.filter(o => o.kind !== 'trans');
    d3.select('#angle-omega-block').property('hidden', !unusual.length);
    if (unusual.length) {
        omegaTable.append('thead').append('tr').selectAll('th')
            .data([[gettext('Kette')], [gettext('Pos.')], [gettext('Aminosäuren')], ['ω', true], [gettext('Art')]])
            .join('th').attr('scope', 'col').attr('class', c => c[1] ? 'text-end' : null).text(c => c[0]);
        const body = omegaTable.append('tbody');
        unusual.forEach(o => {
            const tr = body.append('tr');
            tr.append('td').text(chainId(o.res.chain));
            tr.append('td').text(`${o.prev.pos}–${o.res.pos}`);
            tr.append('td').text(`${o.prev.name}–${o.res.name}`);
            tr.append('td').attr('class', 'text-end').text(formatAngle(o.omega));
            tr.append('td').attr('class', o.kind === 'cis' && ramaClass(o.res) !== 'pro' ? 'z-outlier' : null)
                .text(omegaKindLabel(o));
        });
    }

    renderReferenceStats('angle-ref', rateReferences(angles.bondAngles), 'angle', parentsByName(graphData.residues));
    renderAngleGroups(angles);
    updateSelection(canvas);
}

const formatDegrees = v => v === null || v === undefined ? '–'
    : roundAngle(v).toLocaleString(LOCALE, { minimumFractionDigits: 1, maximumFractionDigits: 1 });
const formatDegreesDelta = v => v.toLocaleString(LOCALE,
    { minimumFractionDigits: 1, maximumFractionDigits: 1, signDisplay: 'exceptZero' });

/**
 * Describes the atoms of a bond angle for tables (like describeBond, with three
 * atoms); repeated chains, positions and residues are listed only once.
 * @param {!Array<!Object>} atoms The three atoms (graph nodes).
 * @return {{chain: string, pos: string, serials: string, residues: string}} Texts.
 */
function describeAngle(atoms) {
    const pos = atom => `${atom.residue_seq}${atom.i_code || ''}`;
    const unique = values => [...new Set(values)];
    return {
        chain: unique(atoms.map(a => chainId(a.chain))).join('–'),
        pos: unique(atoms.map(a => `${a.chain}|${pos(a)}`)).map(v => v.split('|')[1]).join('–'),
        serials: atoms.map(a => a.serial ?? '?').join('–'),
        residues: unique(atoms.map(a => `${a.chain}|${pos(a)}|${a.residue}`)).map(v => v.split('|')[2]).join('–'),
    };
}

// Expanded angle types; kept across re-renders.
const expandedAngleGroups = new Set();
/**
 * Renders bond angles per type (`#angle-bond-types`); each row expands to the
 * individual angles (like the bond groups). The detail list is built lazily on
 * first expansion.
 * @param {!Object} angles Result of analyzeAngles.
 */
function renderAngleGroups(angles) {
    const table = d3.select('#angle-bond-types');
    if (table.empty()) return;
    table.selectAll('*').remove();
    const columns = [
        { label: gettext('Winkel') }, { label: gettext('Anzahl'), num: true }, { label: gettext('Durchschnitt (°)'), num: true },
        { label: gettext('Std.-Abw. (°)'), num: true }, { label: gettext('Min (°)'), num: true, cls: 'length-min' },
        { label: gettext('Max (°)'), num: true, cls: 'length-max' },
    ];
    table.append('thead').append('tr').selectAll('th').data(columns).join('th')
        .attr('scope', 'col')
        .attr('class', c => [c.num && 'text-end', c.cls].filter(Boolean).join(' ') || null)
        .text(c => c.label);

    angles.angleTypes.forEach((group, i) => {
        const body = table.append('tbody').attr('class', 'bond-group');
        const detailsId = `angle-group-${i}`;
        const summary = body.append('tr').attr('class', 'bond-summary');
        const toggle = summary.append('th').attr('scope', 'row')
            .append('button').attr('type', 'button').attr('class', 'bond-toggle')
            .attr('aria-controls', detailsId)
            .html(CHEVRON_SVG);
        toggle.append('span').text(group.type);
        summary.append('td').attr('class', 'text-end').text(formatInt(group.count));
        [[group.mean], [group.sd], [group.min, 'length-min'], [group.max, 'length-max']].forEach(([v, cls]) =>
            summary.append('td').attr('class', cls ? `text-end ${cls}` : 'text-end').text(formatDegrees(v)));

        const details = body.append('tr').attr('class', 'bond-details').attr('id', detailsId);
        const cell = details.append('td').attr('colspan', columns.length);
        const setExpanded = expanded => {
            toggle.attr('aria-expanded', expanded);
            details.property('hidden', !expanded);
            if (expanded) {
                expandedAngleGroups.add(group.type);
                if (cell.select('table').empty()) renderAngleList(cell, group);
            } else {
                expandedAngleGroups.delete(group.type);
            }
        };
        summary.on('click', () => setExpanded(details.property('hidden')));
        setExpanded(expandedAngleGroups.has(group.type));
    });

    const total = table.append('tfoot').append('tr').attr('class', 'share-total');
    total.append('th').attr('scope', 'row').text(gettext('Bindungswinkel gesamt'));
    total.append('td').attr('class', 'text-end').text(formatInt(angles.bondAngles.length));
    total.append('td').attr('colspan', 4);
}

/**
 * Individual angles of one type, ordered by the residue of the centre atom;
 * smallest and largest are marked.
 * @param {!d3.Selection} cell Table cell to render into.
 * @param {!Object} group Angle type from analyzeAngles().angleTypes.
 */
function renderAngleList(cell, group) {
    const list = cell.append('table').attr('class', 'table table-sm bond-list');
    list.append('thead').append('tr').selectAll('th')
        .data([gettext('Winkel (°)'), gettext('Abw. vom Ø (°)'), gettext('Kette'), gettext('Pos.'), gettext('Atom-Nr.'), gettext('Aminosäure'), ''])
        .join('th').attr('scope', 'col')
        .attr('class', (c, i) => i < 2 ? 'text-end' : null)
        .text(c => c);
    const body = list.append('tbody');
    const sorted = [...group.items].sort((a, b) => d3.ascending(a.atoms[1].residue_index, b.atoms[1].residue_index));
    sorted.forEach(angle => {
        const d = describeAngle(angle.atoms);
        const tr = body.append('tr');
        tr.append('td').attr('class', 'text-end').text(formatDegrees(angle.value));
        tr.append('td').attr('class', 'text-end').text(formatDegreesDelta(angle.value - group.mean));
        [d.chain, d.pos, d.serials, d.residues].forEach(v => tr.append('td').text(v));
        const mark = tr.append('td');
        if (group.min === group.max) return;
        if (angle.value === group.min) mark.append('span').attr('class', 'badge length-min').text(gettext('Min'));
        if (angle.value === group.max) mark.append('span').attr('class', 'badge length-max').text(gettext('Max'));
    });
}

const omegaKindLabel = o => o.kind === 'cis'
    ? (ramaClass(o.res) === 'pro' ? gettext('cis (vor Prolin)') : gettext('cis'))
    : gettext('verdreht');

/**
 * Angle data for the JSON/PDF export (the `angles` part of buildStatsExport).
 * @param {!Object} graphData Graph data.
 * @return {!Object} Ramachandran points and evaluation, ω summary, bond angles
 *     per type and their reference rating.
 */
function anglesExport(graphData) {
    const angles = analyzeAngles(graphData);
    const round1 = v => v === null ? null : roundTo(v, 1);
    const residue = r => ({ residue: r.name, position: r.pos, chain: r.chain.trim() });
    return {
        unit: '°',
        ramachandran: angles.ramachandran.map(p => ({
            ...residue(p.res), class: p.cls, phi: round1(p.phi), psi: round1(p.psi),
            type: p.type, value: p.res.rama.value, region: p.region,
        })),
        // Regions after Top8000 as in MolProbity; type: general, glycine, cis-proline, trans-proline, pre-proline, ile-val
        ramachandran_evaluation: (s => ({
            source: structureInfo()?.ramachandran_source ?? null,
            total: s.total, favored: s.favored, allowed: s.allowed, outliers: s.outlier,
            by_type: s.byType.map(t => ({ type: t.type, count: t.count, favored: t.favored, allowed: t.allowed, outliers: t.outlier })),
        }))(ramachandranSummary(angles.ramachandran)),
        omega: {
            cis_below: CIS_MAX_OMEGA, trans_from: TRANS_MIN_OMEGA,
            count: angles.omegas.length,
            trans: angles.omegas.filter(o => o.kind === 'trans').length,
            cis: angles.omegas.filter(o => o.kind === 'cis').length,
            cis_before_proline: angles.omegas.filter(o => o.kind === 'cis' && ramaClass(o.res) === 'pro').length,
            twisted: angles.omegas.filter(o => o.kind === 'twisted').length,
            non_trans: angles.omegas.filter(o => o.kind !== 'trans').map(o => ({
                from: residue(o.prev), to: residue(o.res), omega: round1(o.omega), kind: o.kind,
            })),
        },
        bond_angles: angles.angleTypes.map(t => ({
            angle: t.type, count: t.count, mean: round1(t.mean), sd: round1(t.sd), min: round1(t.min), max: round1(t.max),
        })),
        // Deviation from reference values, Z = (angle − reference) / standard deviation
        reference: referenceExport(rateReferences(angles.bondAngles), '°', 1),
    };
}

// Ramachandran plot as a standalone SVG file (fixed colours, white background).
const RAMA_EXPORT_SIZE = { width: 800, height: 800 };
// A function because HIST_EXPORT_STYLE is defined further down in the file
// (accessing the const at load time would throw).
const ramaExportStyle = () => HIST_EXPORT_STYLE + `
    .rama-frame { stroke: #adb5bd; }
    .rama-minor line { stroke: rgba(0, 0, 0, 0.045); }
    .rama-zero { stroke: #adb5bd; stroke-dasharray: 3 3; }
    .rama-point { stroke: #fff; stroke-width: 0.8; }
    .rama-general { fill: #0d6efd; }
    .rama-gly { fill: #f08c00; }
    .rama-pro { fill: #2f9e44; }
    .rama-outlier { stroke: #e03131; stroke-width: 1.6; }
    .rama-region-allowed { fill: #e7f1fb; }
    .rama-region-favored { fill: #c5dcf3; }
`;
/**
 * Serialises the Ramachandran plot (current category filter applies) as SVG.
 * @param {!Object} graphData Graph data.
 * @return {string} SVG markup.
 */
function ramachandranSvgText(graphData) {
    const { width, height } = RAMA_EXPORT_SIZE;
    const svg = d3.create('svg')
        .attr('width', width).attr('height', height).attr('viewBox', `0 0 ${width} ${height}`);
    svg.append('style').text(ramaExportStyle());
    svg.append('rect').attr('width', width).attr('height', height).attr('fill', '#fff');
    drawRamachandran(svg, analyzeAngles(graphData).ramachandran, width);
    return new XMLSerializer().serializeToString(svg.node());
}

/**
 * Ramachandran plot as a downloadable image.
 * @param {!Object} graphData Graph data.
 * @param {string} format 'svg', 'png' or 'jpeg'.
 * @return {!Promise<!Blob>} The image.
 */
async function ramachandranBlob(graphData, format) {
    const svgText = ramachandranSvgText(graphData);
    if (format === 'svg') return new Blob([svgText], { type: 'image/svg+xml' });
    return rasterizeSvg(svgText, `image/${format}`, RAMA_EXPORT_SIZE);
}

/**
 * Draws the Ramachandran plot as vector graphics into the PDF (axes, grid,
 * points per class, outliers circled).
 * @param {!Object} doc jsPDF document.
 * @param {!Array<!Object>} points Exported Ramachandran points (`class`, `phi`,
 *     `psi`, `region`).
 * @param {number} left Left edge (mm).
 * @param {number} top Top edge (mm).
 * @param {number} size Width and height (mm).
 * @param {!Object<string, !Array<number>>} colors RGB colours: MUTED, LINE,
 *     OUTLIER and one per Ramachandran class.
 */
function drawPdfRamachandran(doc, points, left, top, size, colors) {
    const plotL = left + 12, plotR = left + size, plotT = top + 2, plotB = top + size - 10;
    const x = d3.scaleLinear().domain([-180, 180]).range([plotL, plotR]);
    const y = d3.scaleLinear().domain([-180, 180]).range([plotB, plotT]);
    doc.setFont('helvetica', 'normal').setFontSize(7).setTextColor(...colors.MUTED);
    doc.setDrawColor(242, 243, 245).setLineWidth(0.1);
    d3.range(-180, 181, RAMA_MINOR_STEP).filter(t => t % RAMA_MAJOR_STEP !== 0).forEach(t => {
        doc.line(x(t), plotT, x(t), plotB);
        doc.line(plotL, y(t), plotR, y(t));
    });
    doc.setDrawColor(...colors.LINE).setLineWidth(0.15);
    d3.range(-180, 181, RAMA_MAJOR_STEP).forEach(t => {
        doc.line(x(t), plotT, x(t), plotB);
        doc.line(plotL, y(t), plotR, y(t));
        doc.text(`${t}°`, x(t), plotB + 4, { align: 'center' });
        doc.text(`${t}°`, plotL - 1.5, y(t) + 1, { align: 'right' });
    });
    doc.setDrawColor(...colors.MUTED).setLineWidth(0.25).rect(plotL, plotT, plotR - plotL, plotB - plotT);
    const order = cls => RAMA_CLASSES.indexOf(cls);
    [...points].sort((a, b) => order(a.class) - order(b.class)).forEach(p => {
        doc.setFillColor(...colors[p.class]);
        doc.circle(x(p.phi), y(p.psi), 0.7, 'F');
    });
    doc.setDrawColor(...colors.OUTLIER).setLineWidth(0.3);
    points.filter(p => p.region === 'outlier').forEach(p => doc.circle(x(p.phi), y(p.psi), 1.3, 'S'));
    doc.setFontSize(8);
    doc.text(pdfText('φ'), (plotL + plotR) / 2, top + size - 1, { align: 'center' });
    doc.text(pdfText('ψ'), left + 1, (plotT + plotB) / 2, { angle: 90 });
}

// ===== Statistics: bonds =====
const formatLength = v => v === null || v === undefined ? '–'
    : v.toLocaleString(LOCALE, { minimumFractionDigits: 3, maximumFractionDigits: 3 });
const formatDelta = v => v.toLocaleString(LOCALE,
    { minimumFractionDigits: 3, maximumFractionDigits: 3, signDisplay: 'exceptZero' });
const chainId = c => c.trim() || '–';
const bondKey = (source, target) => `${source}|${target}`;

/**
 * Describes a bond by its two atoms (graph nodes or atoms from the ensemble). If
 * both atoms belong to the same residue, position and residue are named only once.
 * @param {!Array<!Object>} atoms The two atoms.
 * @return {{atoms: string, chain: string, pos: string, serials: string,
 *     residues: string}} Texts for the table columns.
 */
function describeBond([a, b]) {
    const pos = atom => `${atom.residue_seq}${atom.i_code || ''}`;
    const same = a.chain === b.chain && pos(a) === pos(b);
    return {
        atoms: `${a.label}–${b.label}`,
        chain: a.chain === b.chain ? chainId(a.chain) : `${chainId(a.chain)}–${chainId(b.chain)}`,
        pos: same ? pos(a) : `${pos(a)}–${pos(b)}`,
        serials: `${a.serial ?? '?'}–${b.serial ?? '?'}`,
        residues: same ? a.residue : `${a.residue}–${b.residue}`,
    };
}

/**
 * Appends a comparison value and the bond's deviation from it (length − value).
 * @param {!d3.Selection} td Target cell.
 * @param {?number|undefined} value Comparison value (Å), e.g. ensemble mean.
 * @param {number} length Bond length (Å).
 * @param {string} title Tooltip explaining the deviation.
 */
function appendComparison(td, value, length, title) {
    td.append('span').text(formatLength(value));
    if (value === null || value === undefined) return;
    td.append('span').attr('class', 'delta')
        .attr('title', title)
        .text(`(${formatDelta(length - value)})`);
}

// Model label of an extreme bond from the ensemble; files without MODEL records
// contain exactly one (implicit) model, reported as model 1.
const ensembleModelLabel = b => b && interpolate(gettext('Model %s'), [b.model ?? 1]);

/**
 * Renders the shortest and longest bond over all models and in the selected
 * model, each compared with the mean and median of that bond over all models.
 * @param {!Object} graphData Graph data of the selected model.
 * @param {!Object} ensemble Bond ensemble (`bonds`, `min`, `max`, `models`).
 * @param {string} modelLabel Label of the selected model.
 */
function renderBondExtremes(graphData, ensemble, modelLabel) {
    const table = d3.select('#stat-bond-extremes');
    table.selectAll('*').remove();

    const nodeById = new Map(graphData.nodes.map(n => [n.id, n]));
    const fromEdge = e => e && { ...e, atoms: [nodeById.get(e.source), nodeById.get(e.target)] };
    const groups = [
        { title: gettext('Alle Models'), rows: [
            { label: gettext('Kürzeste Bindung'), bond: ensemble.min, model: ensembleModelLabel(ensemble.min) },
            { label: gettext('Längste Bindung'),  bond: ensemble.max, model: ensembleModelLabel(ensemble.max) },
        ] },
        { title: gettext('Ausgewähltes Model'), rows: [
            { label: gettext('Kürzeste Bindung'), bond: fromEdge(d3.least(graphData.edges, e => e.length)), model: modelLabel },
            { label: gettext('Längste Bindung'),  bond: fromEdge(d3.greatest(graphData.edges, e => e.length)), model: modelLabel },
        ] },
    ];
    const columns = [
        { label: '' }, { label: gettext('Model') }, { label: gettext('Länge (Å)'), num: true }, { label: gettext('Atome') },
        { label: gettext('Kette') }, { label: gettext('Pos.') }, { label: gettext('Atom-Nr.') }, { label: gettext('Aminosäure') },
        { label: gettext('Durchschnitt (Å)'), num: true }, { label: gettext('Median (Å)'), num: true },
    ];

    table.append('thead').append('tr').selectAll('th').data(columns).join('th')
        .attr('scope', 'col')
        .attr('class', c => c.num ? 'text-end' : null)
        .text(c => c.label);

    groups.forEach(group => {
        const body = table.append('tbody');
        body.append('tr').attr('class', 'group-head')
            .append('th').attr('scope', 'rowgroup').attr('colspan', columns.length).text(group.title);
        group.rows.forEach(row => {
            const tr = body.append('tr');
            tr.append('th').attr('scope', 'row').text(row.label);
            if (!row.bond) {
                tr.append('td').attr('colspan', columns.length - 1).attr('class', 'text-body-secondary')
                    .text(gettext('Keine Bindungen'));
                return;
            }
            const d = describeBond(row.bond.atoms);
            const [mean, median] = ensemble.bonds[bondKey(row.bond.source, row.bond.target)] || [];
            tr.append('td').text(row.model);
            tr.append('td').attr('class', 'text-end fw-medium').text(formatLength(row.bond.length));
            [d.atoms, d.chain, d.pos, d.serials, d.residues].forEach(v => tr.append('td').text(v));
            appendComparison(tr.append('td').attr('class', 'text-end'), mean, row.bond.length,
                gettext('Abweichung der Bindung vom Durchschnitt'));
            appendComparison(tr.append('td').attr('class', 'text-end'), median, row.bond.length,
                gettext('Abweichung der Bindung vom Median'));
        });
    });
}

// Expanded groups stay open when the model changes.
const expandedBondGroups = new Set();
// Toggle icon of expandable table rows (rotated via CSS when expanded).
const CHEVRON_SVG = '<svg xmlns="http://www.w3.org/2000/svg" width="12" height="12" fill="currentColor" viewBox="0 0 16 16" aria-hidden="true"><path fill-rule="evenodd" d="M4.646 1.646a.5.5 0 0 1 .708 0l6 6a.5.5 0 0 1 0 .708l-6 6a.5.5 0 0 1-.708-.708L10.293 8 4.646 2.354a.5.5 0 0 1 0-.708"/></svg>';

/**
 * Key figures of a set of bonds.
 * @param {!Array<!Object>} bonds Bonds with `length` and `atoms`.
 * @return {{chains: string, count: number, mean: number, median: number,
 *     min: number, max: number, delta: number}} Summary (lengths in Å).
 */
function bondSummary(bonds) {
    const min = d3.min(bonds, b => b.length);
    const max = d3.max(bonds, b => b.length);
    return {
        chains: [...new Set(bonds.map(b => chainId(b.atoms[0].chain)))].join(', '),
        count: bonds.length,
        mean: d3.mean(bonds, b => b.length),
        median: d3.median(bonds, b => b.length),
        min, max, delta: max - min,
    };
}

/**
 * Renders all bonds of the selected model grouped by the atoms involved (e.g. C–N,
 * CB–CG). An expanded group lists every single bond (built lazily on first
 * expansion).
 * @param {!Object} graphData Graph data.
 */
function renderBondGroups(graphData) {
    const table = d3.select('#stat-bond-groups');
    table.selectAll('*').remove();

    const nodeById = new Map(graphData.nodes.map(n => [n.id, n]));
    const bonds = graphData.edges
        .map(e => ({ ...e, atoms: [nodeById.get(e.source), nodeById.get(e.target)] }))
        .filter(b => b.atoms[0] && b.atoms[1]);
    const groups = d3.groups(bonds, b => `${b.atoms[0].label}–${b.atoms[1].label}`)
        .map(([key, items]) => ({ key, items, ...bondSummary(items) }))
        .sort((a, b) => d3.descending(a.count, b.count) || d3.ascending(a.key, b.key));

    const columns = [
        { label: gettext('Atome') }, { label: gettext('Ketten') },
        { label: gettext('Anzahl'), num: true }, { label: gettext('Durchschnitt (Å)'), num: true },
        { label: gettext('Median (Å)'), num: true }, { label: gettext('Min (Å)'), num: true, cls: 'length-min' },
        { label: gettext('Max (Å)'), num: true, cls: 'length-max' }, { label: gettext('Δ Max − Min (Å)'), num: true },
    ];
    table.append('thead').append('tr').selectAll('th').data(columns).join('th')
        .attr('scope', 'col')
        .attr('class', c => [c.num && 'text-end', c.cls].filter(Boolean).join(' ') || null)
        .text(c => c.label);

    const summaryCells = (tr, s) => {
        tr.append('td').text(s.chains);
        tr.append('td').attr('class', 'text-end').text(formatInt(s.count));
        [[s.mean], [s.median], [s.min, 'length-min'], [s.max, 'length-max'], [s.delta]].forEach(([v, cls]) =>
            tr.append('td').attr('class', cls ? `text-end ${cls}` : 'text-end').text(formatLength(v)));
    };

    groups.forEach((group, i) => {
        const body = table.append('tbody').attr('class', 'bond-group');
        const detailsId = `bond-group-${i}`;
        const summary = body.append('tr').attr('class', 'bond-summary');
        const toggle = summary.append('th').attr('scope', 'row')
            .append('button').attr('type', 'button').attr('class', 'bond-toggle')
            .attr('aria-controls', detailsId)
            .html(CHEVRON_SVG);
        toggle.append('span').text(group.key);
        summaryCells(summary, group);

        const details = body.append('tr').attr('class', 'bond-details').attr('id', detailsId);
        const cell = details.append('td').attr('colspan', columns.length);

        const setExpanded = expanded => {
            toggle.attr('aria-expanded', expanded);
            details.property('hidden', !expanded);
            if (expanded) {
                expandedBondGroups.add(group.key);
                if (cell.select('table').empty()) renderBondList(cell, group);
            } else {
                expandedBondGroups.delete(group.key);
            }
        };
        summary.on('click', () => setExpanded(details.property('hidden')));
        setExpanded(expandedBondGroups.has(group.key));
    });

    const total = table.append('tfoot').append('tr').attr('class', 'share-total');
    total.append('th').attr('scope', 'row').text(gettext('Bindungen gesamt'));
    if (bonds.length) summaryCells(total, bondSummary(bonds));
}

/**
 * Individual bonds of a group in file order; shortest and longest are marked.
 * @param {!d3.Selection} cell Table cell to render into.
 * @param {!Object} group Group with `items`, `mean`, `min`, `max`.
 */
function renderBondList(cell, group) {
    const list = cell.append('table').attr('class', 'table table-sm bond-list');
    list.append('thead').append('tr').selectAll('th')
        .data([gettext('Länge (Å)'), gettext('Abw. vom Ø (Å)'), gettext('Kette'), gettext('Pos.'), gettext('Atom-Nr.'), gettext('Aminosäure'), ''])
        .join('th').attr('scope', 'col')
        .attr('class', (c, i) => i < 2 ? 'text-end' : null)
        .text(c => c);
    const body = list.append('tbody');
    group.items.forEach(bond => {
        const d = describeBond(bond.atoms);
        const tr = body.append('tr');
        tr.append('td').attr('class', 'text-end').text(formatLength(bond.length));
        tr.append('td').attr('class', 'text-end').text(formatDelta(bond.length - group.mean));
        [d.chain, d.pos, d.serials, d.residues].forEach(v => tr.append('td').text(v));
        const mark = tr.append('td');
        if (group.min === group.max) return;
        if (bond.length === group.min) mark.append('span').attr('class', 'badge length-min').text(gettext('Min'));
        if (bond.length === group.max) mark.append('span').attr('class', 'badge length-max').text(gettext('Max'));
    });
}

// ===== Histogram of bond lengths =====
const HIST_HEIGHT = 260;
const HIST_MARGIN = { top: 16, right: 16, bottom: 44, left: 64 };
// Last drawn bonds and container width; used to redraw on resize / when the tab is shown.
const histogram = { edges: [], width: 0 };
const HIST_BIN_WIDTH = 0.01; // Å

/**
 * Length bins of 0.01 Å with their share of all bonds; empty bins within the
 * value range are kept. Lengths have three decimals, so binning is done with
 * integers in hundredths of an Å (avoids floating-point bin errors).
 * @param {!Array<{length: number}>} edges Bonds.
 * @return {!Array<{i: number, x0: number, x1: number, count: number,
 *     share: number}>} Bins; `i` = bin index in 0.01 Å, `share` in percent.
 */
function histogramBins(edges) {
    if (!edges.length) return [];
    const binOf = len => Math.floor(Math.round(len * 1000) / 10);
    const counts = d3.rollup(edges, v => v.length, e => binOf(e.length));
    const [lo, hi] = d3.extent(counts.keys());
    return d3.range(lo, hi + 1).map(i => {
        const count = counts.get(i) || 0;
        return { i, x0: i / 100, x1: (i + 1) / 100, count, share: 100 * count / edges.length };
    });
}

const formatAxisLength = v => v.toLocaleString(LOCALE, { minimumFractionDigits: 2, maximumFractionDigits: 2 });

/**
 * SVG path of a bar with rounded top corners, flush with the baseline at the bottom.
 * @param {number} x Left edge.
 * @param {number} y Top edge.
 * @param {number} w Width.
 * @param {number} h Height.
 * @param {number} r Corner radius (clamped to the bar size).
 * @return {string} Path data.
 */
function topRoundedBar(x, y, w, h, r) {
    r = Math.min(r, w / 2, h);
    return `M${x},${y + h}V${y + r}Q${x},${y} ${x + r},${y}H${x + w - r}Q${x + w},${y} ${x + w},${y + r}V${y + h}Z`;
}

/**
 * Renders the histogram of the share of bonds per length bin (0.01 Å) among all
 * bonds of the selected model, with a hover tooltip per bin. Without an argument
 * it redraws with the bonds passed last (e.g. after a resize).
 * @param {!Array<!Object>=} edges Bonds to show.
 */
function renderBondHistogram(edges = histogram.edges) {
    histogram.edges = edges;
    const container = d3.select('#stat-bond-histogram');
    if (container.empty()) return;
    const width = container.node().clientWidth;
    // Tab hidden: reset the width so that the ResizeObserver redraws once it is shown.
    if (!width) { histogram.width = 0; return; }
    histogram.width = width;
    container.selectAll('*').remove();
    if (!edges.length) return;

    const bins = histogramBins(edges);
    const m = HIST_MARGIN;
    const svg = container.append('svg')
        .attr('width', width).attr('height', HIST_HEIGHT)
        .attr('role', 'img')
        .attr('aria-label', interpolate(
            gettext('Histogramm: Anteil der %s Bindungen je Bindungslänge, Klassenbreite 0,01 Å'), [formatInt(edges.length)]));
    const { x, y, bars } = drawHistogram(svg, bins, width, HIST_HEIGHT);

    // Hover: the hit area is the whole column, not just the bar.
    const tip = container.append('div').attr('class', 'hist-tooltip').property('hidden', true);
    svg.append('g').selectAll('rect')
        .data(bins)
        .join('rect')
        .attr('x', b => x(b.x0)).attr('width', b => x(b.x1) - x(b.x0))
        .attr('y', m.top).attr('height', y(0) - m.top)
        .attr('fill', 'transparent')
        .on('pointerenter', (event, b) => {
            bars.classed('dimmed', d => d !== b);
            tip.selectAll('*').remove();
            tip.append('div').attr('class', 'fw-medium')
                .text(`${formatLength(b.i / 100)} – ${formatLength((b.i * 10 + 9) / 1000)} Å`);
            tip.append('div').attr('class', 'text-body-secondary')
                .text(`${formatPercent(b.share)} · ${interpolate(
                    ngettext('%s Bindung', '%s Bindungen', b.count), [formatInt(b.count)])}`);
            const tipWidth = tip.property('hidden', false).node().offsetWidth;
            const cx = (x(b.x0) + x(b.x1)) / 2;
            tip.style('left', `${Math.min(Math.max(cx, tipWidth / 2), width - tipWidth / 2)}px`)
                .style('top', `${y(b.share)}px`);
        })
        .on('pointerleave', () => {
            bars.classed('dimmed', false);
            tip.property('hidden', true);
        });
}

/**
 * Draws axes, labels and bars of the histogram into `svg` (used by the tab and the
 * image export).
 * @param {!d3.Selection} svg Target SVG selection.
 * @param {!Array<!Object>} bins Result of histogramBins (non-empty).
 * @param {number} width Width (px).
 * @param {number} height Height (px).
 * @return {{x: !d3.ScaleLinear, y: !d3.ScaleLinear, bars: !d3.Selection}} Scales
 *     and bar selection (for the hover interaction).
 */
function drawHistogram(svg, bins, width, height) {
    const m = HIST_MARGIN;
    const x = d3.scaleLinear().domain([bins[0].x0, bins[bins.length - 1].x1]).range([m.left, width - m.right]);
    const y = d3.scaleLinear().domain([0, d3.max(bins, b => b.share)]).nice()
        .range([height - m.bottom, m.top]);

    const yAxis = svg.append('g').attr('class', 'hist-axis').attr('transform', `translate(${m.left},0)`)
        .call(d3.axisLeft(y).ticks(5).tickSize(-(width - m.left - m.right)).tickPadding(8)
            .tickFormat(formatPercentTick));
    yAxis.select('.domain').remove();
    yAxis.selectAll('.tick line').attr('class', 'hist-grid');

    svg.append('g').attr('class', 'hist-axis').attr('transform', `translate(0,${y(0)})`)
        .call(d3.axisBottom(x).ticks(Math.max(2, Math.floor((width - m.left - m.right) / 80)))
            .tickFormat(formatAxisLength));

    svg.append('text').attr('class', 'hist-label').attr('text-anchor', 'middle')
        .attr('x', (m.left + width - m.right) / 2).attr('y', height - 6)
        .text(gettext('Bindungslänge (Å)'));
    svg.append('text').attr('class', 'hist-label').attr('text-anchor', 'middle')
        .attr('transform', `translate(14,${(m.top + y(0)) / 2}) rotate(-90)`)
        .text(gettext('Anteil der Bindungen'));

    // Bars with a 2px gap; even very small shares stay visible as 1px bars.
    const bars = svg.append('g').selectAll('path')
        .data(bins.filter(b => b.count))
        .join('path')
        .attr('class', 'hist-bar')
        .attr('d', b => {
            const h = Math.max(1, y(0) - y(b.share));
            return topRoundedBar(x(b.x0) + 1, y(0) - h, Math.max(1, x(b.x1) - x(b.x0) - 2), h, 4);
        });
    return { x, y, bars };
}

// ===== Image export of the histogram (SVG / PNG / JPEG) =====
// Fixed light colours as in the PDF so that the image looks the same regardless of dark mode.
const HIST_EXPORT_SIZE = { width: 1200, height: 480 };
const HIST_EXPORT_SCALE = 2; // Pixel density of raster images
const HIST_EXPORT_STYLE = `
    text { font-family: Roboto, Arial, sans-serif; }
    .hist-bar { fill: #0d6efd; }
    .hist-axis text { fill: #6c757d; font-size: 13px; font-variant-numeric: tabular-nums; }
    .hist-axis .domain, .hist-axis .tick line { stroke: #dee2e6; }
    .hist-axis .tick line.hist-grid { stroke: rgba(0, 0, 0, 0.1); }
    .hist-label { fill: #6c757d; font-size: 14px; }
`;

/**
 * Standalone SVG file of the histogram (white background, embedded styles).
 * @param {!Array<!Object>} edges Bonds.
 * @return {string} SVG markup.
 */
function histogramSvgText(edges) {
    const { width, height } = HIST_EXPORT_SIZE;
    const svg = d3.create('svg')
        .attr('width', width).attr('height', height).attr('viewBox', `0 0 ${width} ${height}`);
    svg.append('style').text(HIST_EXPORT_STYLE);
    svg.append('rect').attr('width', width).attr('height', height).attr('fill', '#fff');
    drawHistogram(svg, histogramBins(edges), width, height);
    return new XMLSerializer().serializeToString(svg.node());
}

/**
 * Rasterises SVG markup via an <img> and a <canvas> at HIST_EXPORT_SCALE.
 * @param {string} svgText SVG markup.
 * @param {string} type 'image/png' or 'image/jpeg'.
 * @param {{width: number, height: number}=} size Image size in CSS px.
 * @return {!Promise<?Blob>} The raster image.
 */
async function rasterizeSvg(svgText, type, size = HIST_EXPORT_SIZE) {
    const { width, height } = size;
    const url = URL.createObjectURL(new Blob([svgText], { type: 'image/svg+xml' }));
    try {
        const img = new Image();
        await new Promise((resolve, reject) => {
            img.onload = resolve;
            img.onerror = () => reject(new Error('SVG konnte nicht gerendert werden'));
            img.src = url;
        });
        const canvas = document.createElement('canvas');
        canvas.width = width * HIST_EXPORT_SCALE;
        canvas.height = height * HIST_EXPORT_SCALE;
        canvas.getContext('2d').drawImage(img, 0, 0, canvas.width, canvas.height);
        return await new Promise(resolve => canvas.toBlob(resolve, type, 0.92));
    } finally {
        URL.revokeObjectURL(url);
    }
}

/**
 * Histogram as a downloadable image.
 * @param {!Array<!Object>} edges Bonds.
 * @param {string} format 'svg', 'png' or 'jpeg'.
 * @return {!Promise<!Blob>} The image.
 */
async function histogramBlob(edges, format) {
    const svgText = histogramSvgText(edges);
    if (format === 'svg') return new Blob([svgText], { type: 'image/svg+xml' });
    return rasterizeSvg(svgText, `image/${format}`);
}

// ===== Statistics: disulfide and peptide bonds =====
/**
 * Disulfide bridges with their atoms (graph nodes), in file order.
 * @param {!Object} graphData Graph data.
 * @return {!Array<!Object>} Edges extended by `atoms`; `via` tells whether the
 *     bridge comes from an SSBOND record ('ssbond') or the SG–SG distance.
 */
function disulfideBonds(graphData) {
    const nodeById = new Map(graphData.nodes.map(n => [n.id, n]));
    return graphData.edges
        .filter(e => e.type === 'disulfide' && nodeById.has(e.source) && nodeById.has(e.target))
        .map(e => ({ ...e, atoms: [nodeById.get(e.source), nodeById.get(e.target)] }));
}

/**
 * Partners of the disulfide bridges per cysteine (sequence-bar tooltips).
 * @param {!Object} graphData Graph data.
 * @return {!Map<string, !Array<string>>} residueKey -> partner labels, e.g.
 *     "CYS 127" or "CYS 214 (Kette L)" (the chain only if it differs).
 */
function disulfidePartners(graphData) {
    const partners = new Map();
    const label = (n, other) => `${n.residue} ${n.residue_seq}${n.i_code || ''}` +
        (n.chain !== other.chain ? ` (${interpolate(gettext('Kette %s'), [chainId(n.chain)])})` : '');
    disulfideBonds(graphData).forEach(({ atoms: [a, b] }) => {
        [[a, b], [b, a]].forEach(([self, other]) => {
            const key = nodeResidueKey(self);
            partners.set(key, [...(partners.get(key) || []), label(other, self)]);
        });
    });
    return partners;
}

const disulfideOrigin = via => via === 'ssbond' ? gettext('SSBOND-Eintrag') : gettext('Abstand');

/**
 * Renders the disulfide table (hidden if there are none).
 * @param {!Object} graphData Graph data.
 */
function renderDisulfides(graphData) {
    const table = d3.select('#stat-disulfides');
    if (table.empty()) return;
    table.selectAll('*').remove();
    const bonds = disulfideBonds(graphData);
    d3.select('#stat-disulfides-block').property('hidden', !bonds.length);
    if (!bonds.length) return;
    table.append('thead').append('tr').selectAll('th')
        .data([[gettext('Kette')], [gettext('Pos.')], [gettext('Länge (Å)'), true], [gettext('Herkunft')]])
        .join('th').attr('scope', 'col').attr('class', c => c[1] ? 'text-end' : null).text(c => c[0]);
    const body = table.append('tbody');
    bonds.forEach(b => {
        const d = describeBond(b.atoms);
        const tr = body.append('tr');
        tr.append('td').text(d.chain);
        tr.append('td').text(d.pos);
        tr.append('td').attr('class', 'text-end').text(formatLength(b.length));
        tr.append('td').text(disulfideOrigin(b.via));
    });
    const total = table.append('tfoot').append('tr').attr('class', 'share-total');
    total.append('th').attr('scope', 'row').attr('colspan', 2).text(gettext('Disulfidbrücken gesamt'));
    total.append('td').attr('class', 'text-end').text(formatInt(bonds.length));
    total.append('td');
}

/**
 * Renders peptide bonds: existing and missing (chain breaks), relative to the
 * theoretically possible ones.
 * @param {!Object} graphData Graph data.
 */
function renderPeptideBonds(graphData) {
    const { possible, formed, breaks } = analyzePeptideBonds(graphData);
    d3.select('#stat-peptide-possible').text(formatInt(possible));
    renderShareTable('#stat-peptide-bonds', [
        { key: 'formed', label: gettext('Existent'), count: formed, color: 'var(--bs-primary)' },
        { key: 'missing', label: gettext('Nicht-existent'), sub: gettext('Kettenlücken'),
          count: breaks.size, color: 'var(--bs-danger)' },
    ], possible || 1, { sort: false });
}

/**
 * Renders the bond part of the statistics tab (references, peptide bonds,
 * disulfides, extremes, histogram, groups).
 * @param {!Object} graphData Graph data of the selected model.
 * @param {?Object} ensemble Bond ensemble; nothing is rendered without it.
 * @param {string} modelLabel Label of the selected model.
 */
function renderBondStats(graphData, ensemble, modelLabel) {
    if (d3.select('#stat-bond-extremes').empty() || !ensemble) return;
    renderReferenceStats('stat-ref', referenceStats(graphData), 'bond', parentsByName(graphData.residues));
    renderPeptideBonds(graphData);
    renderDisulfides(graphData);
    renderBondExtremes(graphData, ensemble, modelLabel);
    renderBondHistogram(graphData.edges);
    renderBondGroups(graphData);
}

// ===== Export of the statistics (JSON / PDF) =====
// Rounds to `digits` decimals and returns a number (null stays null).
const roundTo = (v, digits) => v === null || v === undefined ? null : +v.toFixed(digits);
const sharePercent = (count, total) => roundTo(100 * count / total, 2);

// Atom of a bond in the export (graph nodes and ensemble atoms have the same fields).
const exportAtom = n => ({
    name: n.label, serial: n.serial, residue: n.residue,
    position: `${n.residue_seq}${n.i_code || ''}`, chain: n.chain.trim(),
});
// Reverse mapping for labels in the PDF (describeBond expects fields like the graph nodes).
const nodeLikeAtom = a => ({ label: a.name, serial: a.serial, residue: a.residue, residue_seq: a.position, chain: a.chain });

/**
 * Structure metadata from the PDB header (embedded by the server as
 * `#structure-info`): method, resolution, R values, reference sources.
 * @return {?Object} Parsed JSON, or null if the element is missing.
 */
const structureInfo = () => {
    const el = document.getElementById('structure-info');
    return el ? JSON.parse(el.textContent) : null;
};

/**
 * Collects all values of the statistics tab (plus the angles) into one plain
 * object; the basis for both the JSON and the PDF export. The keys are part of
 * the JSON export format, so keep them stable.
 * @param {!Object} graphData Graph data of the selected model.
 * @param {?Object} ensemble Bond ensemble.
 * @param {string} modelLabel Label of the selected model.
 * @param {string} filename Name of the uploaded PDB file.
 * @return {!Object} Export data.
 */
function buildStatsExport(graphData, ensemble, modelLabel, filename) {
    const { nodes, edges, residues } = graphData;
    const aminoAcids = residues.filter(r => r.one_letter !== '?');
    const oneLetter = new Map(aminoAcids.map(r => [r.name, r.one_letter]));

    const describeElement = el => ({ element: elementSymbol(el) });
    const parents = parentsByName(aminoAcids);
    const describeAminoAcid = code => ({
        code, one_letter: oneLetter.get(code), name: aminoAcidLabel(code, parents.get(code)),
        ...(parents.has(code) ? { parent: parents.get(code) } : {}),
    });
    const composition = (items, keyOf, describe) => d3.rollups(items, v => v.length, keyOf)
        .sort((a, b) => d3.descending(a[1], b[1]) || d3.ascending(a[0], b[0]))
        .map(([key, count]) => ({ ...describe(key), count, share_percent: sharePercent(count, items.length) }));
    const byChain = (items, keyOf, describe, countName) =>
        [...new Set(items.map(i => i.chain))].map(c => {
            const inChain = items.filter(i => i.chain === c);
            return {
                chain: c.trim(), [countName]: inChain.length,
                share_percent: sharePercent(inChain.length, items.length),
                composition: composition(inChain, keyOf, describe),
            };
        });

    const nodeById = new Map(nodes.map(n => [n.id, n]));
    const bonds = edges
        .map(e => ({ ...e, atoms: [nodeById.get(e.source), nodeById.get(e.target)] }))
        .filter(b => b.atoms[0] && b.atoms[1]);
    const exportBond = (bond, model) => {
        if (!bond) return null;
        const [mean, median] = ensemble?.bonds[bondKey(bond.source, bond.target)] || [];
        return {
            length: bond.length, model, atoms: bond.atoms.map(exportAtom),
            ensemble_mean: mean ?? null, ensemble_median: median ?? null,
        };
    };
    const summary = items => {
        const s = bondSummary(items);
        return {
            count: s.count, mean: roundTo(s.mean, 3), median: roundTo(s.median, 3),
            min: s.min, max: s.max, delta: roundTo(s.delta, 3),
        };
    };
    const chainsOf = items => [...new Set(items.map(b => b.atoms[0].chain.trim()))];
    const selectedMin = bonds.length ? bonds[d3.leastIndex(bonds, b => b.length)] : null;
    const selectedMax = bonds.length ? bonds[d3.greatestIndex(bonds, b => b.length)] : null;
    const peptide = analyzePeptideBonds(graphData);
    const reference = referenceStats(graphData);
    const info = structureInfo();
    const exportResidue = r => ({ residue: r.name, position: r.pos, chain: r.chain.trim() });

    return {
        file: filename,
        model: modelLabel,
        created: new Date().toISOString(),
        general: {
            amino_acids: aminoAcids.length,
            atoms: nodes.length,
            chains: new Set(residues.map(r => r.chain)).size,
            models: ensemble?.models?.length || 1,
            bonds: bonds.length,
            method: info?.method ?? null,
            resolution: info?.resolution ?? null,
            r_work: info?.r_work ?? null,
            r_free: info?.r_free ?? null,
            elements: composition(nodes, n => n.element || '', describeElement),
            amino_acid_composition: composition(aminoAcids, r => r.name, describeAminoAcid),
        },
        angles: anglesExport(graphData),
        elements_by_chain: byChain(nodes, n => n.element || '', describeElement, 'atoms'),
        amino_acids_by_chain: byChain(aminoAcids, r => r.name, describeAminoAcid, 'amino_acids'),
        bonds: {
            length_unit: 'Å',
            ...(bonds.length ? summary(bonds) : { count: 0 }),
            // Deviation from reference values, Z = (length − reference) / standard deviation
            reference: referenceExport(reference, 'Å', 3),
            // via: "ssbond" (SSBOND record) or "distance" (SG–SG distance)
            disulfides: disulfideBonds(graphData).map(b => ({ length: b.length, via: b.via, atoms: b.atoms.map(exportAtom) })),
            // Between consecutive amino acids of a chain; missing = chain breaks
            peptide_bonds: {
                possible: peptide.possible,
                formed: peptide.formed,
                missing: peptide.breaks.size,
                chain_breaks: [...peptide.breaks].map(([index, gap]) => ({
                    from: exportResidue(gap.prev), to: exportResidue(residues[index]),
                    distance: roundTo(gap.distance, 3),
                })),
            },
            extremes: {
                all_models: ensemble && {
                    shortest: exportBond(ensemble.min, ensembleModelLabel(ensemble.min)),
                    longest: exportBond(ensemble.max, ensembleModelLabel(ensemble.max)),
                },
                selected_model: {
                    shortest: exportBond(selectedMin, modelLabel),
                    longest: exportBond(selectedMax, modelLabel),
                },
            },
            // Bins are half-open: from <= length < to.
            histogram: {
                bin_width: HIST_BIN_WIDTH,
                bins: histogramBins(bonds).map(b => ({
                    from: roundTo(b.x0, 2), to: roundTo(b.x1, 2), count: b.count, share_percent: roundTo(b.share, 2),
                })),
            },
            by_atoms: d3.groups(bonds, b => `${b.atoms[0].label}–${b.atoms[1].label}`)
                .sort((a, b) => d3.descending(a[1].length, b[1].length) || d3.ascending(a[0], b[0]))
                .map(([, items]) => ({
                    atoms: [items[0].atoms[0].label, items[0].atoms[1].label],
                    chains: chainsOf(items),
                    ...summary(items),
                    bonds: items.map(b => ({ length: b.length, atoms: b.atoms.map(exportAtom) })),
                })),
        },
    };
}

// ===== CSV export =====
// With a decimal comma (German) a semicolon separates the columns, otherwise a
// comma, so that a spreadsheet application in the same language opens the file
// directly. A BOM marks the encoding as UTF-8.
const CSV_DECIMAL = (1.5).toLocaleString(LOCALE).includes(',') ? ',' : '.';
const CSV_SEPARATOR = CSV_DECIMAL === ',' ? ';' : ',';
const csvNumber = (v, digits) => v === null || v === undefined || Number.isNaN(v) ? ''
    : (roundTo(v, digits) || 0).toFixed(digits).replace('.', CSV_DECIMAL);
/**
 * Quotes a CSV field if it contains the separator, quotes or line breaks.
 * @param {*} v Value.
 * @return {string} The field.
 */
function csvField(v) {
    const s = v === null || v === undefined ? '' : String(v);
    return s.includes(CSV_SEPARATOR) || /["\r\n]/.test(s) ?`"${s.replaceAll('"', '""')}"` : s;
}
const csvText = (header, rows) =>
    '\ufeff' + [header, ...rows].map(r => r.map(csvField).join(CSV_SEPARATOR)).join('\r\n') + '\r\n';

// Columns of one atom in bond and angle tables; `i` numbers the atoms.
const csvAtomHeader = i => [gettext('Kette'), gettext('Pos.'), gettext('Aminosäure'), gettext('Atom'), gettext('Atom-Nr.')]
    .map(label => `${label} ${i}`);
const csvAtom = n => [n.chain.trim(), `${n.residue_seq}${n.i_code || ''}`, n.residue, n.label, n.serial];

/**
 * Builds one CSV export of the selected model.
 * @param {string} kind 'atoms', 'bonds', 'angles' or 'dihedrals' (`data-csv` of
 *     the export button).
 * @param {!Object} graphData Graph data.
 * @return {string} CSV text including BOM.
 */
function csvTable(kind, graphData) {
    const nodeById = new Map(graphData.nodes.map(n => [n.id, n]));
    if (kind === 'atoms') {
        return csvText(DATA_COLUMNS.map(c => c.label), graphData.nodes.map(n => DATA_COLUMNS.map(c => {
            const v = c.value(n);
            return c.num && typeof v === 'string' ? v.replace('.', CSV_DECIMAL) : typeof v === 'string' ? v.trim() : v;
        })));
    }
    if (kind === 'bonds') {
        const rows = graphData.edges
            .map(e => [e, nodeById.get(e.source), nodeById.get(e.target)])
            .filter(([, a, b]) => a && b)
            .map(([e, a, b]) => [
                ...csvAtom(a), ...csvAtom(b), BOND_TYPE_LABELS[e.type] || e.type, csvNumber(e.length, 3),
                csvNumber(e.ref?.[0], 3), csvNumber(e.ref?.[1], 3), csvNumber(zScore(e), 2),
            ]);
        return csvText([...csvAtomHeader(1), ...csvAtomHeader(2), gettext('Typ'), gettext('Länge (Å)'),
            gettext('Referenz (Å)'), gettext('σ (Å)'), 'Z'], rows);
    }
    const angles = analyzeAngles(graphData);
    if (kind === 'angles') {
        const rows = angles.bondAngles.map(a => [
            ...a.atoms.flatMap(csvAtom), a.type, csvNumber(a.value, 1),
            csvNumber(a.ref?.[0], 1), csvNumber(a.ref?.[1], 1), csvNumber(a.z, 2),
        ]);
        return csvText([...csvAtomHeader(1), ...csvAtomHeader(2), ...csvAtomHeader(3), gettext('Typ'),
            gettext('Winkel (°)'), gettext('Referenz (°)'), gettext('σ (°)'), 'Z'], rows);
    }
    // 'dihedrals': torsion angles per amino acid; ω belongs to the peptide bond before the amino acid.
    const omegaOf = new Map(angles.omegas.map(o => [residueKey(o.res), o]));
    const rows = graphData.residues
        .filter(res => res.rama || omegaOf.has(residueKey(res)))
        .map(res => {
            const omega = omegaOf.get(residueKey(res));
            return [
                res.chain.trim(), res.pos, res.name,
                csvNumber(res.rama?.phi, 1), csvNumber(res.rama?.psi, 1), csvNumber(omega?.omega, 1),
                res.rama ? ramaTypeLabel(res.rama.type) : '', res.rama ? ramaRegionLabel(res.rama.region) : '',
            ];
        });
    return csvText([gettext('Kette'), gettext('Pos.'), gettext('Aminosäure'), 'φ (°)', 'ψ (°)', 'ω (°)',
        gettext('Kategorie'), gettext('Bereich')], rows);
}

// Translated file name part per CSV kind.
const CSV_FILE_KINDS = {
    atoms: gettext('atome'), bonds: gettext('bindungen'), angles: gettext('bindungswinkel'), dihedrals: gettext('diederwinkel'),
};

// Download file name without extension: "<pdb name>_<kind>_<model slug>".
const exportBaseName = (filename, model, kind = gettext('statistiken')) =>
    `${filename.replace(/\.[^.]+$/, '')}_${kind}_${model.toLowerCase().replace(/[^\p{L}\p{N}]+/gu, '-')}`;

/**
 * Triggers a browser download of a blob via a temporary link.
 * @param {!Blob} blob File content.
 * @param {string} name File name.
 */
function downloadBlob(blob, name) {
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = name;
    document.body.append(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 0);
}

// Scripts (jsPDF) are loaded only on the first PDF export; src -> Promise.
const loadedScripts = new Map();
/**
 * Loads a script once; a failed load is forgotten so that it can be retried.
 * @param {string} src Script URL.
 * @return {!Promise} Resolves when the script has loaded.
 */
function loadScript(src) {
    if (!loadedScripts.has(src)) {
        loadedScripts.set(src, new Promise((resolve, reject) => {
            const script = document.createElement('script');
            script.src = src;
            script.onload = resolve;
            script.onerror = () => {
                loadedScripts.delete(src);
                reject(new Error(`${src} konnte nicht geladen werden`));
            };
            document.head.append(script);
        }));
    }
    return loadedScripts.get(src);
}

// The PDF standard fonts only know WinAnsi characters (no Greek letters, no Δ, no
// minus sign U+2212), so these are spelled out. At module level because
// drawPdfHistogram and drawPdfRamachandran need it as well.
const pdfText = v => String(v).replace(/φ/g, 'phi').replace(/ψ/g, 'psi').replace(/ω/g, 'omega').replace(/Δ/g, 'Delta').replace(/−/g, '-');

/**
 * Builds the PDF from the export data: the same sections and tables as in the
 * statistics and Ramachandran tabs, laid out top to bottom with a running `y`
 * cursor (mm). The list of individual bonds per atom pair is optional because it
 * runs to thousands of rows (many pages) for large proteins; the JSON export
 * always contains it. Requires jsPDF and jspdf-autotable to be loaded.
 * @param {!Object} data Result of buildStatsExport.
 * @param {{bondList: (boolean|undefined)}=} options `bondList`: include the
 *     individual bonds per atom pair.
 * @return {!Object} The jsPDF document.
 */
function buildStatsPdf(data, options = {}) {
    const doc = new window.jspdf.jsPDF({ unit: 'mm', format: 'a4' });
    const M = 15;
    const pageW = doc.internal.pageSize.getWidth();
    const pageH = doc.internal.pageSize.getHeight();
    const contentW = pageW - 2 * M;
    const INK = [33, 37, 41], MUTED = [108, 117, 125], LINE = [222, 226, 230], PRIMARY = [13, 110, 253];
    const LENGTH_MIN = [11, 114, 133], LENGTH_MAX = [217, 72, 15]; // like --length-min/max-color (light theme)
    let y = M;

    // Starts a new page if fewer than `h` mm are left.
    const ensureSpace = h => { if (y + h > pageH - M) { doc.addPage(); y = M; } };
    const sectionTitle = text => {
        ensureSpace(50);
        y += 2;
        doc.setFont('helvetica', 'bold').setFontSize(12).setTextColor(...INK).text(pdfText(text), M, y);
        doc.setDrawColor(...LINE).setLineWidth(0.3).line(M, y + 2, pageW - M, y + 2);
        y += 8;
    };
    // Subheading only where there is still room for the header and a few table rows.
    const subTitle = (text, minSpace = 40) => {
        ensureSpace(minSpace);
        doc.setFont('helvetica', 'bold').setFontSize(9.5).setTextColor(...INK).text(pdfText(text), M, y);
        y += 2.5;
    };
    // autoTable with the common style; all cell texts pass through pdfText. Advances `y`.
    const table = (options) => {
        doc.autoTable({
            startY: y,
            margin: { left: M, right: M },
            theme: 'plain',
            styles: { font: 'helvetica', fontSize: 8, cellPadding: { top: 1.2, bottom: 1.2, left: 1.5, right: 1.5 }, textColor: INK, lineColor: LINE },
            headStyles: { fontStyle: 'bold', lineWidth: { bottom: 0.3 } },
            bodyStyles: { lineWidth: { bottom: 0.1 } },
            footStyles: { fontStyle: 'bold', fillColor: [242, 243, 245], lineWidth: { top: 0.4 } },
            showHead: 'everyPage',
            showFoot: 'lastPage',
            horizontalPageBreak: true,
            horizontalPageBreakRepeat: 0,
            ...options,
            head: options.head.map(row => row.map(c => typeof c === 'object' ? { ...c, content: pdfText(c.content) } : pdfText(c))),
            body: options.body.map(row => row.map(c => typeof c === 'object' && c !== null ? { ...c, content: pdfText(c.content) } : pdfText(c ?? ''))),
            foot: options.foot && options.foot.map(row => row.map(c => typeof c === 'object' ? { ...c, content: pdfText(c.content) } : pdfText(c))),
        });
        y = doc.lastAutoTable.finalY + 7;
    };
    const right = cols => Object.fromEntries(cols.map(i => [i, { halign: 'right' }]));
    const num = (content) => ({ content, styles: { halign: 'right' } });
    const parents = new Map(data.general.amino_acid_composition.filter(a => a.parent).map(a => [a.code, a.parent]));

    // Deviation from reference values (bonds or angles): per section with its groups, outliers, source
    const referenceTables = (ref, kind) => {
        const k = referenceKind(kind);
        const indent = { cellPadding: { top: 1.2, bottom: 1.2, left: 5, right: 1.5 } };
        const body = ref.sections.flatMap(sec => [
            [{ content: referenceSectionLabel(sec.section, parents), colSpan: 3, styles: { fontStyle: 'bold' } },
             num(formatInt(sec.count)), '', '', num(formatRmsz(sec.rmsz)), num(formatInt(sec.outliers))],
            ...sec.groups.map(grp => [{ content: grp.atoms, styles: indent }, grp.class || '',
                `${k.format(grp.reference)} ± ${k.format(grp.sigma)}`, formatInt(grp.count), k.format(grp.mean),
                k.formatDelta(grp.mean - grp.reference), formatRmsz(grp.rmsz), formatInt(grp.outliers)]),
        ]);
        table({
            head: [[k.item, gettext('Klasse'), num(k.reference), num(gettext('Anzahl')),
                    num(k.mean), num(k.delta), num(gettext('RMSZ')), num(gettext('Ausreißer'))]],
            body,
            foot: [[gettext('Gesamt'), '', '', num(formatInt(ref.rated)), '', '', num(formatRmsz(ref.rmsz)), num(formatInt(ref.outlier_count))]],
            columnStyles: right([2, 3, 4, 5, 6, 7]),
        });
        if (ref.outliers.length) {
            subTitle(interpolate(gettext('Ausreißer (|Z| > %s)'), [ref.outlier_threshold]));
            table({
                head: [[gettext('Atome'), gettext('Kette'), gettext('Pos.'), gettext('Aminosäure'), num(k.value), num(k.reference), num('Z')]],
                body: ref.outliers.map(o => {
                    const d = k.describe(o.atoms.map(nodeLikeAtom));
                    return [d.atoms, d.chain, d.pos, d.residues, k.format(o.value),
                            `${k.format(o.reference)} ± ${k.format(o.sigma)}`, formatZ(o.z)];
                }),
                columnStyles: right([4, 5, 6]),
            });
        }
        ensureSpace(8);
        doc.setFont('helvetica', 'normal').setFontSize(7.5).setTextColor(...MUTED)
            .text(pdfText(referenceSourceText(ref.sources)), M, y - 3);
        y += 4;
    };

    // Header
    doc.setFont('helvetica', 'bold').setFontSize(16).setTextColor(...INK).text(pdfText(interpolate(gettext('Statistiken – %s'), [data.file])), M, y + 4);
    y += 10;
    const created = new Date(data.created).toLocaleString(LOCALE, { dateStyle: 'medium', timeStyle: 'short' });
    doc.setFont('helvetica', 'normal').setFontSize(9).setTextColor(...MUTED)
        .text(pdfText(interpolate(gettext('%(model)s · erstellt am %(date)s'), { model: data.model, date: created }, true)), M, y);
    y += 10;

    // General
    const g = data.general;
    sectionTitle(gettext('Allgemein'));
    table({
        head: [[gettext('Models'), gettext('Aminosäuren'), gettext('Atome'), gettext('Ketten'), gettext('Bindungen')]],
        body: [[formatInt(g.models), formatInt(g.amino_acids), formatInt(g.atoms), formatInt(g.chains), formatInt(g.bonds)]],
        tableWidth: contentW * 0.75,
    });
    const info = structureInfo();
    if (info && (info.method || info.resolution)) {
        const fixed = (v, digits) => v === null ? '–' : v.toLocaleString(LOCALE, { minimumFractionDigits: digits, maximumFractionDigits: digits });
        table({
            head: [[gettext('Methode'), gettext('Auflösung'), 'R-work', 'R-free']],
            body: [[info.method_label || '–', info.resolution === null ? '–' : `${fixed(info.resolution, 2)} Å`,
                    fixed(info.r_work, 3), fixed(info.r_free, 3)]],
            tableWidth: contentW * 0.75,
        });
    }
    subTitle(gettext('Elemente'));
    table({
        head: [[gettext('Element'), num(gettext('Anzahl')), num(gettext('Anteil'))]],
        body: g.elements.map(e => [e.element, formatInt(e.count), formatPercent(e.share_percent)]),
        foot: [[gettext('Atome gesamt'), num(formatInt(g.atoms)), num(formatPercent(100))]],
        columnStyles: right([1, 2]),
        tableWidth: contentW * 0.6,
    });
    subTitle(gettext('Aminosäuren'));
    table({
        head: [[gettext('Aminosäure'), gettext('Code'), num(gettext('Anzahl')), num(gettext('Anteil'))]],
        body: g.amino_acid_composition.map(a => [a.name, `${a.code} · ${a.one_letter}`, formatInt(a.count), formatPercent(a.share_percent)]),
        foot: [[gettext('Aminosäuren gesamt'), '', num(formatInt(g.amino_acids)), num(formatPercent(100))]],
        columnStyles: right([2, 3]),
        tableWidth: contentW * 0.6,
    });

    // Elements / amino acids per chain
    const chainSection = ({ title, matrixTitle, shareTitle, rowLabel, rows, keyOf, labelOf, chains, countName, totalLabel, total }) => {
        sectionTitle(title);
        subTitle(matrixTitle);
        const lookup = chains.map(c => new Map(c.composition.map(item => [keyOf(item), item])));
        table({
            head: [
                [{ content: rowLabel, rowSpan: 2, styles: { valign: 'bottom' } },
                    ...chains.map(c => ({ content: chainLabel(c.chain), colSpan: 2, styles: { halign: 'center' } }))],
                chains.flatMap(() => [num(gettext('Anzahl')), num(gettext('Anteil'))]),
            ],
            body: rows.map(row => [labelOf(row), ...lookup.flatMap(map => {
                const item = map.get(keyOf(row));
                return [num(formatInt(item ? item.count : 0)), num(formatPercent(item ? item.share_percent : 0))];
            })]),
            foot: [[gettext('Gesamt'), ...chains.flatMap(c => [num(formatInt(c[countName])), num(formatPercent(100))])]],
        });
        subTitle(shareTitle);
        table({
            head: [[gettext('Kette'), num(gettext('Anzahl')), num(gettext('Anteil'))]],
            body: chains.map(c => [chainLabel(c.chain), formatInt(c[countName]), formatPercent(c.share_percent)]),
            foot: [[totalLabel, num(formatInt(total)), num(formatPercent(100))]],
            columnStyles: right([1, 2]),
            tableWidth: contentW * 0.6,
        });
    };
    chainSection({
        title: gettext('Elemente'), matrixTitle: gettext('Elemente je Kette'), shareTitle: gettext('Atome je Kette'), rowLabel: gettext('Element'),
        rows: g.elements, keyOf: e => e.element, labelOf: e => e.element,
        chains: data.elements_by_chain, countName: 'atoms', totalLabel: gettext('Atome gesamt'), total: g.atoms,
    });
    if (data.amino_acids_by_chain.length) {
        chainSection({
            title: gettext('Aminosäuren'), matrixTitle: gettext('Zusammensetzung je Kette'), shareTitle: gettext('Aminosäuren je Kette'), rowLabel: gettext('Aminosäure'),
            rows: g.amino_acid_composition, keyOf: a => a.code, labelOf: a => a.name,
            chains: data.amino_acids_by_chain, countName: 'amino_acids', totalLabel: gettext('Aminosäuren gesamt'), total: g.amino_acids,
        });
    }

    // Bonds
    const b = data.bonds;
    sectionTitle(gettext('Bindungen'));
    subTitle(gettext('Abweichung von Referenzwerten'));
    referenceTables(b.reference, 'bond');
    subTitle(gettext('Peptidbindungen'));
    const p = b.peptide_bonds;
    table({
        head: [[gettext('Peptidbindung'), num(gettext('Anzahl')), num(gettext('Anteil'))]],
        body: [
            [gettext('Existent'), formatInt(p.formed), formatPercent(sharePercent(p.formed, p.possible || 1))],
            [`${gettext('Nicht-existent')} (${gettext('Kettenlücken')})`, formatInt(p.missing),
             formatPercent(sharePercent(p.missing, p.possible || 1))],
        ],
        foot: [[gettext('Theoretisch möglich'), num(formatInt(p.possible)), num(formatPercent(100))]],
        columnStyles: right([1, 2]),
        tableWidth: contentW * 0.6,
    });
    if (b.disulfides.length) {
        subTitle(gettext('Disulfidbrücken'));
        table({
            head: [[gettext('Kette'), gettext('Pos.'), num(gettext('Länge (Å)')), gettext('Herkunft')]],
            body: b.disulfides.map(ss => {
                const d = describeBond(ss.atoms.map(nodeLikeAtom));
                return [d.chain, d.pos, num(formatLength(ss.length)), disulfideOrigin(ss.via)];
            }),
            tableWidth: contentW * 0.6,
        });
    }
    subTitle(gettext('Kürzeste und längste Bindung'));
    const comparison = (value, length) => value === null || value === undefined ? '–'
        : `${formatLength(value)} (${formatDelta(length - value)})`;
    const extremeRows = (title, pair) => [
        [{ content: title, colSpan: 10, styles: { fontStyle: 'bold', textColor: MUTED } }],
        ...[[gettext('Kürzeste Bindung'), pair?.shortest], [gettext('Längste Bindung'), pair?.longest]].map(([label, bond]) => {
            if (!bond) return [label, { content: gettext('Keine Bindungen'), colSpan: 9 }];
            const d = describeBond(bond.atoms.map(nodeLikeAtom));
            return [label, bond.model, num(formatLength(bond.length)), d.atoms, d.chain, d.pos, d.serials, d.residues,
                    num(comparison(bond.ensemble_mean, bond.length)), num(comparison(bond.ensemble_median, bond.length))];
        }),
    ];
    table({
        styles: { fontSize: 7, cellPadding: { top: 1.2, bottom: 1.2, left: 1.2, right: 1.2 }, textColor: INK, lineColor: LINE },
        head: [['', gettext('Model'), num(gettext('Länge (Å)')), gettext('Atome'), gettext('Kette'), gettext('Pos.'), gettext('Atom-Nr.'), gettext('Aminosäure'), num(gettext('Durchschnitt (Å)')), num(gettext('Median (Å)'))]],
        columnStyles: { 2: { cellPadding: { top: 1.2, bottom: 1.2, left: 1.2, right: 4 } } },
        body: [...extremeRows(gettext('Alle Models'), b.extremes.all_models), ...extremeRows(gettext('Ausgewähltes Model'), b.extremes.selected_model)],
    });

    if (b.histogram.bins.length) {
        subTitle(gettext('Verteilung der Bindungslängen'));
        const h = 62;
        ensureSpace(h + 4);
        drawPdfHistogram(doc, b.histogram.bins, M, y + 2, contentW, h, { INK, MUTED, LINE, PRIMARY });
        y += h + 9;
    }

    subTitle(gettext('Bindungen nach beteiligten Atomen'));
    table({
        head: [[gettext('Atome'), gettext('Ketten'), num(gettext('Anzahl')), num(gettext('Durchschnitt (Å)')), num(gettext('Median (Å)')), num(gettext('Min (Å)')), num(gettext('Max (Å)')), num(gettext('Δ Max − Min (Å)'))]],
        body: b.by_atoms.map(grp => [grp.atoms.join('–'), grp.chains.map(c => c || '–').join(', '), formatInt(grp.count),
            formatLength(grp.mean), formatLength(grp.median), formatLength(grp.min), formatLength(grp.max), formatLength(grp.delta)]),
        foot: b.count ? [[gettext('Bindungen gesamt'), '', num(formatInt(b.count)), num(formatLength(b.mean)), num(formatLength(b.median)),
            num(formatLength(b.min)), num(formatLength(b.max)), num(formatLength(b.delta))]] : undefined,
        columnStyles: {
            ...right([2, 3, 4, 7]),
            5: { halign: 'right', textColor: LENGTH_MIN },
            6: { halign: 'right', textColor: LENGTH_MAX },
        },
    });
    // Individual bonds per atom pair (as expanded in the tab); shortest and longest highlighted in colour
    if (options.bondList) {
        subTitle(gettext('Einzelne Bindungen je Atompaar'));
        table({
            head: [[num(gettext('Länge (Å)')), num(gettext('Abw. vom Ø (Å)')), gettext('Kette'), gettext('Pos.'),
                    gettext('Atom-Nr.'), gettext('Aminosäure'), '']],
            body: b.by_atoms.flatMap(grp => [
                [{ content: `${grp.atoms.join('–')} (${formatInt(grp.count)})`, colSpan: 7, styles: { fontStyle: 'bold' } }],
                ...grp.bonds.map(bond => {
                    const d = describeBond(bond.atoms.map(nodeLikeAtom));
                    const extreme = grp.min === grp.max ? null
                        : bond.length === grp.min ? ['Min', LENGTH_MIN] : bond.length === grp.max ? ['Max', LENGTH_MAX] : null;
                    const mark = extreme ? { content: gettext(extreme[0]), styles: { textColor: extreme[1], fontStyle: 'bold' } } : '';
                    const length = extreme ? { content: formatLength(bond.length), styles: { halign: 'right', textColor: extreme[1] } }
                        : num(formatLength(bond.length));
                    return [length, num(formatDelta(bond.length - grp.mean)), d.chain, d.pos, d.serials, d.residues, mark];
                }),
            ]),
            columnStyles: right([0, 1]),
        });
    }

    // Angles: Ramachandran plot with legend next to it
    const ang = data.angles;
    const RAMA = { general: [13, 110, 253], gly: [240, 140, 0], pro: [47, 158, 68] };
    sectionTitle(gettext('Winkel'));
    subTitle(gettext('Ramachandran-Plot'), 100);
    const ramaSize = 95;
    ensureSpace(ramaSize + 4);
    const OUTLIER = [224, 49, 49];
    drawPdfRamachandran(doc, ang.ramachandran, M, y + 2, ramaSize, { MUTED, LINE, OUTLIER, ...RAMA });
    RAMA_CLASSES.forEach((cls, i) => {
        const ly = y + 8 + i * 6;
        doc.setFillColor(...RAMA[cls]).circle(M + ramaSize + 10, ly - 1, 1, 'F');
        doc.setFont('helvetica', 'normal').setFontSize(8).setTextColor(...INK).text(pdfText(
            `${ramaClassLabel(cls)} (${formatInt(ang.ramachandran.filter(p => p.class === cls).length)})`), M + ramaSize + 13, ly);
    });
    const outlierY = y + 8 + RAMA_CLASSES.length * 6;
    doc.setDrawColor(...OUTLIER).setLineWidth(0.3).circle(M + ramaSize + 10, outlierY - 1, 1.3, 'S');
    doc.setFont('helvetica', 'normal').setFontSize(8).setTextColor(...INK).text(pdfText(
        `${ramaRegionLabel('outlier')} (${formatInt(ang.ramachandran_evaluation.outliers)})`), M + ramaSize + 13, outlierY);
    y += ramaSize + 9;

    // Ramachandran evaluation
    const ev = ang.ramachandran_evaluation;
    subTitle(gettext('Ramachandran-Bewertung'));
    const regionCells = t => [t.favored, t.allowed, t.outliers].map(v => num(formatInt(v)));
    table({
        head: [[gettext('Kategorie'), num(gettext('Anzahl')), ...RAMA_REGIONS.map(r => num(ramaRegionLabel(r)))]],
        body: ev.by_type.map(t => [ramaTypeLabel(t.type), num(formatInt(t.count)), ...regionCells(t)]),
        foot: [[gettext('Gesamt'), num(formatInt(ev.total)), ...regionCells(ev)]],
        tableWidth: contentW * 0.75,
    });
    const ramaOutliers = ang.ramachandran.filter(p => p.region === 'outlier');
    if (ramaOutliers.length) {
        subTitle(gettext('Ausreißer'));
        table({
            head: [[gettext('Kette'), gettext('Pos.'), gettext('Aminosäure'), gettext('Kategorie'), num('φ'), num('ψ')]],
            body: ramaOutliers.map(p => [p.chain || '–', p.position, p.residue, ramaTypeLabel(p.type),
                num(formatAngle(p.phi)), num(formatAngle(p.psi))]),
            tableWidth: contentW * 0.75,
        });
    }
    ensureSpace(8);
    doc.setFont('helvetica', 'normal').setFontSize(7.5).setTextColor(...MUTED)
        .text(pdfText(interpolate(gettext('Ramachandran-Bereiche: %s'), [ev.source || ''])), M, y - 3);
    y += 4;
    subTitle(gettext('cis- und trans-Peptidbindungen'));
    table({
        head: [[gettext('Peptidbindungen'), num('trans'), num('cis'), num(gettext('davon vor Prolin')), num(gettext('verdreht'))]],
        body: [[formatInt(ang.omega.count), formatInt(ang.omega.trans), formatInt(ang.omega.cis),
                formatInt(ang.omega.cis_before_proline), formatInt(ang.omega.twisted)]],
        columnStyles: right([1, 2, 3, 4]),
        tableWidth: contentW * 0.75,
    });
    if (ang.omega.non_trans.length) {
        table({
            head: [[gettext('Kette'), gettext('Pos.'), gettext('Aminosäuren'), num('ω'), gettext('Art')]],
            body: ang.omega.non_trans.map(o => [o.to.chain || '–', `${o.from.position}–${o.to.position}`,
                `${o.from.residue}–${o.to.residue}`, formatAngle(o.omega),
                o.kind === 'cis' ? (o.to.residue === 'PRO' ? gettext('cis (vor Prolin)') : gettext('cis')) : gettext('verdreht')]),
            columnStyles: right([3]),
            tableWidth: contentW * 0.75,
        });
    }
    subTitle(gettext('Bindungswinkel'));
    const fixed1 = v => v === null ? '–' : v.toLocaleString(LOCALE, { minimumFractionDigits: 1, maximumFractionDigits: 1 });
    table({
        head: [[gettext('Winkel'), num(gettext('Anzahl')), num(gettext('Durchschnitt (°)')), num(gettext('Std.-Abw. (°)')),
                num(gettext('Min (°)')), num(gettext('Max (°)'))]],
        body: ang.bond_angles.map(t => [t.angle, formatInt(t.count), fixed1(t.mean), fixed1(t.sd), fixed1(t.min), fixed1(t.max)]),
        columnStyles: right([1, 2, 3, 4, 5]),
    });
    subTitle(gettext('Abweichung von Referenzwerten'));
    referenceTables(ang.reference, 'angle');

    // Footer on every page (added last, when the page count is known)
    const pages = doc.getNumberOfPages();
    for (let i = 1; i <= pages; i++) {
        doc.setPage(i);
        doc.setFont('helvetica', 'normal').setFontSize(7.5).setTextColor(...MUTED);
        doc.text(pdfText(`Bondamino · ${data.file} · ${data.model}`), M, pageH - 8);
        doc.text(pdfText(interpolate(gettext('Seite %(page)s von %(pages)s'), { page: i, pages }, true)),
            pageW - M, pageH - 8, { align: 'right' });
    }
    return doc;
}

/**
 * Draws the histogram as vector graphics into the PDF (same design as in the tab).
 * @param {!Object} doc jsPDF document.
 * @param {!Array<!Object>} bins Exported bins (`from`, `to`, `count`, `share_percent`).
 * @param {number} left Left edge (mm).
 * @param {number} top Top edge (mm).
 * @param {number} width Width (mm).
 * @param {number} height Height (mm).
 * @param {!Object<string, !Array<number>>} colors RGB colours INK, MUTED, LINE, PRIMARY.
 */
function drawPdfHistogram(doc, bins, left, top, width, height, colors) {
    const plotL = left + 16, plotR = left + width, plotT = top + 2, plotB = top + height - 11;
    const x = d3.scaleLinear().domain([bins[0].from, bins[bins.length - 1].to]).range([plotL, plotR]);
    const y = d3.scaleLinear().domain([0, d3.max(bins, b => b.share_percent)]).nice().range([plotB, plotT]);

    doc.setFont('helvetica', 'normal').setFontSize(7).setTextColor(...colors.MUTED);
    doc.setDrawColor(...colors.LINE).setLineWidth(0.15);
    y.ticks(5).forEach(t => {
        doc.line(plotL, y(t), plotR, y(t));
        doc.text(pdfText(formatPercentTick(t)), plotL - 2, y(t) + 1, { align: 'right' });
    });

    doc.setFillColor(...colors.PRIMARY);
    bins.filter(b => b.count).forEach(b => {
        const h = Math.max(0.3, plotB - y(b.share_percent));
        doc.rect(x(b.from) + 0.2, plotB - h, Math.max(0.2, x(b.to) - x(b.from) - 0.4), h, 'F');
    });

    doc.setDrawColor(...colors.MUTED).setLineWidth(0.2).line(plotL, plotB, plotR, plotB);
    x.ticks(Math.max(2, Math.floor((plotR - plotL) / 18))).forEach(t => {
        doc.line(x(t), plotB, x(t), plotB + 1);
        doc.text(formatAxisLength(t), x(t), plotB + 4, { align: 'center' });
    });
    doc.setFontSize(8);
    doc.text(pdfText(gettext('Bindungslänge (Å)')), (plotL + plotR) / 2, top + height - 1, { align: 'center' });
    const yLabel = pdfText(gettext('Anteil der Bindungen'));
    doc.text(yLabel, left + 2.5, (plotT + plotB) / 2 + doc.getTextWidth(yLabel) / 2, { angle: 90 });
}

// ===== Init =====
// Entry point: reads the embedded JSON, renders all tabs once and wires up the
// controls (model/conformation selects, style/layout buttons, exports, tab
// switches, theme and resize observers). Runs only on pages with #graph-data.
(function () {
    const dataEl = document.getElementById('graph-data');
    if (!dataEl) return;
    // Current graph data and bond ensemble; replaced by loadStructure().
    let graphData = JSON.parse(dataEl.textContent);
    const ensembleEl = document.getElementById('bond-ensemble');
    let bondEnsemble = ensembleEl ? JSON.parse(ensembleEl.textContent) : null;

    // Bond colouring and layout of the graph, remembered across visits.
    let activeStyle  = localStorage.getItem('bondStyle') || 'gradient';
    let activeLayout = localStorage.getItem('layout')    || 'schematic';

    const canvas = setupCanvas();

    // Re-renders the graph tab with the current data, style and layout.
    function redraw(animate = true) {
        render(canvas, graphData, activeStyle, activeLayout, animate);
    }

    const selectsWrap = document.getElementById('structure-selects');
    const modelSelect = document.getElementById('model-select');
    const altlocSelect = document.getElementById('altloc-select');
    // Disabled selects (only one model / no alternate locations) offer no choice.
    const selectable = [modelSelect, altlocSelect].filter(s => s && !s.disabled);

    // Label of the selected model, with the conformation id if one is selected.
    function modelLabel() {
        const model = modelSelect
            ? modelSelect.selectedOptions[0].textContent.trim()
            : interpolate(gettext('Model %s'), [bondEnsemble?.models[0] ?? 1]);
        const altloc = altlocSelect && !altlocSelect.disabled ? altlocSelect.value : '';
        return altloc ? `${model} · ${interpolate(gettext('Konformation %s'), [altloc])}` : model;
    }
    // Model of the currently loaded data (the select changes before loading finishes).
    let loadedModelLabel = modelLabel();

    // Loads graph data for the selected combination of model and conformation from
    // the `/pdb/<pk>/models/` endpoint and re-renders all tabs; after a conformation
    // change also the bond ensemble over all models (it depends on the altloc).
    // The selects are disabled while loading.
    async function loadStructure(withEnsemble) {
        const params = new URLSearchParams();
        if (modelSelect && !modelSelect.disabled) params.set('model', modelSelect.value);
        if (altlocSelect && !altlocSelect.disabled) params.set('altloc', altlocSelect.value);
        if (withEnsemble) params.set('ensemble', '1');
        selectable.forEach(s => { s.disabled = true; });
        try {
            const resp = await fetch(`${selectsWrap.dataset.modelUrl}?${params}`);
            if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
            const data = await resp.json();
            graphData = data.graph;
            if (data.ensemble) bondEnsemble = data.ensemble;
            loadedModelLabel = modelLabel();
            redraw();
            renderDataTable(canvas, graphData.nodes);
            renderStats(graphData);
            renderLegend(graphData.nodes);
            renderBondStats(graphData, bondEnsemble, loadedModelLabel);
            renderAngles(canvas, graphData);
        } catch (err) {
            console.error('Model/Konformation konnte nicht geladen werden:', err);
        } finally {
            selectable.forEach(s => { s.disabled = false; });
        }
    }
    selectable.forEach(select => select.addEventListener('change', () => loadStructure(select === altlocSelect)));

    // Export sidebar of the statistics tab: JSON, PDF, CSV and histogram images.
    const exportBar = document.getElementById('stats-export');
    if (exportBar) {
        const { filename, jspdfSrc, autotableSrc } = exportBar.dataset;
        const exportData = () => buildStatsExport(graphData, bondEnsemble, loadedModelLabel, filename);

        document.getElementById('export-json').addEventListener('click', () => {
            const data = exportData();
            downloadBlob(new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' }),
                `${exportBaseName(filename, data.model)}.json`);
        });

        const pdfButton = document.getElementById('export-pdf');
        // "Individual bonds" checkbox for the PDF; the choice is remembered in the browser
        const pdfBondList = document.getElementById('export-pdf-bonds');
        if (pdfBondList) {
            try { pdfBondList.checked = localStorage.getItem('pdfBondList') === '1'; } catch { /* no storage: off */ }
            pdfBondList.addEventListener('change', () => {
                try { localStorage.setItem('pdfBondList', pdfBondList.checked ? '1' : '0'); } catch { /* ignore */ }
            });
        }
        pdfButton.addEventListener('click', async () => {
            pdfButton.disabled = true;
            try {
                await loadScript(jspdfSrc);
                await loadScript(autotableSrc);
                const data = exportData();
                buildStatsPdf(data, { bondList: pdfBondList?.checked }).save(`${exportBaseName(filename, data.model)}.pdf`);
            } catch (err) {
                console.error('PDF konnte nicht erstellt werden:', err);
            } finally {
                pdfButton.disabled = false;
            }
        });

        exportBar.querySelectorAll('[data-csv]').forEach(button => {
            button.addEventListener('click', () => {
                const kind = button.dataset.csv;
                downloadBlob(new Blob([csvTable(kind, graphData)], { type: 'text/csv;charset=utf-8' }),
                    `${exportBaseName(filename, loadedModelLabel, CSV_FILE_KINDS[kind])}.csv`);
            });
        });

        exportBar.querySelectorAll('[data-histogram-format]').forEach(button => {
            button.addEventListener('click', async () => {
                const format = button.dataset.histogramFormat;
                if (!graphData.edges.length) return;
                button.disabled = true;
                try {
                    const extension = format === 'jpeg' ? 'jpg' : format;
                    downloadBlob(await histogramBlob(graphData.edges, format),
                        `${exportBaseName(filename, loadedModelLabel, gettext('histogramm'))}.${extension}`);
                } catch (err) {
                    console.error('Histogramm konnte nicht exportiert werden:', err);
                } finally {
                    button.disabled = false;
                }
            });
        });
    }

    // Syncs the active state of the style/layout buttons and shows the legend
    // matching the bond style.
    function updateButtons() {
        document.querySelectorAll('.style-option[data-style]').forEach(btn => {
            btn.classList.toggle('active', btn.dataset.style === activeStyle);
        });
        const referenceLegend = document.getElementById('legend-reference');
        if (referenceLegend) referenceLegend.hidden = activeStyle !== 'reference';
        const bFactorLegend = document.getElementById('legend-bfactor');
        if (bFactorLegend) bFactorLegend.hidden = activeStyle !== 'bfactor';
        const elementLegend = document.getElementById('legend-elements');
        if (elementLegend) elementLegend.hidden = activeStyle === 'bfactor';
        document.querySelectorAll('.style-option[data-layout]').forEach(btn => {
            btn.classList.toggle('active', btn.dataset.layout === activeLayout);
        });
    }

    // Global handlers for the onclick attributes of the toolbar buttons in graph.html.
    window.setBondStyle = function (style) {
        activeStyle = style;
        localStorage.setItem('bondStyle', style);
        updateButtons();
        redraw();
    };

    window.setLayout = function (layout) {
        activeLayout = layout;
        localStorage.setItem('layout', layout);
        updateButtons();
        redraw();
    };

    window.resetZoom = function () {
        canvas.svg.transition().duration(400).call(canvas.zoom.transform, d3.zoomIdentity);
    };

    // Initial render of all tabs (hidden tabs are rendered too; the histogram
    // catches up once its tab is visible, see below).
    updateButtons();
    redraw(false);
    renderDataTable(canvas, graphData.nodes);
    renderStats(graphData);
    renderLegend(graphData.nodes);
    renderBondStats(graphData, bondEnsemble, modelLabel());
    renderAngles(canvas, graphData);

    const anglesExportBar = document.getElementById('angles-export');
    // Category filter of the Ramachandran plot; loads the region grids on first use.
    const ramaTypeSelect = document.getElementById('rama-type');
    if (ramaTypeSelect) {
        d3.select(ramaTypeSelect).selectAll('option.rama-type-option').data(RAMA_TYPES).join('option')
            .attr('class', 'rama-type-option').attr('value', t => t).text(t => ramaTypeLabel(t));
        ramaTypeSelect.addEventListener('change', async () => {
            ramaView.type = ramaTypeSelect.value;
            try {
                if (ramaView.type) await loadRamaRegions();
            } catch (err) {
                console.error('Ramachandran-Bereiche konnten nicht geladen werden:', err);
            }
            if (canvas.ramaAngles) renderRamachandranPlot(canvas, canvas.ramaAngles);
        });
    }
    if (anglesExportBar) {
        anglesExportBar.querySelectorAll('[data-ramachandran-format]').forEach(button => {
            button.addEventListener('click', async () => {
                const format = button.dataset.ramachandranFormat;
                button.disabled = true;
                try {
                    const extension = format === 'jpeg' ? 'jpg' : format;
                    downloadBlob(await ramachandranBlob(graphData, format),
                        `${exportBaseName(anglesExportBar.dataset.filename, loadedModelLabel, gettext('ramachandran'))}.${extension}`);
                } catch (err) {
                    console.error('Ramachandran-Plot konnte nicht exportiert werden:', err);
                } finally {
                    button.disabled = false;
                }
            });
        });
    }

    // Scrolling the table updates the highlight in the sequence bar (at most once per frame).
    let tableScrollPending = false;
    canvas.tableScroll.addEventListener('scroll', () => {
        if (tableScrollPending) return;
        tableScrollPending = true;
        requestAnimationFrame(() => {
            tableScrollPending = false;
            if (canvas.state.activeTab === 'data') updateTableView(canvas);
        });
    });

    // The sequence bar belongs to the graph and data tabs only. When a tab has
    // been shown, catch up on a pending jump to the selected residue and refresh
    // the sequence highlight.
    const sequenceCard = document.getElementById('sequence-card');
    document.querySelectorAll('[data-bs-toggle="tab"][data-tab]').forEach(btn => {
        btn.addEventListener('show.bs.tab', () => {
            sequenceCard.hidden = btn.dataset.tab !== 'graph' && btn.dataset.tab !== 'data';
        });
        btn.addEventListener('shown.bs.tab', () => {
            const tab = btn.dataset.tab;
            const { state } = canvas;
            state.activeTab = tab;
            const res = state.pendingJump.has(tab) ? state.selectedResidue : null;
            state.pendingJump.delete(tab);
            if (tab === 'graph') {
                if (res) jumpGraphToResidue(canvas, res, 0);
                updateSequenceView(canvas, d3.zoomTransform(canvas.svg.node()));
            } else if (tab === 'data') {
                if (res) scrollTableToResidue(canvas, res, false);
                updateTableView(canvas);
            } else {
                markSequenceInView(canvas, () => false);
            }
            if (tab === 'stats') {
                // Catch up on the histogram skipped while the tab was hidden.
                const histogramEl = document.getElementById('stat-bond-histogram');
                if (histogramEl && histogramEl.clientWidth !== histogram.width) renderBondHistogram();
            }
        });
    });

    // Theme switch (data-bs-theme): redraw, because the graph reads colours from
    // CSS variables at render time.
    const observer = new MutationObserver(() => redraw(false));
    observer.observe(document.documentElement, { attributes: true, attributeFilter: ['data-bs-theme'] });

    // Content width and zoom bounds depend on the container size.
    // Do not redraw in a hidden tab (size 0), otherwise the zoom bounds get shifted.
    new ResizeObserver(() => {
        if (canvas.svg.node().clientWidth > 0) redraw(false);
    }).observe(canvas.svg.node());

    // Fit the histogram to its width; in a hidden tab draw it only once it is shown.
    const histogramEl = document.getElementById('stat-bond-histogram');
    if (histogramEl) {
        new ResizeObserver(() => {
            if (histogramEl.clientWidth !== histogram.width) renderBondHistogram();
        }).observe(histogramEl);
    }
})();
