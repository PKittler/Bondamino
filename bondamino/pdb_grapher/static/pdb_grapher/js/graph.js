'use strict';

// ===== Constants =====
const BAND_HEIGHT = 200;
const LABEL_W = 90;
const BACKBONE_X = LABEL_W + 80;
const SC_COL_W = 82;
const NODE_R = 16;
const LEAF_H = 46;
const PROP_SCALE = 55; // px per Å in proportional layout
const SC_RING_SIDE = 60; // Kantenlänge der Ringe im schematischen Layout (px)
const RESIDUE_SEPARATOR_W = 100;

const BB_Y = { N: 30, CA: 90, C: 150, O: 160 };
const BACKBONE_ATOMS = new Set(['N', 'CA', 'C', 'O', 'OXT']);

// Ringförmige Seitenketten: Atomnamen in Ringreihenfolge (letzter schließt zum ersten).
// TRP hat zwei fusionierte Ringe, die sich die Kante CD2-CE2 teilen.
// PRO's Ring bindet zurück in die Backbone-Atome N/CA.
const SIDECHAIN_RINGS = {
    HIS: [['CG', 'ND1', 'CE1', 'NE2', 'CD2']],
    PHE: [['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']],
    TYR: [['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']],
    TRP: [['CG', 'CD1', 'NE1', 'CE2', 'CD2'], ['CD2', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3']],
    PRO: [['N', 'CA', 'CB', 'CG', 'CD']],
};

// ===== Helpers =====
// Sprache der Seite (<html lang>) für Zahlen- und Datumsformate; Texte kommen aus dem
// JavaScript-Katalog von Django (gettext, ngettext, interpolate).
const LOCALE = document.documentElement.lang || 'de';

function cssVar(name) {
    return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

// Farbe aus der CSS-Variable --atom-<Element>; Elemente ohne eigene Farbe erhalten --atom-default.
function getElementColor(element) {
    return cssVar(`--atom-${(element || '').toUpperCase()}`) || cssVar('--atom-default');
}

// ===== Ring-Geometrie =====
// Im Proportional-Layout mit echten Bindungslängen, im schematischen Layout als
// regelmäßige Vielecke (`bondPx` liefert dann eine konstante Kantenlänge).
// Berechnet die Eckpunkte eines geschlossenen Vielecks mit vorgegebenen Kantenlängen,
// bei dem alle Ecken auf einem gemeinsamen Kreis liegen (zyklisches Polygon). So bekommt
// jede Kante — inklusive der ringschließenden — exakt ihre echte, proportionale Länge,
// ohne dass eine Kante gestreckt werden muss, um den Ring zu schließen.
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

function rotatePoint(p, cosA, sinA) {
    return { x: p.x * cosA - p.y * sinA, y: p.x * sinA + p.y * cosA };
}

// Hängt einen einzelnen Ring-Eckpunkt an eine bereits feste Weltposition, orientiert so,
// dass der Ring von dort aus in Richtung `outDir` "aufgeht" (weg vom Backbone).
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

// Hängt zwei benachbarte Ring-Eckpunkte an zwei bereits feste Weltpositionen
// (z. B. Backbone N/CA bei PRO, oder die gemeinsame Kante zweier fusionierter Ringe bei TRP).
// Da zwei Punkte den Ring nur bis auf Spiegelung festlegen, wählt `chooseMirror` die Seite.
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

// Legt alle Ringe einer Seitenkette fest (in Reihenfolge, damit ein zweiter, fusionierter
// Ring bereits platzierte Atome des ersten als Anker nutzen kann, siehe TRP).
function layoutSidechainRings(ringDefs, positions, scAdj, bondPx, pfx) {
    const priorRingPositions = []; // Atome vorheriger Ringe dieser Seitenkette (z. B. TRP-Ring 1)

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
                // Bei einem fusionierten Ring (z. B. TRP) so wählen, dass er vom bereits
                // platzierten ersten Ring wegzeigt statt ihn zu überlappen; sonst (z. B. PRO)
                // einfach die Seite mit größerem x (weiter weg vom Backbone).
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

// Dreht das Ringsystem einer Seitenkette um sein erstes Atom (z. B. CG), sodass sein
// Schwerpunkt auf der Höhe dieses Atoms liegt. Fusionierte Ringe (TRP) knicken sonst
// schräg nach unten ab und ragen weit in die nächste Aminosäure hinein.
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

// Platziert Atome, die von einem Ring abzweigen (z. B. das OH bei TYR), radial vom
// (angenäherten) Ringzentrum weg, im Abstand `bondPx`.
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

// ===== Layout =====
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
    function getBL(a, b) { return (bl[a] && bl[a][b]) ? bl[a][b] : 1.5; }

    const leafH = isProp ? 1.5 * PROP_SCALE : LEAF_H;
    let yPropCursor = 20;
    const bandBounds = [];

    graphData.residues.forEach(res => {
        const ri = res.index;
        const pfx = s => `${res.chain}_${res.pos}_${s}`;
        const nId = pfx('N'), caId = pfx('CA'), cId = pfx('C');
        const oId = pfx('O'), oxtId = pfx('OXT');
        const bandY = ri * BAND_HEIGHT;

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
            // Von O wegversetzt (Carboxylat-Sauerstoffe liegen sonst fast deckungsgleich).
            // Im Proportional-Layout per Winkel statt Pixel-Versatz, damit die
            // Kantenlänge weiterhin der echten Bindungslänge entspricht.
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

                // Ring-Atome werden separat mit Ring-Geometrie platziert und daher von der
                // normalen Baum-Traversierung ausgeklammert; als Baum gezeichnet würden sich
                // die ringschließenden Bindungen sonst überkreuzen (z. B. bei TRP).
                // Im schematischen Layout bleibt PRO beim Baum: Sein Ring läuft über die
                // im Raster fest platzierten Backbone-Atome N/CA.
                // Modifizierte Aminosäuren (z. B. PTR, HYP) übernehmen die Ringe der Ausgangs-
                // Aminosäure, sofern die Atome gleich heißen.
                const ringDefs = (SIDECHAIN_RINGS[res.name] || SIDECHAIN_RINGS[res.parent] || [])
                    .filter(ring => ring.every(name => atomById[pfx(name)]))
                    .filter(ring => isProp || !ring.some(name => BACKBONE_ATOMS.has(name)));
                const bondPx = isProp ? (a, b) => getBL(a, b) * PROP_SCALE : () => SC_RING_SIDE;
                ringDefs.forEach(ringNames => ringNames.forEach(name => scVisited.add(pfx(name))));

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

                // Ring-Atome erst platzieren, nachdem ihr Anker (z. B. CB) oben eine Position bekommen hat.
                if (ringDefs.length > 0) {
                    layoutSidechainRings(ringDefs, positions, scAdj, bondPx, pfx);
                    // PRO-Ring nicht drehen: Er hängt an den festen Backbone-Atomen N/CA.
                    const anchoredToBackbone = ringDefs.some(ring => ring.some(name => BACKBONE_ATOMS.has(name)));
                    if (!anchoredToBackbone) balanceRingSystem(ringDefs, positions, pfx);
                    placeRingBranches(ringDefs, positions, atomById, scAdj, bondPx, pfx, resPrefix);
                }
            }
        }

        // --- Band bounds & cursor advance ---
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
const TRANSITION_MS = 600;
const VIEW_PAD = 40;

const residueKey = r => `${r.chain}_${r.pos}_${r.name}`;
const edgeKey = e => `${e.source}|${e.target}`;

// D3s Standard-Begrenzung, zusätzlich darf die Ansicht nie rechts vom linken Inhaltsrand
// beginnen (translate-x ≤ 0) — auch nicht beim Herauszoomen, wo D3 sonst zentriert.
const d3DefaultConstrain = d3.zoom().constrain();
function constrainLeftEdge(transform, extent, translateExtent) {
    const t = d3DefaultConstrain(transform, extent, translateExtent);
    return t.x > 0 ? d3.zoomIdentity.translate(0, t.y).scale(t.k) : t;
}

// Einmaliger Aufbau: Viewport-Gruppe (wird per Zoom transformiert) mit den Layern.
// Mausrad scrollt wie bisher, Strg + Mausrad (bzw. Pinch auf dem Trackpad) zoomt,
// Ziehen verschiebt die Ansicht.
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
        .filter(event => (event.type === 'wheel' ? event.ctrlKey : !event.ctrlKey) && !event.button)
        // D3 verzehnfacht bei gedrückter Strg-Taste das Delta (für Pinch-Gesten); hier
        // ein einheitlicher Faktor, damit eine Mausrad-Raste nicht zu stark zoomt.
        .wheelDelta(event => -event.deltaY * (event.deltaMode === 1 ? 0.05 : event.deltaMode ? 1 : 0.005))
        .constrain(constrainLeftEdge)
        .on('zoom', event => {
            viewport.attr('transform', event.transform);
            if (canvas && canvas.state.activeTab === 'graph') updateSequenceView(canvas, event.transform);
        });

    svg.call(zoom)
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
        // pendingJump: Ansichten, die beim nächsten Einblenden noch zur Auswahl springen müssen.
        state: { residues: [], bandBounds: [], selectedKey: null, activeTab: 'graph', pendingJump: new Set() },
    };
    return canvas;
}

function render(canvas, graphData, bondStyle, layout, animate) {
    const { positions, bandBounds } = computePositions(graphData, layout);

    const allPos = Object.values(positions);
    if (allPos.length === 0) return;

    const { svg, viewport, zoom } = canvas;
    const viewW = svg.node().clientWidth;
    const viewH = svg.node().clientHeight;
    // Spur für Disulfidbrücken rechts neben allen Atomen, eine Bahn je Brücke
    const disulfideCount = graphData.edges.filter(e => e.type === 'disulfide').length;
    const gutterX = d3.max(allPos, p => p.x) + NODE_R + SS_GUTTER_GAP;
    const maxX = d3.max(allPos, p => p.x) + NODE_R + 60 + (disulfideCount ? SS_GUTTER_GAP + disulfideCount * SS_LANE_W : 0);
    const maxY = d3.max(allPos, p => p.y) + NODE_R + 20;
    const contentW = Math.max(maxX, 500, viewW);
    const contentH = layout === 'proportional'
        ? maxY + 20
        : graphData.residues.length * BAND_HEIGHT + 10;

    canvas.state.residues = graphData.residues;
    canvas.state.bandBounds = bandBounds;
    renderSequence(canvas, graphData);

    zoom.translateExtent([[-VIEW_PAD, -VIEW_PAD], [contentW + VIEW_PAD, Math.max(contentH, viewH) + VIEW_PAD]]);
    // Aktuelle Transformation an die neuen Grenzen anpassen (z. B. nach Layoutwechsel).
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

    // Kurzer Trennstrich links an jeder Grenze zwischen zwei Aminosäuren (statt
    // farbiger Bänder, in die Seitenketten benachbarter Residuen hineinragen können).
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
    // Disulfidbrücken verbinden weit entfernte Residuen und werden als Bogen gezeichnet; die
    // Farbskala bleibt bei den übrigen Bindungen (SG–SG ist deutlich länger).
    const straightEdges = visibleEdges.filter(e => e.type !== 'disulfide');
    const arcEdges = visibleEdges.filter(e => e.type === 'disulfide');
    const [minLen, maxLen] = d3.extent(straightEdges, e => e.length);
    const hue = d3.scaleLinear().domain([minLen, maxLen]).range([0, 120]).clamp(true);
    const monoColor = cssVar('--bond-mono');

    const zColors = Object.fromEntries(['ok', 'warn', 'outlier', 'none'].map(c => [c, cssVar(`--z-${c}`)]));
    function bondColor(e) {
        if (bondStyle === 'mono') return monoColor;
        if (bondStyle === 'reference') return zColors[zClass(zScore(e))];
        return d3.hsl(minLen === maxLen ? 60 : hue(e.length), 0.8, 0.42).toString();
    }

    // --- Edges ---
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
        .attr('stroke-width', e => e.type === 'sidechain' ? 1.8 : 2.5)
        .attr('stroke-dasharray', e => e.type === 'sidechain' ? '5,3' : null)
        .call(sel => sel.select('title').text(e =>
            `${atomLabel(e.source)}–${atomLabel(e.target)}: ${e.length} Å (${BOND_TYPE_LABELS[e.type] || e.type})\n${referenceText(e)}`))
        .transition(t)
        .attr('stroke', e => bondColor(e))
        .attr('x1', e => positions[e.source].x).attr('y1', e => positions[e.source].y)
        .attr('x2', e => positions[e.target].x).attr('y2', e => positions[e.target].y);

    // Disulfidbrücke als Klammer: vom SG waagerecht nach rechts in die eigene Bahn der Spur,
    // dort senkrecht zum Partner und waagerecht zurück (abgerundete Ecken). So kreuzt sie keine
    // Seitenketten anderer Residuen. `mid` ist die Mitte des senkrechten Abschnitts.
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

    // Längenbeschriftung neben der Kantenmitte; bei zu kurzen Kanten ausgeblendet.
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
    // Bindungswinkel mit dem Atom in der Mitte (für den Tooltip)
    const anglesByCenter = analyzeAngles(graphData).byCenter;
    const nodeTransform = n => `translate(${positions[n.id].x},${positions[n.id].y})`;
    gNodes.selectAll('g.atom-group')
        .data(graphData.nodes.filter(n => positions[n.id]), n => n.id)
        .join(
            enter => {
                const g = enter.append('g').attr('class', 'atom-group').attr('transform', nodeTransform);
                g.append('title');
                g.append('circle').attr('class', 'atom-circle').attr('r', NODE_R)
                    .attr('fill', n => getElementColor(n.element));
                g.append('text').attr('class', 'atom-text');
                return g.call(fadeIn);
            },
            update => update,
            fadeOut,
        )
        .call(g => g.select('title').text(n => [
            `${n.label} · ${n.residue} ${n.residue_seq}${n.i_code || ''} (${interpolate(gettext('Chain %s'), [n.chain])})`,
            ...(anglesByCenter.get(n.id) || []).map(a => `${a.type}: ${formatAngle(a.value)}`),
        ].join('\n')))
        .call(g => g.select('text').text(n => n.label))
        .call(g => g.select('circle').transition(t).attr('fill', n => getElementColor(n.element)))
        .transition(t)
        .attr('transform', nodeTransform);
}

// ===== Sequenzleiste =====
// Spur der Disulfidbrücken im Graphen: Abstand zu den Atomen und Breite je Bahn
const SS_GUTTER_GAP = 40;
const SS_LANE_W = 16;

const SEQ_CELL_W = 16;
const SEQ_CELL_H = 22;
const SEQ_TICK_H = 14;
const SEQ_CHAIN_GAP = 10;
const SEQ_BREAK_GAP = 4;

// Peptidbindungen zwischen aufeinanderfolgenden Aminosäuren einer Chain: `possible` zählt alle
// solchen Paare, `breaks` die Kettenlücken darunter, also Paare ohne Peptidbindung, weil der
// Server den C–N-Abstand als zu groß verworfen hat (z. B. in der Struktur fehlende Residuen).
// `breaks` liefert je Lücke (Schlüssel: Index des Residuums danach) das Residuum davor und den
// Abstand. Liganden (ohne N/CA/C) zählen nicht.
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

// Breite eines Zeichens der Nummernbeschriftung (9px Roboto Mono) plus Mindestabstand
const SEQ_TICK_CHAR_W = 5.5;
const SEQ_TICK_GAP = 6;

// Nummern über der Sequenzleiste: Kettenanfang ("A 1") und jede Zehnerposition. Überlappt
// eine Beschriftung die vorherige, entfällt die Zehnerposition; Kettenanfänge bleiben immer.
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

// Mikroheterogenität: Die Konformationen (altLoc) des Residuums sind verschiedene Aminosäuren.
// Liefert z. B. "PRO (A) · SER (B, C)" oder null, wenn alle Konformationen dieselbe sind.
function residueVariants(res) {
    const byName = d3.groups(Object.entries(res.conformations || {}), ([, name]) => name);
    if (byName.length < 2) return null;
    return byName.map(([name, items]) => `${name} (${items.map(([altloc]) => altloc).join(', ')})`).join(' · ');
}

// Ein-Buchstaben-Sequenz als Zellenreihe; Chain-Wechsel bekommen eine Lücke, Kettenlücken
// innerhalb einer Chain einen roten Strich. Klick auf eine Zelle zentriert die zugehörige
// Aminosäure im Graphen.
function renderSequence(canvas, graphData) {
    const svg = d3.select('#sequence-svg');
    if (svg.empty()) return;

    const residues = graphData.residues;
    const { breaks } = analyzePeptideBonds(graphData);
    const bridges = disulfidePartners(graphData);
    // cis-Peptidbindungen, Schlüssel: Index des Residuums nach der Bindung
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
            // Mitte des Zwischenraums zwischen den beiden Zellen (Zellen sind SEQ_CELL_W - 2 breit)
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

    // Zusätzliche Breite, damit eine Tick-Beschriftung an der letzten Zelle nicht abgeschnitten wird.
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
            // Punkt oben rechts: modifizierte Aminosäure; oben links: Cystein mit Disulfidbrücke
            g.append('circle').attr('class', 'seq-mod-mark').attr('cx', SEQ_CELL_W - 5).attr('cy', 3.5).attr('r', 2);
            g.append('circle').attr('class', 'seq-ss-mark').attr('cx', 3).attr('cy', 3.5).attr('r', 2);
            // Markierung am unteren Rand: Konformationen sind verschiedene Aminosäuren
            g.append('rect').attr('class', 'seq-conf-mark')
                .attr('x', 2).attr('y', SEQ_CELL_H - 3).attr('width', SEQ_CELL_W - 6).attr('height', 2.5).attr('rx', 1.25);
            g.append('text').attr('x', (SEQ_CELL_W - 2) / 2).attr('y', SEQ_CELL_H / 2);
            return g;
        })
        .attr('class', c => `seq-cell${c.res.one_letter === '?' ? ' het' : ''}${c.res.parent ? ' modified' : ''}${residueVariants(c.res) ? ' variant' : ''}${bridges.has(residueKey(c.res)) ? ' disulfide' : ''}`)
        .attr('transform', c => `translate(${c.x},${SEQ_TICK_H})`)
        .on('click', (event, c) => jumpToResidue(canvas, c.res))
        .call(g => g.select('title').text(c => {
            const lines = [`${c.res.name} ${c.res.pos} (${interpolate(gettext('Chain %s'), [c.res.chain])})`];
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
            // Breitere, unsichtbare Fläche, damit der Tooltip leicht zu treffen ist
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

// Markiert die Zellen, deren Band gerade im Graph-Ausschnitt liegt.
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

// Hebt die Zellen hervor, für die `inView` zutrifft, und scrollt die Leiste mit,
// falls dieser Bereich aus ihrem sichtbaren Teil herausläuft.
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

function updateSelection(canvas) {
    const key = canvas.state.selectedKey;
    canvas.viewport.selectAll('g.residue-label').classed('selected', r => residueKey(r) === key);
    if (canvas.seqCells) canvas.seqCells.classed('selected', c => residueKey(c.res) === key);
    if (canvas.tableRows) canvas.tableRows.classed('selected', n => nodeResidueKey(n) === key);
    if (canvas.ramaPoints) canvas.ramaPoints.classed('selected', p => residueKey(p.res) === key).filter('.selected').raise();
}

// Springt in der sichtbaren Ansicht (Graph oder Tabelle) zur Aminosäure; die jeweils
// andere Ansicht holt den Sprung nach, sobald ihr Tab eingeblendet wird.
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

function jumpGraphToResidue(canvas, res, duration) {
    const b = canvas.state.bandBounds[res.index];
    if (!b) return;
    // Horizontale Position und Zoomstufe beibehalten, nur vertikal zum Band springen.
    const node = canvas.svg.node();
    const t = d3.zoomTransform(node);
    const mid = (b.yStart + b.yEnd) / 2;
    const dy = (node.clientHeight / 2 - (mid * t.k + t.y)) / t.k;
    canvas.svg.transition().duration(duration)
        .call(canvas.zoom.translateBy, 0, dy);
}

// ===== Datentabelle =====
const fixed = digits => v => (v === null || v === undefined ? '' : v.toFixed(digits));
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

const nodeResidueKey = n => `${n.chain}_${n.residue_seq}${n.i_code || ''}_${n.residue}`;

// Eine Zeile pro Atom, in der Reihenfolge der PDB-Datei; die erste Zeile jeder
// Aminosäure bekommt eine Trennlinie.
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

function scrollTableToResidue(canvas, res, smooth) {
    const row = canvas.firstRowByResidue && canvas.firstRowByResidue.get(res.index);
    if (!row) return;
    const headH = canvas.tableScroll.querySelector('thead').offsetHeight;
    canvas.tableScroll.scrollTo({ top: row.offsetTop - headH, behavior: smooth ? 'smooth' : 'auto' });
}

// Markiert in der Sequenzleiste die Aminosäuren, deren Zeilen gerade sichtbar sind.
function updateTableView(canvas) {
    const rows = canvas.tableRows ? canvas.tableRows.nodes() : [];
    if (rows.length === 0) return;
    const scroller = canvas.tableScroll;
    const top = scroller.scrollTop + scroller.querySelector('thead').offsetHeight;
    const bottom = scroller.scrollTop + scroller.clientHeight;

    // Zeilen sind nach offsetTop sortiert: erste Zeile, deren Unterkante unter `y` liegt
    // (1 px Toleranz, damit eine Zeile, die genau an der Kante endet, nicht mitzählt).
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

// ===== Statistiken =====
const formatInt = v => v.toLocaleString(LOCALE);
// Prozentwerte (0–100) im Format der Sprache, z. B. „12,5 %“ bzw. „12.5%“.
const formatPercent = v => (v / 100).toLocaleString(LOCALE, { style: 'percent', minimumFractionDigits: 1, maximumFractionDigits: 1 });
const formatPercentTick = v => (v / 100).toLocaleString(LOCALE, { style: 'percent', maximumFractionDigits: 2 });
// Elementsymbol in üblicher Schreibweise (CA → Ca); die Farbe kommt aus den Atom-Variablen.
const elementSymbol = el => el ? el[0].toUpperCase() + el.slice(1).toLowerCase() : '?';
const elementColor = el => `var(--atom-${(el || '').toUpperCase()}, var(--atom-default))`;

const AMINO_ACID_NAMES = {
    ALA: gettext('Alanin'), ARG: gettext('Arginin'), ASN: gettext('Asparagin'), ASP: gettext('Asparaginsäure'),
    CYS: gettext('Cystein'), GLN: gettext('Glutamin'), GLU: gettext('Glutaminsäure'), GLY: gettext('Glycin'),
    HIS: gettext('Histidin'), ILE: gettext('Isoleucin'), LEU: gettext('Leucin'), LYS: gettext('Lysin'),
    MET: gettext('Methionin'), PHE: gettext('Phenylalanin'), PRO: gettext('Prolin'), SER: gettext('Serin'),
    THR: gettext('Threonin'), TRP: gettext('Tryptophan'), TYR: gettext('Tyrosin'), VAL: gettext('Valin'),
};

// Name einer Aminosäure; modifizierte (z. B. MSE) über ihre Ausgangs-Aminosäure, falls bekannt.
function aminoAcidLabel(name, parent) {
    if (AMINO_ACID_NAMES[name]) return AMINO_ACID_NAMES[name];
    if (parent && AMINO_ACID_NAMES[parent]) return interpolate(gettext('%s (modifiziert)'), [AMINO_ACID_NAMES[parent]]);
    return name;
}
// Ausgangs-Aminosäure je Residuenname (nur modifizierte)
const parentsByName = residues => new Map(residues.filter(r => r.parent).map(r => [r.name, r.parent]));

const BOND_TYPE_LABELS = {
    backbone: gettext('Backbone'), sidechain: gettext('Seitenkette'), peptide: gettext('Peptidbindung'),
    disulfide: gettext('Disulfidbrücke'),
};

// Tabelle "Anzahl und Anteil", standardmäßig absteigend nach Anzahl sortiert.
// rows: [{ key, label, sub, count, color, dot }]; der Anteil bezieht sich auf `total`.
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

function renderStats(graphData) {
    if (d3.select('#stat-atoms').empty()) return;
    const { residues, nodes } = graphData;

    // Nur Standard-Aminosäuren zählen, keine Liganden, Ionen oder Wasser (HETATM-Residuen).
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

// Abschnitt "je Kette": Matrix (Zeilen × Ketten) mit Anteilen bezogen auf die jeweilige
// Kette sowie die Verteilung der gezählten Objekte auf die Ketten (absolut und relativ).
// items: gezählte Objekte (Atome bzw. Aminosäuren) mit `chain`; keyOf: Zeilenschlüssel;
// rowLabel: Kopf der ersten Spalte; renderLabel(th, key): Zeilenbeschriftung.
function renderChainBreakdown({ matrix, share, shareTotal, items, keyOf, rowLabel, renderLabel }) {
    const table = d3.select(matrix);
    if (table.empty()) return;

    const chains = [...new Set(items.map(n => n.chain))]; // Reihenfolge wie in der Datei
    const chainTotals = d3.rollup(items, v => v.length, n => n.chain);
    const byChain = d3.rollup(items, v => v.length, n => n.chain, keyOf);
    const keys = d3.rollups(items, v => v.length, keyOf)
        .sort((a, b) => d3.descending(a[1], b[1]) || d3.ascending(a[0], b[0]))
        .map(([key]) => key);

    // Tabelle 1: Anteile beziehen sich auf die jeweilige Kette.
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

    // Tabelle 2: Verteilung auf die Ketten.
    renderShareTable(share,
        chains.map(c => ({ key: c, label: chainLabel(c), count: chainTotals.get(c), color: 'var(--bs-primary)' })),
        items.length,
        { sort: false });
    d3.select(shareTotal).text(formatInt(items.length));
}

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
    // Nur Ketten mit Aminosäuren erscheinen (reine Liganden-Ketten hätten keine Anteile).
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

// ===== Legende =====
// Elemente des angezeigten Models: zuerst die Hauptelemente der Proteine, dann alphabetisch.
const LEGEND_ORDER = ['C', 'N', 'O', 'S', 'H', 'P'];

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

// ===== Referenzwerte (Engh & Huber 1991) =====
// e.ref = [Mittelwert, Standardabweichung, Gruppe] vom Server oder null (keine Referenz).
const Z_WARN = 2;
const Z_OUTLIER = 4;
const zScore = e => e.ref ? (e.length - e.ref[0]) / e.ref[1] : null;
const zClass = z => z === null ? 'none' : Math.abs(z) > Z_OUTLIER ? 'outlier' : Math.abs(z) > Z_WARN ? 'warn' : 'ok';
const formatZ = z => z.toLocaleString(LOCALE, { minimumFractionDigits: 1, maximumFractionDigits: 1, signDisplay: 'exceptZero' });
const formatRmsz = v => v === null ? '–' : v.toLocaleString(LOCALE, { minimumFractionDigits: 2, maximumFractionDigits: 2 });

// Bezeichnung einer Referenzgruppe, z. B. "N:CA:GLY" -> ["N–CA", "Gly"], "C:N:PRO" -> ["C–N", "vor Pro"].
function referenceGroupLabel(group) {
    const [a, b, cls] = group.split(':');
    const qualifier = {
        GLY: 'Gly', ALA: 'Ala', IVT: 'Ile, Thr, Val', all: gettext('alle'), other: gettext('übrige'),
        PRO: a === 'C' && b === 'N' ? gettext('vor Pro') : 'Pro',
    }[cls] || cls;
    return [`${a}–${b}`, qualifier];
}

function referenceText(e) {
    if (!e.ref) return gettext('Ohne Referenzwert');
    return interpolate(gettext('Referenz %(ref)s ± %(sigma)s Å · Z = %(z)s'), {
        ref: formatLength(e.ref[0]), sigma: formatLength(e.ref[1]), z: formatZ(zScore(e)),
    }, true);
}

// Bewertung aller Bindungen mit Referenz: Gesamtwerte, je Gruppe und Ausreißer (nach |Z| absteigend).
function referenceStats(graphData) {
    const nodeById = new Map(graphData.nodes.map(n => [n.id, n]));
    const bonds = graphData.edges
        .filter(e => nodeById.has(e.source) && nodeById.has(e.target))
        .map(e => ({ ...e, atoms: [nodeById.get(e.source), nodeById.get(e.target)], z: zScore(e) }));
    const rated = bonds.filter(b => b.z !== null);
    const rmsz = items => items.length ? Math.sqrt(d3.mean(items, b => b.z * b.z)) : null;
    const groups = d3.groups(rated, b => b.ref[2])
        .map(([group, items]) => ({
            group, ref: items[0].ref[0], sigma: items[0].ref[1], count: items.length,
            mean: d3.mean(items, b => b.length), rmsz: rmsz(items),
            outliers: items.filter(b => Math.abs(b.z) > Z_OUTLIER).length,
        }))
        .sort((a, b) => d3.descending(a.count, b.count) || d3.ascending(a.group, b.group));
    return {
        total: bonds.length, rated: rated.length, rmsz: rmsz(rated), groups,
        outliers: rated.filter(b => Math.abs(b.z) > Z_OUTLIER).sort((a, b) => d3.descending(Math.abs(a.z), Math.abs(b.z))),
    };
}

const REFERENCE_OUTLIER_ROWS = 50; // Im Tab und im PDF; alle stehen im JSON-Export

function renderReferenceStats(graphData) {
    const facts = d3.select('#stat-ref-facts');
    if (facts.empty()) return;
    const s = referenceStats(graphData);

    facts.selectAll('div')
        .data([
            [gettext('Bewertete Bindungen'), interpolate(gettext('%(rated)s von %(total)s (%(share)s)'), {
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

    const groups = d3.select('#stat-ref-groups');
    groups.selectAll('*').remove();
    const columns = [
        [gettext('Bindung')], [gettext('Klasse')], [gettext('Referenz (Å)'), true], [gettext('Anzahl'), true],
        [gettext('Durchschnitt (Å)'), true], [gettext('Abw. (Å)'), true], [gettext('RMSZ'), true], [gettext('Ausreißer'), true],
    ];
    groups.append('thead').append('tr').selectAll('th').data(columns).join('th')
        .attr('scope', 'col').attr('class', c => c[1] ? 'text-end' : null).text(c => c[0]);
    const rows = groups.append('tbody').selectAll('tr').data(s.groups).join('tr');
    rows.append('th').attr('scope', 'row').text(g => referenceGroupLabel(g.group)[0]);
    rows.append('td').text(g => referenceGroupLabel(g.group)[1]);
    rows.append('td').attr('class', 'text-end').text(g => `${formatLength(g.ref)} ± ${formatLength(g.sigma)}`);
    rows.append('td').attr('class', 'text-end').text(g => formatInt(g.count));
    rows.append('td').attr('class', 'text-end').text(g => formatLength(g.mean));
    rows.append('td').attr('class', 'text-end').text(g => formatDelta(g.mean - g.ref));
    rows.append('td').attr('class', g => `text-end z-${zClass(g.rmsz)}`).text(g => formatRmsz(g.rmsz));
    rows.append('td').attr('class', g => `text-end${g.outliers ? ' z-outlier' : ''}`).text(g => formatInt(g.outliers));

    // Ausreißer einzeln, nach |Z| absteigend
    d3.select('#stat-ref-outliers-title').text(interpolate(gettext('Ausreißer (|Z| > %s)'), [Z_OUTLIER]));
    const list = d3.select('#stat-ref-outliers');
    list.selectAll('*').remove();
    d3.select('#stat-ref-outliers-block').property('hidden', !s.outliers.length);
    if (!s.outliers.length) return;
    list.append('thead').append('tr').selectAll('th')
        .data([[gettext('Atome')], [gettext('Kette')], [gettext('Pos.')], [gettext('Aminosäure')],
               [gettext('Länge (Å)'), true], [gettext('Referenz (Å)'), true], ['Z', true]])
        .join('th').attr('scope', 'col').attr('class', c => c[1] ? 'text-end' : null).text(c => c[0]);
    const body = list.append('tbody');
    s.outliers.slice(0, REFERENCE_OUTLIER_ROWS).forEach(b => {
        const d = describeBond(b.atoms);
        const tr = body.append('tr');
        [d.atoms, d.chain, d.pos, d.residues].forEach(v => tr.append('td').text(v));
        tr.append('td').attr('class', 'text-end').text(formatLength(b.length));
        tr.append('td').attr('class', 'text-end').text(`${formatLength(b.ref[0])} ± ${formatLength(b.ref[1])}`);
        tr.append('td').attr('class', 'text-end z-outlier').text(formatZ(b.z));
    });
    if (s.outliers.length > REFERENCE_OUTLIER_ROWS) {
        body.append('tr').append('td').attr('colspan', 7).attr('class', 'small text-body-secondary')
            .text(interpolate(gettext('… und %s weitere (vollständig im JSON-Export)'), [formatInt(s.outliers.length - REFERENCE_OUTLIER_ROWS)]));
    }
}

// ===== Winkel: Bindungswinkel, Torsionswinkel φ/ψ und ω =====
// Berechnet im Browser aus den Atomkoordinaten des angezeigten Models (bzw. Durchschnitts).
const CIS_MAX_OMEGA = 30;     // |ω| < 30°: cis-Peptidbindung
const TRANS_MIN_OMEGA = 150;  // |ω| ≥ 150°: trans; dazwischen verdreht
const BACKBONE_ORDER = ['N', 'CA', 'C', 'O', 'OXT'];
const DEG = 180 / Math.PI;

const vec = (a, b) => [b.x - a.x, b.y - a.y, b.z - a.z];
const dot = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
const cross = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];

// Winkel a–b–c in Grad (b ist das mittlere Atom).
function bondAngle(a, b, c) {
    const u = vec(b, a), v = vec(b, c);
    const cos = dot(u, v) / Math.sqrt(dot(u, u) * dot(v, v));
    return Math.acos(Math.max(-1, Math.min(1, cos))) * DEG;
}

// Torsionswinkel p0–p1–p2–p3 in Grad (−180 … 180, IUPAC-Konvention).
function dihedral(p0, p1, p2, p3) {
    const b0 = vec(p1, p0), b2 = vec(p2, p3);
    let b1 = vec(p1, p2);
    const len = Math.sqrt(dot(b1, b1));
    b1 = b1.map(v => v / len);
    const v = b0.map((c, i) => c - dot(b0, b1) * b1[i]);
    const w = b2.map((c, i) => c - dot(b2, b1) * b1[i]);
    return Math.atan2(dot(cross(b1, v), w), dot(v, w)) * DEG;
}

// Auf 0,1° gerundet; `|| 0` macht aus −0 eine 0 (sonst "−0,0°")
const roundAngle = v => Math.round(v * 10) / 10 || 0;
const formatAngle = v => `${roundAngle(v).toLocaleString(LOCALE, { minimumFractionDigits: 1, maximumFractionDigits: 1 })}°`;

// Klasse im Ramachandran-Plot; modifizierte Aminosäuren zählen wie ihre Ausgangs-Aminosäure.
function ramaClass(res) {
    const name = res.parent || res.name;
    return name === 'GLY' ? 'gly' : name === 'PRO' ? 'pro' : 'general';
}
const RAMA_CLASSES = ['general', 'gly', 'pro'];
const ramaClassLabel = cls => ({ general: gettext('Übrige Aminosäuren'), gly: gettext('Glycin'), pro: gettext('Prolin') })[cls];

// φ/ψ je Aminosäure, ω je Peptidbindung und alle Bindungswinkel. Torsionen über eine
// Peptidbindung gibt es nur, wo der Server eine Bindung erkannt hat (keine Kettenlücke).
function analyzeAngles(graphData) {
    const nodes = new Map(graphData.nodes.map(n => [n.id, n]));
    const peptides = new Set(graphData.edges.filter(e => e.type === 'peptide').map(edgeKey));
    const atom = (res, name) => nodes.get(`${res.chain}_${res.pos}_${name}`);
    const linked = (a, b) => {
        const c = atom(a, 'C'), n = atom(b, 'N');
        return Boolean(c && n && peptides.has(`${c.id}|${n.id}`));
    };

    const residues = graphData.residues;
    const ramachandran = [];
    const omegas = [];
    residues.forEach((res, i) => {
        const prev = residues[i - 1], next = residues[i + 1];
        const n = atom(res, 'N'), ca = atom(res, 'CA'), c = atom(res, 'C');
        if (!n || !ca || !c) return;
        const phi = prev && linked(prev, res) ? dihedral(atom(prev, 'C'), n, ca, c) : null;
        const psi = next && linked(res, next) ? dihedral(n, ca, c, atom(next, 'N')) : null;
        if (phi !== null && psi !== null) ramachandran.push({ res, phi, psi, cls: ramaClass(res) });
        if (next && linked(res, next) && atom(next, 'CA')) {
            const omega = dihedral(ca, c, atom(next, 'N'), atom(next, 'CA'));
            const abs = Math.abs(omega);
            const kind = abs < CIS_MAX_OMEGA ? 'cis' : abs < TRANS_MIN_OMEGA ? 'twisted' : 'trans';
            omegas.push({ prev: res, res: next, omega, kind });
        }
    });

    // Bindungswinkel: je Atom alle Paare seiner Bindungspartner. Reihenfolge der Endatome:
    // früheres Residuum zuerst (C–N–CA über die Peptidbindung), im Residuum nach Rückgrat
    // (N, CA, C, O); der Partner einer Disulfidbrücke steht immer am Ende (CB–SG–SG).
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
                const angle = { atoms: [a, center, c], type: `${a.label}–${center.label}–${c.label}`, value: bondAngle(a, center, c) };
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

// Ramachandran-Plot in `svg` (Tab und Bild-Export). Farben kommen über CSS-Klassen
// (.rama-general, .rama-gly, .rama-pro). Liefert die Auswahl der Punkte.
const RAMA_SIZE = 480;
const RAMA_MAJOR_STEP = 60; // Hauptgitter und Achsenbeschriftung (°)
const RAMA_MINOR_STEP = 20; // Hilfsgitter (°)
const RAMA_MARGIN = { top: 12, right: 14, bottom: 46, left: 56 };
function drawRamachandran(svg, points, size) {
    const m = RAMA_MARGIN;
    const x = d3.scaleLinear().domain([-180, 180]).range([m.left, size - m.right]);
    const y = d3.scaleLinear().domain([-180, 180]).range([size - m.bottom, m.top]);
    const ticks = d3.range(-180, 181, RAMA_MAJOR_STEP);
    const degrees = d => `${d}°`;

    // Hilfsgitter (feiner als das Hauptgitter), unter allem anderen
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

    // Glycin und Prolin über den übrigen zeichnen, damit sie sichtbar bleiben
    const order = cls => RAMA_CLASSES.indexOf(cls);
    return svg.append('g').selectAll('circle')
        .data([...points].sort((a, b) => order(a.cls) - order(b.cls)))
        .join('circle')
        .attr('class', p => `rama-point rama-${p.cls}`)
        .attr('cx', p => x(p.phi)).attr('cy', p => y(p.psi)).attr('r', 3.2);
}

function renderAngles(canvas, graphData) {
    const container = d3.select('#angle-rama');
    if (container.empty()) return;
    const angles = analyzeAngles(graphData);

    // Ramachandran-Plot (skaliert über viewBox mit der Breite)
    container.selectAll('*').remove();
    const svg = container.append('svg')
        .attr('viewBox', `0 0 ${RAMA_SIZE} ${RAMA_SIZE}`)
        .attr('role', 'img')
        .attr('aria-label', interpolate(gettext('Ramachandran-Plot mit %s Aminosäuren'), [formatInt(angles.ramachandran.length)]));
    canvas.ramaPoints = drawRamachandran(svg, angles.ramachandran, RAMA_SIZE)
        .on('click', (event, p) => jumpToResidue(canvas, p.res));
    canvas.ramaPoints.append('title').text(p =>
        `${p.res.name} ${p.res.pos} (${chainLabel(p.res.chain)})\nφ = ${formatAngle(p.phi)}, ψ = ${formatAngle(p.psi)}`);

    const counts = d3.rollup(angles.ramachandran, v => v.length, p => p.cls);
    d3.select('#angle-rama-legend').selectAll('div')
        .data(RAMA_CLASSES)
        .join(enter => {
            // Punkt an der ersten Textzeile ausrichten (lange Bezeichnungen brechen in der Leiste um)
            const div = enter.append('div').attr('class', 'd-flex align-items-start gap-2');
            div.append('span').attr('class', 'legend-dot rama-legend-dot');
            div.append('span').attr('class', 'rama-legend-label');
            return div;
        })
        .call(div => div.select('.legend-dot').style('background', cls => `var(--rama-${cls})`))
        .call(div => div.select('.rama-legend-label').text(cls => `${ramaClassLabel(cls)} (${formatInt(counts.get(cls) || 0)})`));

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
    // Liste nur, wenn es Peptidbindungen gibt, die nicht trans sind
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

    renderAngleGroups(angles);
    updateSelection(canvas);
}

const formatDegrees = v => v === null || v === undefined ? '–'
    : roundAngle(v).toLocaleString(LOCALE, { minimumFractionDigits: 1, maximumFractionDigits: 1 });
const formatDegreesDelta = v => v.toLocaleString(LOCALE,
    { minimumFractionDigits: 1, maximumFractionDigits: 1, signDisplay: 'exceptZero' });

// Atome eines Bindungswinkels für die Tabellen (wie describeBond, mit drei Atomen).
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

// Bindungswinkel je Typ; jede Zeile klappt die einzelnen Winkel auf (wie bei den Bindungen).
const expandedAngleGroups = new Set();
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

// Einzelne Winkel eines Typs in Dateireihenfolge; kleinster und größter sind markiert.
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

// Winkel für den JSON-/PDF-Export
function anglesExport(graphData) {
    const angles = analyzeAngles(graphData);
    const round1 = v => v === null ? null : roundTo(v, 1);
    const residue = r => ({ residue: r.name, position: r.pos, chain: r.chain.trim() });
    return {
        unit: '°',
        ramachandran: angles.ramachandran.map(p => ({ ...residue(p.res), class: p.cls, phi: round1(p.phi), psi: round1(p.psi) })),
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
    };
}

// Ramachandran-Plot als eigenständige SVG-Datei (feste Farben, weißer Hintergrund).
const RAMA_EXPORT_SIZE = { width: 800, height: 800 };
// Als Funktion, da HIST_EXPORT_STYLE weiter unten in der Datei definiert ist
const ramaExportStyle = () => HIST_EXPORT_STYLE + `
    .rama-frame { stroke: #adb5bd; }
    .rama-minor line { stroke: rgba(0, 0, 0, 0.045); }
    .rama-zero { stroke: #adb5bd; stroke-dasharray: 3 3; }
    .rama-point { stroke: #fff; stroke-width: 0.8; }
    .rama-general { fill: #0d6efd; }
    .rama-gly { fill: #f08c00; }
    .rama-pro { fill: #2f9e44; }
`;
function ramachandranSvgText(graphData) {
    const { width, height } = RAMA_EXPORT_SIZE;
    const svg = d3.create('svg')
        .attr('width', width).attr('height', height).attr('viewBox', `0 0 ${width} ${height}`);
    svg.append('style').text(ramaExportStyle());
    svg.append('rect').attr('width', width).attr('height', height).attr('fill', '#fff');
    drawRamachandran(svg, analyzeAngles(graphData).ramachandran, width);
    return new XMLSerializer().serializeToString(svg.node());
}

async function ramachandranBlob(graphData, format) {
    const svgText = ramachandranSvgText(graphData);
    if (format === 'svg') return new Blob([svgText], { type: 'image/svg+xml' });
    return rasterizeSvg(svgText, `image/${format}`, RAMA_EXPORT_SIZE);
}

// Ramachandran-Plot als Vektorgrafik im PDF (Achsen, Raster, Punkte je Klasse).
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
    doc.setFontSize(8);
    doc.text(pdfText('φ'), (plotL + plotR) / 2, top + size - 1, { align: 'center' });
    doc.text(pdfText('ψ'), left + 1, (plotT + plotB) / 2, { angle: 90 });
}

// ===== Abschnitt "Bindungen" =====
const formatLength = v => v === null || v === undefined ? '–'
    : v.toLocaleString(LOCALE, { minimumFractionDigits: 3, maximumFractionDigits: 3 });
const formatDelta = v => v.toLocaleString(LOCALE,
    { minimumFractionDigits: 3, maximumFractionDigits: 3, signDisplay: 'exceptZero' });
const chainId = c => c.trim() || '–';
const bondKey = (source, target) => `${source}|${target}`;

// Beschriftung einer Bindung aus ihren beiden Atomen (Graph-Knoten oder Atome aus dem Ensemble).
// Liegen beide Atome in derselben Aminosäure, werden Position und Aminosäure nur einmal genannt.
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

// Vergleichswert mit der Abweichung der Bindung davon (Länge − Wert); `title` erklärt die Abweichung.
function appendComparison(td, value, length, title) {
    td.append('span').text(formatLength(value));
    if (value === null || value === undefined) return;
    td.append('span').attr('class', 'delta')
        .attr('title', title)
        .text(`(${formatDelta(length - value)})`);
}

// Model einer Extrem-Bindung aus dem Ensemble; Dateien ohne MODEL-Records enthalten genau
// ein (implizites) Model.
const ensembleModelLabel = b => b && interpolate(gettext('Model %s'), [b.model ?? 1]);

// Kürzeste und längste Bindung über alle Models und im ausgewählten Model, jeweils
// verglichen mit Durchschnitt und Median dieser Bindung über alle Models.
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

// Aufgeklappte Gruppen bleiben beim Model-Wechsel offen.
const expandedBondGroups = new Set();
const CHEVRON_SVG = '<svg xmlns="http://www.w3.org/2000/svg" width="12" height="12" fill="currentColor" viewBox="0 0 16 16" aria-hidden="true"><path fill-rule="evenodd" d="M4.646 1.646a.5.5 0 0 1 .708 0l6 6a.5.5 0 0 1 0 .708l-6 6a.5.5 0 0 1-.708-.708L10.293 8 4.646 2.354a.5.5 0 0 1 0-.708"/></svg>';

// Kennzahlen einer Menge von Bindungen.
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

// Alle Bindungen des ausgewählten Models, gruppiert nach den beteiligten Atomen (z. B. C–N,
// CB–CG). Eine aufgeklappte Gruppe listet jede einzelne Bindung.
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

// Einzelne Bindungen einer Gruppe in Dateireihenfolge; kürzeste und längste sind markiert.
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

// ===== Histogramm der Bindungslängen =====
const HIST_HEIGHT = 260;
const HIST_MARGIN = { top: 16, right: 16, bottom: 44, left: 64 };
const histogram = { edges: [], width: 0 };
const HIST_BIN_WIDTH = 0.01; // Å

// Längenklassen zu 0,01 Å mit Anteil an allen Bindungen; leere Klassen innerhalb des
// Wertebereichs bleiben erhalten. Die Längen haben drei Nachkommastellen, gerechnet wird
// daher ganzzahlig in Hundertstel-Å.
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

// Balken mit oben abgerundeten Ecken, unten bündig auf der Grundlinie.
function topRoundedBar(x, y, w, h, r) {
    r = Math.min(r, w / 2, h);
    return `M${x},${y + h}V${y + r}Q${x},${y} ${x + r},${y}H${x + w - r}Q${x + w},${y} ${x + w},${y + r}V${y + h}Z`;
}

// Anteil der Bindungen je Längenklasse (0,01 Å) an allen Bindungen des ausgewählten Models.
// Ohne Argument wird mit den zuletzt übergebenen Bindungen neu gezeichnet (z. B. nach Größenänderung).
function renderBondHistogram(edges = histogram.edges) {
    histogram.edges = edges;
    const container = d3.select('#stat-bond-histogram');
    if (container.empty()) return;
    const width = container.node().clientWidth;
    // Tab ausgeblendet: Breite zurücksetzen, damit der ResizeObserver beim Anzeigen neu zeichnet.
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

    // Hover: Trefferfläche ist die ganze Spalte, nicht nur der Balken.
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

// Achsen, Beschriftungen und Balken des Histogramms in `svg` (im Tab und für den Bild-Export).
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

    // Balken mit 2px Abstand; auch sehr kleine Anteile bleiben als 1px-Balken sichtbar.
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

// ===== Export des Histogramms (SVG / PNG / JPEG) =====
// Feste helle Farben wie im PDF, damit das Bild unabhängig vom Dark Mode gleich aussieht.
const HIST_EXPORT_SIZE = { width: 1200, height: 480 };
const HIST_EXPORT_SCALE = 2; // Pixeldichte der Rasterbilder
const HIST_EXPORT_STYLE = `
    text { font-family: Roboto, Arial, sans-serif; }
    .hist-bar { fill: #0d6efd; }
    .hist-axis text { fill: #6c757d; font-size: 13px; font-variant-numeric: tabular-nums; }
    .hist-axis .domain, .hist-axis .tick line { stroke: #dee2e6; }
    .hist-axis .tick line.hist-grid { stroke: rgba(0, 0, 0, 0.1); }
    .hist-label { fill: #6c757d; font-size: 14px; }
`;

// Eigenständige SVG-Datei des Histogramms (weißer Hintergrund, eingebettete Styles).
function histogramSvgText(edges) {
    const { width, height } = HIST_EXPORT_SIZE;
    const svg = d3.create('svg')
        .attr('width', width).attr('height', height).attr('viewBox', `0 0 ${width} ${height}`);
    svg.append('style').text(HIST_EXPORT_STYLE);
    svg.append('rect').attr('width', width).attr('height', height).attr('fill', '#fff');
    drawHistogram(svg, histogramBins(edges), width, height);
    return new XMLSerializer().serializeToString(svg.node());
}

// SVG-Text als Rasterbild (`image/png` oder `image/jpeg`).
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

async function histogramBlob(edges, format) {
    const svgText = histogramSvgText(edges);
    if (format === 'svg') return new Blob([svgText], { type: 'image/svg+xml' });
    return rasterizeSvg(svgText, `image/${format}`);
}

// Disulfidbrücken mit ihren Atomen (Graph-Knoten), in Reihenfolge der Datei.
function disulfideBonds(graphData) {
    const nodeById = new Map(graphData.nodes.map(n => [n.id, n]));
    return graphData.edges
        .filter(e => e.type === 'disulfide' && nodeById.has(e.source) && nodeById.has(e.target))
        .map(e => ({ ...e, atoms: [nodeById.get(e.source), nodeById.get(e.target)] }));
}

// Je Cystein (residueKey) die Partner seiner Disulfidbrücken, z. B. "CYS 127" oder "CYS 214 (Kette L)".
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

// Peptidbindungen: existente und nicht-existente (Kettenlücken), bezogen auf die theoretisch möglichen.
function renderPeptideBonds(graphData) {
    const { possible, formed, breaks } = analyzePeptideBonds(graphData);
    d3.select('#stat-peptide-possible').text(formatInt(possible));
    renderShareTable('#stat-peptide-bonds', [
        { key: 'formed', label: gettext('Existent'), count: formed, color: 'var(--bs-primary)' },
        { key: 'missing', label: gettext('Nicht-existent'), sub: gettext('Kettenlücken'),
          count: breaks.size, color: 'var(--bs-danger)' },
    ], possible || 1, { sort: false });
}

function renderBondStats(graphData, ensemble, modelLabel) {
    if (d3.select('#stat-bond-extremes').empty() || !ensemble) return;
    renderReferenceStats(graphData);
    renderPeptideBonds(graphData);
    renderDisulfides(graphData);
    renderBondExtremes(graphData, ensemble, modelLabel);
    renderBondHistogram(graphData.edges);
    renderBondGroups(graphData);
}

// ===== Export der Statistiken (JSON / PDF) =====
const roundTo = (v, digits) => v === null || v === undefined ? null : +v.toFixed(digits);
const sharePercent = (count, total) => roundTo(100 * count / total, 2);

// Atom einer Bindung im Export (Graph-Knoten und Ensemble-Atome haben dieselben Felder).
const exportAtom = n => ({
    name: n.label, serial: n.serial, residue: n.residue,
    position: `${n.residue_seq}${n.i_code || ''}`, chain: n.chain.trim(),
});
// Rückweg für die Beschriftung im PDF (describeBond erwartet Felder wie die Graph-Knoten).
const nodeLikeAtom = a => ({ label: a.name, serial: a.serial, residue: a.residue, residue_seq: a.position, chain: a.chain });

// Strukturmetadaten aus dem Header (vom Server in die Seite eingebettet).
const structureInfo = () => {
    const el = document.getElementById('structure-info');
    return el ? JSON.parse(el.textContent) : null;
};

// Alle Werte des Statistik-Tabs als ein Objekt; Grundlage für JSON- und PDF-Export.
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
            // Abweichung von Referenzwerten (Rückgrat und CA–CB), Z = (Länge − Referenz) / Standardabweichung
            reference: {
                source: info?.reference_source ?? null,
                outlier_threshold: Z_OUTLIER,
                rated: reference.rated,
                rmsz: roundTo(reference.rmsz, 3),
                outlier_count: reference.outliers.length,
                groups: reference.groups.map(g => ({
                    bond: referenceGroupLabel(g.group)[0], class: referenceGroupLabel(g.group)[1],
                    reference: g.ref, sigma: g.sigma, count: g.count, mean: roundTo(g.mean, 3),
                    rmsz: roundTo(g.rmsz, 3), outliers: g.outliers,
                })),
                outliers: reference.outliers.map(b => ({
                    length: b.length, reference: b.ref[0], sigma: b.ref[1], z: roundTo(b.z, 2),
                    atoms: b.atoms.map(exportAtom),
                })),
            },
            // Zwischen aufeinanderfolgenden Aminosäuren einer Kette; missing = Kettenlücken
            // via: "ssbond" (SSBOND-Eintrag) oder "distance" (SG–SG-Abstand)
            disulfides: disulfideBonds(graphData).map(b => ({ length: b.length, via: b.via, atoms: b.atoms.map(exportAtom) })),
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
            // Klassen sind halboffen: from <= Länge < to.
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

const exportBaseName = (filename, model, kind = gettext('statistiken')) =>
    `${filename.replace(/\.[^.]+$/, '')}_${kind}_${model.toLowerCase().replace(/[^\p{L}\p{N}]+/gu, '-')}`;

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

// Skripte (jsPDF) erst beim ersten PDF-Export nachladen.
const loadedScripts = new Map();
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

// Die PDF-Standardschriften kennen nur WinAnsi-Zeichen (kein Δ, kein Minuszeichen U+2212).
// Auf Modulebene, da auch drawPdfHistogram sie braucht.
const pdfText = v => String(v).replace(/φ/g, 'phi').replace(/ψ/g, 'psi').replace(/ω/g, 'omega').replace(/Δ/g, 'Delta').replace(/−/g, '-');

// PDF aus den Exportdaten: dieselben Abschnitte und Tabellen wie im Tab; die Liste der
// einzelnen Bindungen je Atompaar steht nur im JSON (bei großen Proteinen Tausende Zeilen).
function buildStatsPdf(data) {
    const doc = new window.jspdf.jsPDF({ unit: 'mm', format: 'a4' });
    const M = 15;
    const pageW = doc.internal.pageSize.getWidth();
    const pageH = doc.internal.pageSize.getHeight();
    const contentW = pageW - 2 * M;
    const INK = [33, 37, 41], MUTED = [108, 117, 125], LINE = [222, 226, 230], PRIMARY = [13, 110, 253];
    const LENGTH_MIN = [11, 114, 133], LENGTH_MAX = [217, 72, 15]; // wie --length-min/max-color
    let y = M;

    const ensureSpace = h => { if (y + h > pageH - M) { doc.addPage(); y = M; } };
    const sectionTitle = text => {
        ensureSpace(50);
        y += 2;
        doc.setFont('helvetica', 'bold').setFontSize(12).setTextColor(...INK).text(pdfText(text), M, y);
        doc.setDrawColor(...LINE).setLineWidth(0.3).line(M, y + 2, pageW - M, y + 2);
        y += 8;
    };
    // Unterüberschrift nur dort, wo noch Platz für Kopf und einige Tabellenzeilen ist.
    const subTitle = (text, minSpace = 40) => {
        ensureSpace(minSpace);
        doc.setFont('helvetica', 'bold').setFontSize(9.5).setTextColor(...INK).text(pdfText(text), M, y);
        y += 2.5;
    };
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

    // Kopf
    doc.setFont('helvetica', 'bold').setFontSize(16).setTextColor(...INK).text(pdfText(interpolate(gettext('Statistiken – %s'), [data.file])), M, y + 4);
    y += 10;
    const created = new Date(data.created).toLocaleString(LOCALE, { dateStyle: 'medium', timeStyle: 'short' });
    doc.setFont('helvetica', 'normal').setFontSize(9).setTextColor(...MUTED)
        .text(pdfText(interpolate(gettext('%(model)s · erstellt am %(date)s'), { model: data.model, date: created }, true)), M, y);
    y += 10;

    // Allgemein
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

    // Elemente / Aminosäuren je Kette
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

    // Bindungen
    const b = data.bonds;
    sectionTitle(gettext('Bindungen'));
    const ref = b.reference;
    subTitle(gettext('Abweichung von Referenzwerten'));
    table({
        head: [[gettext('Bindung'), gettext('Klasse'), num(gettext('Referenz (Å)')), num(gettext('Anzahl')),
                num(gettext('Durchschnitt (Å)')), num(gettext('Abw. (Å)')), num(gettext('RMSZ')), num(gettext('Ausreißer'))]],
        body: ref.groups.map(grp => [grp.bond, grp.class, `${formatLength(grp.reference)} ± ${formatLength(grp.sigma)}`,
            formatInt(grp.count), formatLength(grp.mean), formatDelta(grp.mean - grp.reference), formatRmsz(grp.rmsz), formatInt(grp.outliers)]),
        foot: [[gettext('Gesamt'), '', '', num(formatInt(ref.rated)), '', '', num(formatRmsz(ref.rmsz)), num(formatInt(ref.outlier_count))]],
        columnStyles: right([2, 3, 4, 5, 6, 7]),
    });
    if (ref.outliers.length) {
        subTitle(interpolate(gettext('Ausreißer (|Z| > %s)'), [ref.outlier_threshold]));
        table({
            head: [[gettext('Atome'), gettext('Kette'), gettext('Pos.'), gettext('Aminosäure'), num(gettext('Länge (Å)')), num(gettext('Referenz (Å)')), num('Z')]],
            body: ref.outliers.slice(0, REFERENCE_OUTLIER_ROWS).map(o => {
                const d = describeBond(o.atoms.map(nodeLikeAtom));
                return [d.atoms, d.chain, d.pos, d.residues, formatLength(o.length),
                        `${formatLength(o.reference)} ± ${formatLength(o.sigma)}`, formatZ(o.z)];
            }),
            columnStyles: right([4, 5, 6]),
        });
    }
    ensureSpace(8);
    doc.setFont('helvetica', 'normal').setFontSize(7.5).setTextColor(...MUTED)
        .text(pdfText(interpolate(gettext('Referenzwerte: %s'), [ref.source || ''])), M, y - 3);
    y += 4;
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
    ensureSpace(8);
    doc.setFont('helvetica', 'normal').setFontSize(7.5).setTextColor(...MUTED)
        .text(pdfText(gettext('Die einzelnen Bindungen je Atompaar sind im JSON-Export enthalten.')), M, y - 3);

    // Winkel
    const ang = data.angles;
    const RAMA = { general: [13, 110, 253], gly: [240, 140, 0], pro: [47, 158, 68] };
    sectionTitle(gettext('Winkel'));
    subTitle(gettext('Ramachandran-Plot'), 100);
    const ramaSize = 95;
    ensureSpace(ramaSize + 4);
    drawPdfRamachandran(doc, ang.ramachandran, M, y + 2, ramaSize, { MUTED, LINE, ...RAMA });
    RAMA_CLASSES.forEach((cls, i) => {
        const ly = y + 8 + i * 6;
        doc.setFillColor(...RAMA[cls]).circle(M + ramaSize + 10, ly - 1, 1, 'F');
        doc.setFont('helvetica', 'normal').setFontSize(8).setTextColor(...INK).text(pdfText(
            `${ramaClassLabel(cls)} (${formatInt(ang.ramachandran.filter(p => p.class === cls).length)})`), M + ramaSize + 13, ly);
    });
    y += ramaSize + 9;
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

    // Fußzeile auf jeder Seite
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

// Histogramm als Vektorgrafik im PDF (gleiche Gestaltung wie im Tab).
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
(function () {
    const dataEl = document.getElementById('graph-data');
    if (!dataEl) return;
    let graphData = JSON.parse(dataEl.textContent);
    const ensembleEl = document.getElementById('bond-ensemble');
    let bondEnsemble = ensembleEl ? JSON.parse(ensembleEl.textContent) : null;

    let activeStyle  = localStorage.getItem('bondStyle') || 'gradient';
    let activeLayout = localStorage.getItem('layout')    || 'schematic';

    const canvas = setupCanvas();

    function redraw(animate = true) {
        render(canvas, graphData, activeStyle, activeLayout, animate);
    }

    const selectsWrap = document.getElementById('structure-selects');
    const modelSelect = document.getElementById('model-select');
    const altlocSelect = document.getElementById('altloc-select');
    // Deaktivierte Felder (nur ein Model / keine Alternativen) bieten keine Auswahl.
    const selectable = [modelSelect, altlocSelect].filter(s => s && !s.disabled);

    // Name des angezeigten Models, bei gewählter Konformation mit deren Kennung.
    function modelLabel() {
        const model = modelSelect
            ? modelSelect.selectedOptions[0].textContent.trim()
            : interpolate(gettext('Model %s'), [bondEnsemble?.models[0] ?? 1]);
        const altloc = altlocSelect && !altlocSelect.disabled ? altlocSelect.value : '';
        return altloc ? `${model} · ${interpolate(gettext('Konformation %s'), [altloc])}` : model;
    }
    // Model der aktuell geladenen Daten (das Auswahlfeld wechselt schon vor dem Laden).
    let loadedModelLabel = modelLabel();

    // Lädt Graph-Daten für die gewählte Kombination aus Model und Konformation; nach einem
    // Konformationswechsel auch die Bindungsstatistik über alle Models.
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
        pdfButton.addEventListener('click', async () => {
            pdfButton.disabled = true;
            try {
                await loadScript(jspdfSrc);
                await loadScript(autotableSrc);
                const data = exportData();
                buildStatsPdf(data).save(`${exportBaseName(filename, data.model)}.pdf`);
            } catch (err) {
                console.error('PDF konnte nicht erstellt werden:', err);
            } finally {
                pdfButton.disabled = false;
            }
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

    function updateButtons() {
        document.querySelectorAll('.style-option[data-style]').forEach(btn => {
            btn.classList.toggle('active', btn.dataset.style === activeStyle);
        });
        const referenceLegend = document.getElementById('legend-reference');
        if (referenceLegend) referenceLegend.hidden = activeStyle !== 'reference';
        document.querySelectorAll('.style-option[data-layout]').forEach(btn => {
            btn.classList.toggle('active', btn.dataset.layout === activeLayout);
        });
    }

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

    updateButtons();
    redraw(false);
    renderDataTable(canvas, graphData.nodes);
    renderStats(graphData);
    renderLegend(graphData.nodes);
    renderBondStats(graphData, bondEnsemble, modelLabel());
    renderAngles(canvas, graphData);

    const anglesExportBar = document.getElementById('angles-export');
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

    // Tabellen-Scrollen aktualisiert die Hervorhebung in der Sequenzleiste (einmal pro Frame).
    let tableScrollPending = false;
    canvas.tableScroll.addEventListener('scroll', () => {
        if (tableScrollPending) return;
        tableScrollPending = true;
        requestAnimationFrame(() => {
            tableScrollPending = false;
            if (canvas.state.activeTab === 'data') updateTableView(canvas);
        });
    });

    // Die Sequenzleiste gehört nur zu den Tabs Graph und Daten.
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
                // Im ausgeblendeten Tab übersprungenes Histogramm nachholen.
                const histogramEl = document.getElementById('stat-bond-histogram');
                if (histogramEl && histogramEl.clientWidth !== histogram.width) renderBondHistogram();
            }
        });
    });

    const observer = new MutationObserver(() => redraw(false));
    observer.observe(document.documentElement, { attributes: true, attributeFilter: ['data-bs-theme'] });

    // Bandbreite und Zoom-Grenzen hängen von der Containergröße ab.
    // Im ausgeblendeten Tab (Größe 0) nicht neu zeichnen, sonst verrutschen die Zoom-Grenzen.
    new ResizeObserver(() => {
        if (canvas.svg.node().clientWidth > 0) redraw(false);
    }).observe(canvas.svg.node());

    // Histogramm an die Breite anpassen; im ausgeblendeten Tab erst beim Anzeigen zeichnen.
    const histogramEl = document.getElementById('stat-bond-histogram');
    if (histogramEl) {
        new ResizeObserver(() => {
            if (histogramEl.clientWidth !== histogram.width) renderBondHistogram();
        }).observe(histogramEl);
    }
})();
