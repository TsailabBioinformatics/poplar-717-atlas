// Route: #/karyotype -- navigate by the ancestral genome instead of the modern one.
//
// The whole duplication half of this atlas is about descent: which ancestral chromosome a gene
// came from, which gamma copy, which salicoid copy. This page is that frame, in three views:
//
//   1. Painted chromosomes: the 19 modern chromosomes of each haplotype, to scale, every gene
//      painted by its ancestral chromosome (AEK1-7) -- or, within one chosen AEK, by gamma copy
//      or by gamma x salicoid copy. HAP1 and HAP2 are separate figures, never merged.
//   2. The 42-cell grid (7 AEK x 3 gamma x 2 salicoid, plus "no copy assigned"), each cell a
//      Browse query.
//   3. Per-AEK breakdown: for each copy slot of one ancestral chromosome, which modern
//      chromosomes carry it, per haplotype.
//
// NO NEW DATA. Counts are reduced in the browser from the `aek` facet column (chromosome*100 +
// gamma*10 + salicoid). Positions come from data/index/pos.json (gene model start and end)
// (data/locus/index/<hap>/<chr>.json), joined to the facet row by gene ID. Precomputing any of
// this would create an aggregate a build stage owns and can silently leave stale (TODO 25).
//
// EXPLORATORY, and labelled so. Copy numbers are Wang et al.'s completeness RANKS within one
// ancestral chromosome, not subgenomes, and "no genes" in a cell means no collinear anchor
// carries that slot in this haplotype -- never a proven loss.
import { el, fmt, clear } from '../core/dom.js';
import { loadFacets, facetRow, loadPositions as loadAllPositions } from '../core/data.js';
import { browseHref } from '../core/browse.js';

const CHRS = [1, 2, 3, 4, 5, 6, 7];
const GAMMA = [1, 2, 3];
const SALICOID = [1, 2];
const HAPS = ['hap1', 'hap2'];
const HAPNAME = { hap1: 'HAP1 (TreH)', hap2: 'HAP2 (AlbH)' };
const MODERN = Array.from({ length: 19 }, (_, i) => 'Chr' + String(i + 1).padStart(2, '0'));
const SVGNS = 'http://www.w3.org/2000/svg';

const dec = (v) => ({ c: Math.floor(v / 100), g: Math.floor(v / 10) % 10, s: v % 10 });

function svg(tag, attrs = {}, ...kids) {
  const n = document.createElementNS(SVGNS, tag);
  for (const [k, v] of Object.entries(attrs)) if (v != null) n.setAttribute(k, v);
  for (const c of kids.flat()) if (c != null) n.append(c instanceof Node ? c : document.createTextNode(String(c)));
  return n;
}

/** Reduce the facet index to the 42 cells, per haplotype and pooled. */
function grid(rows) {
  const cells = new Map();
  const bump = (key, r) => {
    if (!cells.has(key)) cells.set(key, { n: 0, hap1: 0, hap2: 0 });
    const c = cells.get(key);
    c.n += 1;
    c[r.hap] += 1;
  };
  let placed = 0;
  let noCopy = 0;
  let unplaced = 0;
  for (const r of rows) {
    const v = r.aek;
    if (v == null || v < 0) { unplaced += 1; continue; }
    const { c, g, s } = dec(v);
    if (!c) { unplaced += 1; continue; }
    placed += 1;
    if (!g || !s) { noCopy += 1; bump(`${c}|0|0`, r); continue; }
    bump(`${c}|${g}|${s}`, r);
  }
  return { cells, placed, noCopy, unplaced };
}

/** Browse filter value for an AEK, optionally narrowed to a gamma copy and a salicoid copy. */
function aekValue(c, g, s) {
  if (!g) return `AEK${c}`;
  if (!s) return `AEK${c}_gamma_copy_${g}`;
  return `AEK${c}_gamma_copy_${g}_salicoid_copy_${s}`;
}
const cellHref = (c, g, s) => browseHref({ aek: aekValue(c, g, s) });

function cellBox(c, g, s, cell, max) {
  const n = cell ? cell.n : 0;
  // Tint by count on a log scale, capped at 50% so the ink numbers stay readable on top.
  const t = n ? 0.10 + 0.40 * (Math.log10(n + 1) / Math.log10(max + 1)) : 0;
  const bg = n ? `background:color-mix(in srgb, var(--accent) ${Math.round(t * 100)}%, transparent)` : '';
  const inner = el('div', { class: 'kp-cell' + (n ? '' : ' empty'), style: bg },
    el('div', { class: 'kp-total' }, n ? fmt(n) : '–'),
    n ? el('div', { class: 'kp-haps' }, `HAP1 ${fmt(cell.hap1)}`, el('br'), `HAP2 ${fmt(cell.hap2)}`) : null);
  if (!n) {
    return el('div', { title: 'No collinear anchor in either haplotype carries this slot. '
      + 'That is not evidence the gene was lost.' }, inner);
  }
  return el('a', {
    class: 'kp-celllink', href: cellHref(c, g, s), 'data-kp-cell': `${c}|${g}|${s}`, 'data-n': n,
    title: `AEK${c}, aWGT copy ${g || 'none'}, sWGD copy ${s || 'none'}: ${fmt(n)} genes `
      + `(${fmt(cell.hap1)} HAP1 + ${fmt(cell.hap2)} HAP2). Click to browse them.`,
  }, inner);
}

/* ---------------- painted chromosomes ---------------- */

/** The colour category of one gene under the current mode, or null if it is not painted. */
function category(v, mode, focus) {
  if (v == null || v < 0) return null;
  const { c, g, s } = dec(v);
  if (!c) return null;
  if (mode === 'aek') return `a${c}`;
  if (c !== focus) return 'other';
  if (mode === 'gamma') return `g${g}`;
  return g && s ? `g${g}s${s}` : 'g0s0';
}

function categoriesFor(mode, focus) {
  if (mode === 'aek') {
    return CHRS.map((c) => ({ key: `a${c}`, label: `AEK${c}`, aek: aekValue(c) }));
  }
  if (mode === 'gamma') {
    return [...GAMMA.map((g) => ({ key: `g${g}`, label: `AEK${focus} aWGT copy ${g}`, aek: aekValue(focus, g) })),
      { key: 'g0', label: `AEK${focus}, no copy assigned`, aek: null }];
  }
  return [...GAMMA.flatMap((g) => SALICOID.map((s) => ({
    key: `g${g}s${s}`, label: `AEK${focus} gamma ${g} · salicoid ${s}`, aek: aekValue(focus, g, s) }))),
  { key: 'g0s0', label: `AEK${focus}, no copy assigned`, aek: null }];
}

/** Merge consecutive genes (sorted by start) of the same category into painted segments. */
function segments(genes, mode, focus) {
  const out = [];
  let cur = null;
  for (const gn of genes) {
    const k = category(gn.aek, mode, focus);
    if (k == null) { cur = null; continue; }
    if (cur && cur.k === k) { cur.e = Math.max(cur.e, gn.e); cur.n += 1; continue; }
    cur = { k, s: gn.s, e: gn.e, n: 1, v: gn.aek };
    out.push(cur);
  }
  return out;
}

// Positions come from one aligned file (data/index/pos.json, ~290 KB gzipped) rather than the
// 38 per-chromosome locus index files (~1.5 MB) this page used to fetch for the same numbers.
async function loadPositions(rows, onProgress) {
  const pos = await loadAllPositions();
  const set = new Set(MODERN);
  const out = { hap1: {}, hap2: {} };
  for (const hap of HAPS) for (const chr of MODERN) out[hap][chr] = { genes: [], len: 0 };
  for (let i = 0; i < rows.length; i++) {
    const r = rows[i];
    if (!set.has(r.chr) || !out[r.hap]) continue;
    const c = out[r.hap][r.chr];
    c.genes.push({ s: pos.start[i], e: pos.end[i], aek: r.aek });
    if (pos.end[i] > c.len) c.len = pos.end[i];
  }
  for (const hap of HAPS) for (const chr of MODERN) out[hap][chr].genes.sort((a, b) => a.s - b.s);
  onProgress(1, 1);
  return out;
}

/* ---------------- the view ---------------- */

export async function karyotypeView() {
  const facets = await loadFacets();
  const rows = facets.rows.map((_, i) => facetRow(facets, i));
  const { cells, placed, noCopy, unplaced } = grid(rows);
  const max = Math.max(...[...cells.values()].map((c) => c.n), 1);

  // count(hap, chr|null, aek value filter) straight from the facets, the same column Browse filters
  const counts = new Map();
  for (const r of rows) {
    if (r.aek < 0 || !dec(r.aek).c) continue;
    const { c, g, s } = dec(r.aek);
    for (const chr of [null, r.chr]) {
      for (const key of [aekValue(c), g ? aekValue(c, g) : null, g && s ? aekValue(c, g, s) : null]) {
        if (!key) continue;
        const k = `${r.hap}|${chr ?? '*'}|${key}`;
        counts.set(k, (counts.get(k) || 0) + 1);
      }
    }
  }
  const count = (hap, chr, aek) => counts.get(`${hap}|${chr ?? '*'}|${aek}`) || 0;

  const state = { mode: 'aek', focus: 3, pos: null, sel: null };

  /* ----- painted chromosomes card ----- */
  const controls = el('div', { class: 'kp-controls' });
  const figs = el('div', { class: 'kp-figs' }, el('p', { class: 'muted', 'data-kp-loading': '' },
    'Loading gene positions for 38 chromosomes…'));
  const detail = el('div', { class: 'kp-detail', 'data-kp-detail': '' });
  const brk = el('div', { class: 'kp-brk' });
  const brkHead = el('h2', {});

  function renderControls() {
    clear(controls);
    const modes = [['aek', 'Ancestral chromosome'], ['gamma', 'aWGT copy'], ['salicoid', 'aWGT × sWGD copy']];
    controls.append(el('span', { class: 'muted', style: 'font-size:13px' }, 'Colour by'),
      ...modes.map(([m, label]) => el('button', {
        type: 'button', class: 'btn' + (state.mode === m ? ' on' : ''), 'data-kp-mode': m,
        'aria-pressed': state.mode === m ? 'true' : 'false',
        onclick: () => { state.mode = m; state.sel = null; renderAll(); },
      }, label)));
    const pick = el('div', { class: 'kp-controls kp-aekpick', style: 'margin:0' },
      el('span', { class: 'muted', style: 'font-size:13px' }, 'Ancestral chromosome in focus'),
      ...CHRS.map((c) => el('button', {
        type: 'button', class: 'btn' + (state.focus === c ? ' on' : ''), 'data-kp-focus': c,
        'aria-pressed': state.focus === c ? 'true' : 'false',
        onclick: () => { state.focus = c; state.sel = null; renderAll(); },
      }, `AEK${c}`)));
    controls.append(pick);
  }

  function legendFor(hap) {
    const cats = categoriesFor(state.mode, state.focus);
    return el('ul', { class: 'kp-legend', 'data-kp-legend-hap': hap },
      ...cats.map((cat) => {
        const n = cat.aek ? count(hap, null, cat.aek)
          : rows.reduce((a, r) => a + (r.hap === hap && r.aek >= 0 && dec(r.aek).c === state.focus
            && !(dec(r.aek).g && (state.mode === 'gamma' || dec(r.aek).s)) ? 1 : 0), 0);
        const sw = el('span', { class: 'kp-sw', style: `background:${swatch(cat.key)}` });
        const body = [sw, el('span', {}, cat.label), el('span', { class: 'kp-n' }, fmt(n))];
        if (!cat.aek) return el('li', {}, el('a', { href: browseHref({ aek: aekValue(state.focus), hap }), title: 'Browse every gene of this ancestral chromosome; filter the no-copy ones in the table' }, ...body));
        return el('li', {}, el('a', {
          href: browseHref({ aek: cat.aek, hap }), 'data-kp-legend': `${hap}|${cat.aek}`, 'data-n': n,
          title: `${HAPNAME[hap]}: ${fmt(n)} genes in ${cat.label}. Click to browse them.`,
        }, ...body));
      }));
  }

  function swatch(key) {
    if (key.startsWith('a')) return `var(--aek${key.slice(1)})`;
    const m = /^g(\d)(?:s(\d))?$/.exec(key);
    if (!m || m[1] === '0') return 'var(--ink-3)';
    return m[2] === '2' ? `color-mix(in srgb, var(--kp-g${m[1]}) 45%, var(--panel))` : `var(--kp-g${m[1]})`;
  }

  function figure(hap) {
    const data = state.pos[hap];
    const maxLen = Math.max(...HAPS.flatMap((h) => MODERN.map((c) => state.pos[h][c].len)));
    const L = 52; const R = 12; const T = 26; const ROW = 22; const BAR = 12; const W = 760;
    const sx = (bp) => L + (bp / maxLen) * (W - L - R);
    const H = T + MODERN.length * ROW + 6;
    const root = svg('svg', { class: 'kp-svg', viewBox: `0 0 ${W} ${H}`, role: 'img',
      'aria-label': `${HAPNAME[hap]}: 19 chromosomes painted by ancestral origin` });
    for (let mb = 0; mb * 1e6 <= maxLen; mb += 10) {
      root.append(svg('line', { class: 'kp-tick', x1: sx(mb * 1e6), x2: sx(mb * 1e6), y1: T - 6, y2: H - 4 }),
        svg('text', { x: sx(mb * 1e6), y: T - 10, 'text-anchor': 'middle' }, `${mb}`));
    }
    root.append(svg('text', { x: W - R, y: 10, 'text-anchor': 'end' }, 'Mb'));
    let painted = 0;
    MODERN.forEach((chr, i) => {
      const y = T + i * ROW + (ROW - BAR) / 2;
      const { genes, len } = data[chr];
      const segs = segments(genes, state.mode, state.focus);
      const g = svg('g', { 'data-kp-row': `${hap}|${chr}` });
      g.append(svg('text', { x: L - 6, y: y + BAR - 2, 'text-anchor': 'end' }, chr),
        svg('rect', { class: 'kp-bg', x: sx(0), y, width: sx(len) - sx(0), height: BAR, rx: 3 }));
      for (const s of segs) {
        painted += s.n;
        g.append(svg('rect', { class: `kp-${s.k}`, x: sx(s.s), y, width: Math.max(1.2, sx(s.e) - sx(s.s)), height: BAR }));
      }
      if (state.sel && state.sel.hap === hap && state.sel.chr === chr) {
        const s = state.sel.seg;
        g.append(svg('rect', { class: 'kp-sel', x: sx(s.s) - 1, y: y - 2, width: Math.max(3, sx(s.e) - sx(s.s)) + 2, height: BAR + 4 }));
      }
      // one hit area per chromosome row; the click picks the nearest painted segment
      const hit = svg('rect', { class: 'kp-hit', x: sx(0), y: T + i * ROW, width: W - L - R, height: ROW,
        'data-kp-hit': `${hap}|${chr}` });
      hit.addEventListener('click', (ev) => {
        const box = root.getBoundingClientRect();
        const x = ((ev.clientX - box.left) / box.width) * W;
        const bp = ((x - L) / (W - L - R)) * maxLen;
        let best = null; let bd = Infinity;
        for (const s of segs) {
          const d = bp < s.s ? s.s - bp : (bp > s.e ? bp - s.e : 0);
          if (d < bd) { bd = d; best = s; }
        }
        if (!best) return;
        state.sel = { hap, chr, seg: best };
        renderFigs();
        showDetail();
      });
      g.append(hit);
      root.append(g);
    });
    const unpainted = rows.filter((r) => r.hap === hap).length - painted;
    return el('div', { class: 'kp-hap', 'data-kp-fig': hap },
      el('h3', {}, HAPNAME[hap]),
      el('div', { class: 'kp-svgbox' }, root),
      el('p', { class: 'muted', style: 'font-size:12px;margin:6px 0 0', 'data-kp-painted': painted },
        `n = ${fmt(painted)} genes painted on 19 chromosomes; ${fmt(unpainted)} ${HAPNAME[hap]} genes are `
        + 'not painted (not placed on the ancestral genome, or on an unplaced scaffold). '
        + 'Click a chromosome to inspect the segment under the pointer.'),
      el('p', { style: 'font-size:12.5px;margin:10px 0 0;color:var(--ink-2)' },
        state.mode === 'aek' ? `Genes per ancestral chromosome, ${HAPNAME[hap]} (whole haplotype, scaffolds included):`
          : `AEK${state.focus} genes by copy, ${HAPNAME[hap]}; other ancestral chromosomes are drawn in pale grey:`),
      legendFor(hap));
  }

  function renderFigs() {
    clear(figs);
    figs.append(...HAPS.map(figure));
  }

  function showDetail() {
    clear(detail);
    const sel = state.sel;
    if (!sel) return;
    const { c, g, s } = dec(sel.seg.v);
    const g2 = state.mode === 'aek' ? 0 : g;
    const s2 = state.mode === 'salicoid' ? s : 0;
    const what = state.mode === 'aek' || c === state.focus
      ? `AEK${c}${g2 ? ` aWGT copy ${g2}` : ''}${s2 ? ` sWGD copy ${s2}` : ''}${state.mode !== 'aek' && !g ? ' (no copy assigned)' : ''}`
      : 'genes from ancestral chromosomes other than the one in focus';
    const q = state.mode === 'aek' || (c === state.focus && g) ? aekValue(c, g2 || 0, g2 ? s2 : 0) : aekValue(c);
    const n = count(sel.hap, sel.chr, q);
    const mb = (bp) => (bp / 1e6).toFixed(2);
    detail.append(
      el('p', { style: 'margin:0 0 8px;font-size:14px' },
        el('strong', {}, `${HAPNAME[sel.hap]} ${sel.chr}, ${mb(sel.seg.s)}–${mb(sel.seg.e)} Mb`),
        `: a run of ${fmt(sel.seg.n)} consecutive genes painted ${what}.`),
      el('div', { class: 'kp-controls', style: 'margin:0' },
        el('a', { class: 'btn', href: `#/map/${sel.hap}/${sel.chr}`, 'data-kp-map': '' },
          `Open ${sel.chr} (${sel.hap.toUpperCase()}) in the genome map`),
        el('a', { class: 'btn primary', href: browseHref({ aek: q, chr: sel.chr, hap: sel.hap }),
          'data-kp-browse': q, 'data-n': n },
        `Browse all ${fmt(n)} ${q.replace(/_/g, ' ')} genes on ${sel.hap.toUpperCase()} ${sel.chr}`)));
  }

  /* ----- per-AEK breakdown ----- */
  function renderBreakdown() {
    clear(brk);
    const c = state.focus;
    brkHead.replaceChildren(`Where AEK${c}'s copies sit today`,
      el('span', { class: 'tag' }, 'per haplotype, by modern chromosome'));
    for (const g of GAMMA) {
      for (const s of SALICOID) {
        const q = aekValue(c, g, s);
        const bars = HAPS.map((hap) => {
          const parts = MODERN.map((chr) => ({ chr, n: count(hap, chr, q) })).filter((p) => p.n);
          const tot = count(hap, null, q);
          const onChr = parts.reduce((a, p) => a + p.n, 0);
          parts.sort((a, b) => b.n - a.n);
          const top = parts[0];
          const bar = el('div', { class: 'kp-bar', 'data-kp-bar': `${hap}|${q}` },
            ...parts.map((p, k) => el('a', {
              class: `${hap === 'hap1' ? 'h1' : 'h2'}${k % 2 ? 'b' : 'a'}`,
              style: `flex:${p.n} 0 0`,
              href: browseHref({ aek: q, chr: p.chr, hap }), 'data-n': p.n,
              title: `${HAPNAME[hap]} ${p.chr}: ${fmt(p.n)} of ${fmt(tot)} genes (${Math.round((100 * p.n) / tot)}%)`,
            }, p.n / (onChr || 1) >= 0.14 ? p.chr.replace('Chr', 'Chr ') : '')));
          const say = tot
            ? `${hap.toUpperCase()}: n = ${fmt(tot)}; mostly ${top.chr} (${fmt(top.n)}, ${Math.round((100 * top.n) / tot)}%)`
              + `${parts.length > 1 ? `, then ${parts.slice(1, 3).map((p) => `${p.chr} (${fmt(p.n)})`).join(', ')}` : ''}`
              + `${tot > onChr ? `; ${fmt(tot - onChr)} on unplaced scaffolds` : ''}`
            : `${hap.toUpperCase()}: no genes carry this slot (absence of an anchor, not a proven loss)`;
          return [el('div', { class: 'kp-barrow' }, el('span', { class: 'kp-hl' }, hap.toUpperCase()), bar),
            el('div', { class: 'kp-say' }, say)];
        });
        brk.append(el('div', { class: 'kp-slot' },
          el('h4', {}, `AEK${c} · aWGT copy ${g} · sWGD copy ${s}`), ...bars.flat()));
      }
    }
  }

  function renderAll() {
    renderControls();
    if (state.pos) renderFigs();
    showDetail();
    renderBreakdown();
  }

  /* ----- grid ----- */
  const header = el('tr', {}, el('th', { style: 'text-align:left' }, ''),
    ...CHRS.map((c) => el('th', { style: 'text-align:center' },
      el('span', { class: 'kp-sw', style: `background:var(--aek${c});vertical-align:-2px;margin-right:4px` }), `AEK${c}`)));
  const body = [];
  for (const g of GAMMA) {
    for (const s of SALICOID) {
      body.push(el('tr', {},
        el('th', { style: 'text-align:left;font-weight:500;white-space:nowrap;font-size:12.5px' },
          `gamma ${g} · salicoid ${s}`),
        ...CHRS.map((c) => el('td', {}, cellBox(c, g, s, cells.get(`${c}|${g}|${s}`), max)))));
    }
  }
  body.push(el('tr', {},
    el('th', { style: 'text-align:left;font-weight:500;font-size:12.5px;color:var(--ink-2)' },
      'no copy assigned'),
    ...CHRS.map((c) => el('td', {}, cellBox(c, 0, 0, cells.get(`${c}|0|0`), max)))));

  renderAll();
  const view = el('div', {},
    el('h1', {}, 'Ancestral karyotype', el('span', { class: 'badge', style: 'margin-left:10px;font-size:12px;padding:1px 7px' }, 'exploratory')),
    el('p', { class: 'sub' },
      'Which ancestral chromosome each stretch of the genome is collinear with. Seven ancestral '
      + 'eudicot chromosomes (AEK1 to AEK7), each in three aWGT copies, each in two sWGD copies. '
      + 'This is a positional label for the region, not the descent of each gene: a gene younger '
      + 'than these events, including most lineage-specific genes, still carries the label of the '
      + 'segment it sits in, and some genes inherit it from their nearest anchors. '
      + 'Below: the modern chromosomes painted by ancestry, the 42-cell table, and where each '
      + 'copy landed. Every count is a link to the genes behind it.'),

    el('div', { class: 'stat-row' },
      el('div', { class: 'stat' }, el('div', { class: 'n' }, fmt(placed)), el('div', { class: 'l' }, 'genes placed on the ancestral genome')),
      el('div', { class: 'stat' }, el('div', { class: 'n' }, fmt(noCopy)), el('div', { class: 'l' }, 'placed, but no copy assigned')),
      el('div', { class: 'stat' }, el('div', { class: 'n' }, fmt(unplaced)), el('div', { class: 'l' }, 'not placed')),
    ),

    el('div', { class: 'card', 'data-kp-painted-card': '' },
      el('h2', {}, 'Modern chromosomes, painted by ancestry', el('span', { class: 'tag' }, 'HAP1 and HAP2 separately')),
      el('p', { style: 'margin:0 0 8px;font-size:13.5px;color:var(--ink-2)' },
        'Each bar is one modern chromosome drawn to scale; each coloured stretch is a run of '
        + 'consecutive genes placed in segments collinear with the same ancestral chromosome. Switch to a copy view '
        + 'to see, within one ancestral chromosome, which aWGT or sWGD copy each stretch is.'),
      controls, figs, detail),

    el('div', { class: 'card' },
      el('h2', {}, 'The 42 cells', el('span', { class: 'tag' }, 'both haplotypes, split below')),
      el('p', { style: 'margin:0 0 10px;font-size:13.5px;color:var(--ink-2)' },
        'A cell counts the genes (HAP1 and HAP2 together, then each haplotype on its own line) '
        + 'whose ancestral chromosome is the column and whose aWGT and sWGD copy ranks are '
        + 'the row. Click a cell to browse exactly those genes.'),
      el('div', { class: 'scroll-x' },
        el('table', { class: 'data kp-grid', style: 'width:100%' },
          el('thead', {}, header), el('tbody', {}, ...body))),
      el('p', { class: 'muted', style: 'font-size:11.5px;margin:12px 0 0' },
        'Shading is the gene count on a log scale, not confidence. A dash means no collinear '
        + 'anchor in either haplotype carries that slot, that is an absence of evidence, '
        + 'not evidence the gene was lost. The HAP1 and HAP2 lines under each count show where '
        + 'the haplotypes disagree without opening anything.')),

    el('div', { class: 'card', 'data-kp-breakdown': '' },
      brkHead,
      el('p', { style: 'margin:0 0 10px;font-size:13.5px;color:var(--ink-2)' },
        'For the ancestral chromosome in focus (pick it above), each copy slot as a bar per '
        + 'haplotype, split by the modern chromosome its genes sit on now. Each piece is a link.'),
      brk),

    el('div', { class: 'card' },
      el('h2', {}, 'What these numbers are, and are not'),
      el('p', { style: 'margin:0;font-size:13.5px;color:var(--ink-2)' },
        'aWGT copy (1 to 3) and sWGD copy (1 to 2) are completeness ',
        el('strong', {}, 'ranks'),
        ' from Wang et al.’s karyotype projection, ranked within one ancestral chromosome. ',
        el('strong', {}, 'They are not subgenomes'),
        ', and copy 1 of AEK3 has no relationship to copy 1 of AEK5. A copy names the ancestral '
        + 'copy a gene\u2019s segment is collinear with, not whether its duplicate survived; that '
        + 'is the whole-genome-duplication call on each gene page.'),
      el('p', { class: 'muted', style: 'font-size:12px;margin:10px 0 0' },
        'Exploratory. Ancestral chromosomes come from collinearity between the 717 '
        + 'genome and the ancestral karyotype genes of Wang et al. 2022, BMC Biology 20:216; '
        + 'copy names come from collinear blocks shared with P. trichocarpa v3.1. Painted '
        + 'positions are gene start and end coordinates of the v5.1 gene models.')));

  loadPositions(rows, (d, n) => {
    const l = figs.querySelector('[data-kp-loading]');
    if (l) l.textContent = `Loading gene positions… ${d} of ${n} chromosomes`;
  }).then((pos) => { state.pos = pos; renderFigs(); })
    .catch((e) => { clear(figs).append(el('p', { class: 'muted' }, `Gene positions failed to load (${e.message}).`)); });
  return view;
}

export default { id: 'karyotype', label: 'Karyotype' };
