// Shared helpers for the array and pair pages: tile lookups, the matching region on the
// other haplotype, the SVG helper for the pair page's gene models, and the tissue order.
// The locus drawing itself is core/locusbrowser.js.
//
// Tile records (data/locus/tile/<hap>/<chr>/<Mb>.json) carry, per gene: g id, s/e coordinates,
// d strand (+1/-1), al allele id, y synteny class, and f = [[offset, length, isCDS], ...]
// relative to s in genomic orientation. f is the only exon structure in the release.
import { loadLocusTile, getGene } from './data.js';

const NS = 'http://www.w3.org/2000/svg';
export function svgEl(tag, attrs = {}, ...kids) {
  const n = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) if (v != null) n.setAttribute(k, v);
  for (const c of kids.flat()) if (c != null) n.append(c instanceof Node ? c : document.createTextNode(String(c)));
  return n;
}

/** Tile genes overlapping [from, to] on one chromosome, sorted by start. Missing tiles
 *  (chromosome ends, scaffolds without tiles) are skipped, not errors. */
export async function genesIn(hap, chr, from, to) {
  const a = Math.max(0, Math.floor(from / 1e6));
  const b = Math.floor(to / 1e6);
  const out = new Map();
  for (let m = a; m <= b; m += 1) {
    try {
      const t = await loadLocusTile(hap, chr, m);
      for (const r of t || []) if (r.e >= from && r.s <= to) out.set(r.g, r);
    } catch { /* no tile for this Mb */ }
  }
  return [...out.values()].sort((x, y) => x.s - y.s || (x.g < y.g ? -1 : 1));
}

/** One tile record for a gene (with its exon features), or null. */
export async function tileRecord(g) {
  const rows = await genesIn(g.hap, g.chr, g.start, g.end);
  return rows.find((r) => r.g === g.id) || null;
}

/** The region on the OTHER haplotype that matches a set of genes. Uses the genes' own alleles
 *  first; if none has one, the alleles of up to `flank` flanking genes on each side. Returns
 *  {hap, chr, from, to, via: 'members'|'flanks', anchors: [[ownId, otherId]]} or null. */
export async function matchingRegion(own, ownRows, focusIds, flank = 12, keepPerSide = Infinity) {
  const focus = new Set(focusIds);
  let anchors = ownRows.filter((r) => focus.has(r.g) && r.al).map((r) => [r.g, r.al]);
  let via = 'members';
  if (!anchors.length) {
    via = 'flanks';
    const idx = ownRows.map((r, i) => (focus.has(r.g) ? i : -1)).filter((i) => i >= 0);
    const lo = Math.min(...idx);
    const hi = Math.max(...idx);
    // the nearest `keepPerSide` flanks with an allele on each side
    const left = ownRows.filter((r, i) => r.al && i < lo && i >= lo - flank).slice(-keepPerSide);
    const right = ownRows.filter((r, i) => r.al && i > hi && i <= hi + flank).slice(0, keepPerSide);
    anchors = [...left, ...right].map((r) => [r.g, r.al]);
  }
  if (!anchors.length) return null;
  const others = (await Promise.all(anchors.map(([, id]) => getGene(id).catch(() => null)))).filter(Boolean);
  if (!others.length) return null;
  // Majority chromosome: a stray allele call on another chromosome must not stretch the window.
  const count = new Map();
  others.forEach((o) => count.set(o.chr, (count.get(o.chr) || 0) + 1));
  const chr = [...count.entries()].sort((x, y) => y[1] - x[1] || (x[0] < y[0] ? -1 : 1))[0][0];
  const on = others.filter((o) => o.chr === chr);
  const from = Math.min(...on.map((o) => o.start));
  const to = Math.max(...on.map((o) => o.end));
  // Median coordinate distance, to reject one far-away allele from spanning megabases.
  const mids = on.map((o) => (o.start + o.end) / 2).sort((x, y) => x - y);
  const med = mids[Math.floor(mids.length / 2)];
  const keep = on.filter((o) => Math.abs((o.start + o.end) / 2 - med) < 2e6);
  const f2 = Math.min(...keep.map((o) => o.start));
  const t2 = Math.max(...keep.map((o) => o.end));
  return { hap: on[0].hap, chr, from: f2, to: t2, via,
    anchors: anchors.filter(([, id]) => keep.some((o) => o.id === id)), spanAll: [from, to] };
}

/** Tissue order shared by both expression figures: fixed, never alphabetical by value. */
export const TISSUES = ['leaf_young', 'leaf', 'leaf_old', 'bud', 'stem', 'bark', 'xylem', 'root', 'catkin', 'callus'];
export function tissueList(genes) {
  const seen = new Set();
  genes.forEach((g) => Object.keys(((g && g.x) || {}).t || {}).forEach((k) => seen.add(k)));
  return [...TISSUES.filter((t) => seen.has(t)), ...[...seen].filter((t) => !TISSUES.includes(t)).sort()];
}
/** Per-tissue sample n, from the records; returns {tissue: n} and whether records disagree. */
export function tissueN(genes, tissues) {
  const n = {};
  let disagree = false;
  for (const t of tissues) {
    const vals = [...new Set(genes.map((g) => (((g || {}).x || {}).n || {})[t]).filter((v) => v != null))];
    n[t] = vals.length ? vals[0] : null;
    if (vals.length > 1) disagree = true;
  }
  return { n, disagree };
}
