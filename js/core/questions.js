// Question-shaped entry points. The home page offered a search box and one example gene,
// which asks the reader to already know what they are looking for. These are the questions
// the shipped layers can actually answer, each one a real Browse query.
//
// TWO RULES, both learned the hard way in this repo:
//
//   1. Counts are computed IN THE BROWSER from the facet index, never written into the
//      manifest. A count baked into an aggregate that a build stage owns goes stale the
//      moment a patcher runs without it -- TODO 25 documents exactly that failure, and the
//      home page has already shipped a withdrawn number once (0.3.11, the 1,725 LSG pool).
//      Computing here means the number cannot disagree with the Browse page it links to,
//      because it IS the Browse page's own filter function.
//
//   2. Every cutoff appears in the question text. "Expressed" and "tissue-specific" are not
//      absolute properties, and a question that hides its threshold invites reading it as
//      one. If you cannot state the cutoff in the sentence, the question is wrong.
//
// `params` are literally Browse's own FILTER ids, so browseHref(params) round-trips into a
// working, shareable, deep-linked view with the widgets pre-filled.
export const QUESTIONS = [
  {
    id: 'young-and-duplicated',
    q: 'Which tandem-array genes are xylem-specific?',
    sub: 'in a curated tandem array, peak expression in xylem, τ ≥ 0.8 (no age or Ks cut)',
    params: { tandem: 'yes', top: 'xylem', tau_min: '0.8' },
  },
  {
    id: 'lsg-with-protein',
    q: 'Which lineage-specific candidates have a peptide match in public Populus mass spectra?',
    sub: 'phylostratum ≥ 18, a peptide in at least one of 5 public datasets; mostly single, shared or cross-species (P. tomentosa, P. × canescens) peptides',
    params: { ps_min: '18', ms: 'yes' },
  },
  {
    id: 'hemizygous',
    q: 'Which genes exist on only one haplotype?',
    sub: 'Strict PAV: no syntelog and no protein-level match on the other haplotype',
    params: { rel: 'hemizygous_strict_PAV' },
  },
  {
    id: 'three-layers',
    q: 'Which tandem-array genes are also called for both whole-genome events?',
    sub: 'a curated tandem array plus sWGD and aWGT anchors, under the primary call',
    params: { tandem: 'yes', wgd9: 'both' },
  },
  {
    id: 'never-expressed',
    q: 'Which genes have a median WT-control TPM of 1 or less in every tissue?',
    sub: 'per-tissue medians over WT-control samples only, 10 tissues; single samples may exceed it',
    params: { xin: '' },     // handled specially below: xin is a tissue mask, not a value
    test: (r) => !r.xin || r.xin.length === 0,
  },
  {
    id: 'fast-evolving',
    q: 'Which genes look fast-evolving against P. trichocarpa?',
    sub: 'ω ≥ 1 compared with the reciprocal-best ortholog, a flag to inspect, not a selection claim',
    params: { w_min: '1.0' },
  },
  {
    id: 'private',
    q: 'Which genes are classed private (no ortholog in the panel)?',
    sub: 'synteny class across the 12-genome panel, which omits the parent P. alba',
    params: { cls: 'private_no_ortholog' },
  },
  {
    id: 'disputed-wgd',
    q: 'Where do the two whole-genome-duplication calls disagree?',
    // NOT { wgd: 'salicoid', wgd9: 'gamma' }. That was the first draft and it reported
    // 7,574, because a gene carrying BOTH events under BOTH calls matches it while agreeing
    // about everything. The disagreement is its own fact and needs its own column.
    sub: 'the exploratory topology re-check compared with the primary call, gene by gene',
    params: { wgd_disputed: 'yes' },
  },
];

/** Count matches for a question over the facet rows, using Browse's OWN filter functions so
 *  the number here and the number on the linked page cannot drift apart. */
export function countQuestion(q, rows, FILTERS) {
  if (q.test) return rows.reduce((n, r) => n + (q.test(r) ? 1 : 0), 0);
  const active = FILTERS.filter((f) => q.params[f.id] !== undefined && q.params[f.id] !== '');
  if (!active.length) return null;
  let n = 0;
  for (const r of rows) {
    let ok = true;
    for (const f of active) {
      if (!f.test(r, q.params[f.id])) { ok = false; break; }
    }
    if (ok) n += 1;
  }
  return n;
}
