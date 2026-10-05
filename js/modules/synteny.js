// Module: synteny / syntelog.
// Reference implementation of the module contract -- copy this shape for new modules.
import { el, pct, fmt } from '../core/dom.js';
import { browseHref } from '../core/browse.js';

export const CLASSES = {
  core_syntenic:         { label: 'Core syntenic',    dot: 'core',        blurb: 'Syntenic placement, orthogroup spans ≥9 of 12 genomes.' },
  dispensable_syntenic:  { label: 'Dispensable',      dot: 'dispensable', blurb: 'Syntenic placement, orthogroup spans 2–8 genomes.' },
  non_syntenic_ortholog: { label: 'Non-syntenic',     dot: 'nonsyn',      blurb: 'Has orthologs but no syntenic placement.' },
  private_no_ortholog:   { label: 'Private',          dot: 'private',     blurb: 'No ortholog found elsewhere in the panel (orthogroup spans 1 genome).' },
  tandem_array_member:   { label: 'Tandem array',     dot: 'tandem',      blurb: 'Array-only placement; a non-representative member of a tandem array.' },
};

export function classBadge(cls) {
  const c = CLASSES[cls];
  if (!c) return el('span', { class: 'badge' }, cls || '–');
  // Swatch ALWAYS ships with its text label -- three light-mode hues sit under 3:1
  // contrast, so identity may never be carried by color alone.
  return el('span', { class: 'badge', title: c.blurb }, el('i', { class: 'dot ' + c.dot }), c.label);
}

function geneCard(g) {
  const s = g.syn;
  if (!s) return null;
  const kv = el('dl', { class: 'kv' });
  const add = (k, ...v) => { kv.append(el('dt', {}, k), el('dd', {}, ...v)); };

  add('Class', classBadge(s.cls));
  add('Orthogroup', s.og ? el('span', { class: 'mono' }, s.og) : el('span', { class: 'muted' }, 'unassigned by OrthoFinder'));
  add('Genomes in orthogroup', s.n == null ? '–' : `${s.n} of 12`);

  const f = s.flags || {};
  const on = Object.entries({ 'syntenic (PASS)': f.pass, 'non-syntenic ortholog': f.ns, 'tandem array': f.ar })
    .filter(([, v]) => v).map(([k]) => k);
  add('Placement flags', on.length ? on.join(' · ') : el('span', { class: 'muted' }, 'none'));

  if (s.stable === false) {
    kv.append(el('dt', {}, 'Stability'),
      el('dd', {}, el('span', { class: 'badge' }, '⚠ sensitivity-unstable'),
        el('div', { class: 'muted', style: 'font-size:12px;margin-top:3px' },
          'This call differs between the default and --ultra-sensitive DIAMOND runs.')));
  }

  const card = el('div', { class: 'card' },
    el('h2', {}, 'Synteny', el('span', { class: 'tag' }, 'GENESPACE')), kv);
  if (s.cls === 'private_no_ortholog') {
    card.append(el('div', { class: 'warn' },
      'Private calls are panel-dependent. P. alba, a parent species of this hybrid, is not in ' +
      'the panel, and private genes are far more likely than core genes to have a homolog in P. alba.'));
  }
  return card;
}

function overview(man) {
  const counts = man.summary.synteny_class || {};
  const total = Object.values(counts).reduce((a, b) => a + b, 0);
  const order = Object.keys(CLASSES).filter((k) => counts[k]);
  const max = Math.max(...order.map((k) => counts[k]));

  // Each class row is now the way into Browse for that class. Browse's own filter id for
  // this is `cls` (this module's `syn_class` filter is the gene-page one), so the link uses
  // the id Browse actually reads.
  const rows = order.map((k) => el('tr', {},
    el('td', {}, el('a', { href: browseHref({ cls: k }),
                           title: `Browse the ${fmt(counts[k])} ${CLASSES[k].label} genes`,
                           style: 'text-decoration:none' }, classBadge(k))),
    el('td', { class: 'num' }, fmt(counts[k])),
    el('td', { class: 'num' }, pct(counts[k], total)),
    el('td', { style: 'width:44%' },
      el('div', { class: 'bar', style: `width:${(counts[k] / max) * 100}%;background:var(--c-${
        { core_syntenic: 'core', dispensable_syntenic: 'disp', non_syntenic_ortholog: 'nonsyn',
          private_no_ortholog: 'priv', tandem_array_member: 'tandem' }[k]})` })),
  ));

  return el('div', {},
    el('h1', {}, 'Synteny'),
    el('p', { class: 'sub' },
      'GENESPACE 1.3.1 over a 12-genome panel (10 syntenic + O. sativa / P. patens as ' +
      'orthogroup-only outgroups). Classes follow the original five-way vocabulary. The panel '
      + 'lacks P. alba, one parent of this hybrid, so "Private" means no ortholog in the panel, '
      + 'not absent from every other genome.'),
    el('div', { class: 'stat-row' },
      el('div', { class: 'stat' }, el('div', { class: 'n' }, fmt(total)), el('div', { class: 'l' }, 'genes classified')),
      el('div', { class: 'stat' }, el('div', { class: 'n' }, fmt(man.summary.unstable)), el('div', { class: 'l' }, 'DIAMOND-sensitivity unstable')),
    ),
    el('div', { class: 'card' },
      el('h2', {}, 'Class distribution',
        el('span', { class: 'tag' }, 'click a class to browse it')),
      el('div', { class: 'scroll-x' }, el('table', { class: 'data' },
        el('thead', {}, el('tr', {}, el('th', {}, 'Class'), el('th', { class: 'num' }, 'Genes'),
          el('th', { class: 'num' }, 'Share'), el('th', {}, ''))),
        el('tbody', {}, ...rows))),
      el('div', { class: 'muted', style: 'margin-top:10px;font-size:12px' },
        Object.entries(CLASSES).map(([, c]) => c.label + ', ' + c.blurb).join('  ·  ')),
    ),
  );
}

export default {
  id: 'synteny',
  label: 'Synteny',
  geneCard,
  overview,
  filters: [
    { id: 'syn_class', label: 'Synteny class', options: Object.keys(CLASSES),
      test: (g, v) => g.syn && g.syn.cls === v },
  ],
};
