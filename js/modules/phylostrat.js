// Module: phylostratigraphy (genEra). The gene's evolutionary age, its gene family's origin,
// and the caveats the method itself reports.
import { el, fmt, pct } from '../core/dom.js';
import { browseHref } from '../core/browse.js';

const LSG_MIN = 18;

const ABSENSE_VERDICT = {'a_retract_young_call': 'abSENSE retracts the young call', 'b_supports_lineage_specific': 'abSENSE supports lineage-specific', 'b_weak_low_sensitivity': 'abSENSE supports it, weakly (low sensitivity)', 'ancestral_sequence_present_unannotated': 'an ancestral sequence appears present but unannotated', 'ambiguous': 'abSENSE result is ambiguous', 'not_testable': 'not testable with abSENSE'};

// Two QC checks, and only two (Chen, 2026-09-21 evening; ms2-lsg CORRECTIONS C5): abSENSE and
// the synteny window search. A flag never removes a gene. `post` is the chapter's post-QC set,
// the genes with neither flag, used for its group comparisons. It is one line UNDER the flags,
// and is never called "clean" -- a gene can be in that set and still carry a synteny contest, a
// short ORF or a sensitivity note, and until 2026-09-21 a "clean" line hid those on 389 genes.
// The C4 flags are gone: the Idesia rerun is the pool itself now, and WGD-era Ks stopped being
// a QC check (duplication is its own analysis, on the Duplication card).
function qcSummary(qc) {
  const flags = [];
  if (qc.abs) {
    flags.push(el('div', {}, el('span', { class: 'badge' }, 'abSENSE'),
      ' rejects the young call: the absence of older homologs is plausibly a detection failure.'));
  }
  if (qc.ret) {
    flags.push(el('div', {}, el('span', { class: 'badge' }, 'synteny window'),
      ' a tblastn search of the syntenic window finds a strong hit inside an annotated gene of '
      + 'an outgroup, which contradicts the young call.'));
  }
  if (!flags.length) {
    flags.push(el('div', { class: 'muted' }, 'Neither QC check flags this gene.'));
  }
  const verdictLabel = ABSENSE_VERDICT[qc.verdict] || qc.verdict;
  flags.push(el('div', {}, el('span', { class: 'muted' }, 'abSENSE: '), verdictLabel,
    qc.hdf === 'detection_failure_plausible'
      ? el('div', { class: 'muted', style: 'font-size:12px;margin-top:2px' },
          'This is about the reliability of the age call, not the gene: an ancestral homolog '
          + 'could exist and be undetectable at this sensitivity. It does not mean the gene is old.')
      : null,
    el('div', { class: 'muted', style: 'font-size:12px;margin-top:2px' },
      'How well calibrated is that? Measured for the first time in September 2026: abSENSE\u2019s '
      + 'nominal 99% prediction intervals contain 87 to 88% of true bitscores overall, but only '
      + '59 to 61% at Salix purpurea, and every abSENSE flag on this pool is a prediction at '
      + 'Salix purpurea, because that is the boundary species a Populus-specific call is made '
      + 'against. The bias there is about +19 bitscore units, in the direction that makes an '
      + 'absence look more surprising than it is. So when abSENSE supports the young call, that '
      + 'support is too generous; when it retracts the young call, the retraction is conservative. '
      + 'The coverage figure is also the optimistic bound: it is measured on '
      + 'genes that do have a Salix homolog, which are by definition not the genes in this pool.')));
  // Other evidence: shown on every gene that carries it, never folded into `post`. The
  // synteny-contest routes belong to the synteny checks upstream but are not one of the two
  // QC flags (ms2-lsg s00c keeps them out of `unflagged`).
  const other = [];
  if (qc.syn_any) {
    other.push(el('div', {},
      'A 12-genome OrthoFinder rerun disputes this gene\'s orthogroup slot.'));
  }
  if (qc.short) {
    other.push(el('div', {}, `Short ORF (${qc.plen} aa), a known annotation-artifact class.`));
  }
  // `stable` is the source's sens_stable column: whether this gene's synteny call is the same
  // in the default and the --ultra-sensitive DIAMOND runs of the 12-genome GENESPACE rerun
  // (ms2-lsg analyses/repro_audit_20260917, section E3). The synteny card describes the same
  // field. Until 2026-09-21 this line said the chapter's headline result depended on the gene,
  // which the column never measured.
  if (qc.stable === false) {
    other.push(el('div', {}, el('span', { class: 'badge' }, 'sensitivity'),
      ' the synteny call differs between the default and --ultra-sensitive DIAMOND runs.'));
  }
  return el('div', {},
    el('div', {}, ...flags),
    el('div', { style: 'margin-top:6px' },
      el('span', { class: 'muted' }, 'In the post-QC set (neither QC flag): '),
      el('strong', {}, qc.post ? 'yes' : 'no')),
    other.length
      ? el('div', { style: 'margin-top:6px' },
          el('div', { class: 'muted' }, 'Other evidence, not QC flags:'), ...other)
      : null,
    el('div', { class: 'muted', style: 'font-size:11px;margin-top:4px' },
      'Flags, not filters: nothing here removes or renumbers a gene, and the pool size is '
      + 'never set by these fields; every claim about this gene should cite the flag alongside it. '
      + 'The two QC checks are abSENSE and the synteny window search.'));
}

function geneCard(g) {
  const p = g.ps;
  if (!p) return null;
  const kv = el('dl', { class: 'kv' });
  const add = (k, ...v) => { kv.append(el('dt', {}, k), el('dd', {}, ...v)); };

  add('Stratum', p.name
    ? el('span', {}, p.name, p.rank != null ? el('span', { class: 'muted' }, ` · PS${p.rank} of 19`) : null,
        p.rank === 19
          ? el('div', { class: 'muted', style: 'font-size:12px;margin-top:2px', 'data-ps19-caveat': '' },
              'PS19 means no homolog was detected outside this assembly in the species searched. '
              + 'That panel lacks P. alba, one parent of this hybrid, so a PS19 gene is most likely '
              + 'inherited from one parent species rather than born in the hybrid. It is a detection '
              + 'and sampling result, not proof of a new gene.')
          : null)
    : el('span', { class: 'muted' }, 'unplaced: genEra could not assign a stratum'));
  if (p.tr != null) {
    add('Taxonomic representativeness', `${p.tr}%`,
      el('div', { class: 'muted', style: 'font-size:12px;margin-top:2px' },
        'Share of the stratum’s lineages with a detectable homolog. Low values mean the '
        + 'assignment rests on few species.'));
  }
  if (p.rank != null) {
    add('LSG pool', p.rank >= LSG_MIN
      ? el('strong', {}, `yes, PS${p.rank}`)
      : el('span', { class: 'muted' }, `no (pool is PS ≥ ${LSG_MIN})`));
  }
  // The frozen 2026-04 genEra call, where its rank differs from the pool's call shown above.
  // Provenance only (C5, Chen 2026-09-21 evening): it drives no filter, count or default view,
  // and it stays so a reader holding an older figure can see why a gene moved. Named "the frozen
  // genEra call", never "withdrawn" (its numbers are correct statements about that call) and
  // never a "12-species panel" call (MS2 audit A38: it searched NCBI nr only).
  if (p.frz && p.frz.rank !== undefined) {
    const left = p.frz.lsg && !p.lsg;
    add('Frozen call',
      el('div', {},
        el('div', {}, 'The frozen 2026-04 genEra call, superseded as the pool, placed this gene at ',
          el('strong', {}, p.frz.rank != null ? `PS${p.frz.rank}` : 'no rank'),
          p.frz.name ? ` (${p.frz.name})` : '', '.'),
        left
          ? el('div', { style: 'margin-top:3px' },
              el('span', { class: 'badge' }, 'left the pool'),
              ' this gene was a lineage-specific candidate in the frozen call and is not one in '
              + 'the current pool. It left by a rank change, not by failing a QC check, so it '
              + 'carries no QC flags here.')
          : null,
        (p.frz.rank === 19 && p.rank === 18)
          ? el('div', { style: 'margin-top:3px' },
              'PS19 (no homolog detected outside this assembly) in the frozen call, PS18 (Populus) '
              + 'in the current one. 179 pool genes do this; none moves the other way.')
          : null,
        el('div', { class: 'muted', style: 'font-size:12px;margin-top:3px' },
          'The current call is the genEra rerun that added Idesia polycarpa, the only Salicaceae '
          + 'genome available from outside Saliceae. That rerun was also the first to search the '
          + 'full comparison panel (the frozen call searched NCBI nr only), so a rank change is '
          + 'not always due to Idesia.')));
  }
  if (g.qc) {
    kv.append(el('dt', {}, 'Evidence quality'), el('dd', {}, qcSummary(g.qc)));
  }
  // The founder family comes from the same genEra run as the rank above (the rerun that added
  // Idesia; sources/founder_families_rerun.tsv.gz, from ms2-lsg). Until 2026-09-26 it came from
  // the frozen run, and the join made a family look younger than its own member on 419 genes.
  // Not shown on 45 genes, all with no rank in the pool's run: 6 have no family assignment and
  // 39 have a family but no gene age. Neither is "no family" for all 45, so the card gives the
  // one reason that is true of every one of them.
  if (p.fnd && p.rank != null && p.fnd.rank <= p.rank) {
    const s = p.fnd.fam, o = p.fnd.rank;
    add('Gene family',
      el('span', {}, s === 1 ? el('strong', {}, 'singleton') : `${fmt(s)} members`,
        o != null && o !== p.rank
          ? el('span', { class: 'muted' }, ` · family originates at PS${o}`)
          : null),
      // A lineage-specific gene in an old family is a diverged or newly added member of that
      // family, not a new family: 49 pool genes, 16 of them in families from cellular organisms.
      (p.lsg && o != null && o < LSG_MIN)
        ? el('div', { class: 'muted', style: 'font-size:12px;margin-top:2px' },
            `A young gene in an older family: its own oldest detectable homolog is at PS${p.rank}, `
            + `but its family in this genome goes back to PS${o}. That makes it a diverged or newly added member `
            + 'of an old family, not a lineage-specific family.')
        : null);
  } else {
    add('Gene family', el('span', { class: 'muted' },
      'family origin not shown: the current pool run gives this gene no rank'));
  }

  const card = el('div', { class: 'card' },
    el('h2', {}, 'Phylostratigraphy', el('span', { class: 'tag' }, 'genEra')), kv);

  if (p.amb) {
    card.append(el('div', { class: 'warn' },
      el('strong', {}, 'genEra flagged this gene as ambiguous. '),
      'Its stratum could not be resolved to one node; the candidates are ',
      el('span', { class: 'mono', style: 'font-size:12px' }, p.amb.join(' · ')),
      '. The rank shown above is genEra’s pick among these, not a settled call.'));
  }
  if (p.rank == null) {
    card.append(el('div', { class: 'warn' },
      'genEra returned no rank for this gene. It is excluded from the LSG pool and from every '
      + 'stratum count, rather than being treated as youngest or oldest.'));
  }
  return card;
}

// One ladder per haplotype, on a SHARED log scale so the two can be read against each other.
// They used to be summed into one "both haplotypes" ladder; the lab rule is that HAP1 and HAP2
// are reported separately, and the pool is 776 vs 835, not one number.
function ladderChart(rows, hap, max) {
  // Counts span three orders of magnitude (cellular organisms down to fabids). A linear axis
  // would render most strata as invisible slivers, so the bar length is log-scaled and the
  // true count is printed beside every bar.
  const lg = (n) => (n <= 0 ? 0 : Math.log10(n) / Math.log10(max));
  const ns = rows.map((r) => r.n);
  return el('div', {},
    el('div', { class: 'barchart' }, ...rows.slice().reverse().map((r) => {
      const lsg = r.rank >= LSG_MIN;
      return el('div', { class: 'barrow', style: 'grid-template-columns:150px 1fr 60px' },
        // the bar was a dead end; it is now the way into Browse for that stratum
        el('a', { class: 'lbl', href: browseHref({ hap: hap.toLowerCase(), ps: r.rank }),
                  title: `Browse the ${fmt(r.n)} ${hap} genes in PS${r.rank}` },
          lsg ? el('b', {}, r.name) : r.name,
          el('span', { class: 'muted' }, ` PS${r.rank}`)),
        el('div', { class: 'bartrack' },
          el('div', {
            class: 'barfill',
            style: `width:${Math.max(lg(r.n) * 100, 1)}%;background:${lsg ? 'var(--c-priv)' : 'var(--accent)'}`,
          })),
        el('span', { class: 'val' }, fmt(r.n)));
    })),
    el('p', { class: 'muted', style: 'font-size:11.5px;margin:9px 0 0' },
      `Bar length is log-scaled, on the same scale for both haplotypes: ${hap} counts run from `
      + `${fmt(Math.max(...ns))} to ${fmt(Math.min(...ns))}, and a linear axis would hide most `
      + 'strata. The number beside each bar is the true count. Amber marks the LSG pool '
      + '(PS ≥ 18). Click a stratum to browse its genes.'));
}

async function overview(man) {
  const ps = man.phylostrat;
  if (!ps) return el('div', { class: 'empty' }, 'Phylostratigraphy data not built.');
  const h1 = ps.ladder.hap1, h2 = ps.ladder.hap2;
  const sum = (rows, f) => rows.filter(f).reduce((a, r) => a + r.n, 0);
  const at = (rows, k) => sum(rows, (r) => r.rank === k);
  const placed = { HAP1: sum(h1, () => true), HAP2: sum(h2, () => true) };
  const pool = { HAP1: sum(h1, (r) => r.rank >= LSG_MIN), HAP2: sum(h2, (r) => r.rank >= LSG_MIN) };
  const max = Math.max(...h1.map((r) => r.n), ...h2.map((r) => r.n));
  const naH1 = (ps.by_rank.hap1 || {}).NA || 0;
  const naH2 = (ps.by_rank.hap2 || {}).NA || 0;
  // The pool is the genEra rerun that added Idesia (C5, Chen 2026-09-21 evening, reconfirmed
  // 2026-09-26). The frozen call is named as superseded, with its count, so a reader arriving
  // from an older figure finds it; it never drives a number on this page's charts.
  const frz = ((man.lsg_pool || {}).frozen || {});
  const frzN = frz.at_rank_18_or_above || {};
  const left = frz.left_the_pool || {};
  const post = ((man.lsg_qc || {}).per_haplotype || {}).in_the_chapter_post_qc_set || {};
  const poolStat = (hap, rows) => el('div', { class: 'stat' },
    el('div', { class: 'n' }, fmt(pool[hap])),
    el('div', { class: 'l' }, `${hap} LSG pool: ${fmt(at(rows, 18))} PS18 · ${fmt(at(rows, 19))} PS19`));

  return el('div', {},
    el('h1', {}, 'Phylostratigraphy'),
    el('p', { class: 'sub' },
      'genEra assigns each gene to the oldest phylostratum in which a homolog is detectable, '
      + 'over 19 strata from cellular organisms to this hybrid. The '
      + 'lineage-specific candidate pool is PS ≥ 18 in the genEra rerun that added Idesia polycarpa: '
      + `${fmt(pool.HAP1)} HAP1 and ${fmt(pool.HAP2)} HAP2 genes. PS17 Saliceae is flanking `
      + 'context, not pool. Two QC checks flag genes without removing them, abSENSE and the '
      + `synteny window search; ${fmt(post.HAP1)} HAP1 and ${fmt(post.HAP2)} HAP2 carry neither `
      + 'flag, the post-QC set. The frozen 2026-04 genEra call, superseded as the '
      + `pool, had ${fmt(frzN.hap1)} HAP1 and ${fmt(frzN.hap2)} HAP2; the ${fmt(left.hap1)} HAP1 `
      + `and ${fmt(left.hap2)} HAP2 genes that left say so on their own page.`),
    el('div', { class: 'stat-row' },
      poolStat('HAP1', h1),
      poolStat('HAP2', h2),
      el('div', { class: 'stat' }, el('div', { class: 'n' }, fmt(placed.HAP1 + placed.HAP2)),
        el('div', { class: 'l' }, `genes placed: ${fmt(placed.HAP1)} HAP1 · ${fmt(placed.HAP2)} HAP2`)),
      el('div', { class: 'stat' }, el('div', { class: 'n' }, fmt(ps.n_ambiguous)), el('div', { class: 'l' }, 'ambiguous assignments')),
      el('div', { class: 'stat' }, el('div', { class: 'n' }, fmt(naH1 + naH2)),
        el('div', { class: 'l' }, `unplaced by genEra: ${fmt(naH1)} HAP1 · ${fmt(naH2)} HAP2`)),
    ),
    el('div', { class: 'hap-pair' },
      el('div', { class: 'card' },
        el('h2', {}, 'Genes per stratum', el('span', { class: 'tag' }, 'HAP1')),
        ladderChart(h1, 'HAP1', max)),
      el('div', { class: 'card' },
        el('h2', {}, 'Genes per stratum', el('span', { class: 'tag' }, 'HAP2')),
        ladderChart(h2, 'HAP2', max))),
    el('div', { class: 'card' },
      el('h2', {}, 'Reading this carefully'),
      el('ul', { style: 'margin:0;padding-left:18px;font-size:13.5px;color:var(--ink-2);line-height:1.6' },
        el('li', {}, el('strong', {}, `${fmt(ps.n_ambiguous)} genes carry an ambiguity flag. `),
          'genEra could not resolve them to one node and lists several candidate strata. They '
          + 'still get a rank, so the flag is the only thing separating a confident call from a '
          + 'coin flip, every affected gene page says so.'),
        el('li', {}, el('strong', {}, `${fmt(naH1)} HAP1 and ${fmt(naH2)} HAP2 genes are unplaced. `),
          'No rank at all in the current call. They are excluded from the pool and from these '
          + 'counts rather than defaulted in either direction.'),
        el('li', {}, el('strong', {}, 'A young stratum is not proof of a young gene. '),
          'Phylostratigraphy under-detects fast-evolving homologs, so PS ≥ 18 is a candidate '
          + 'pool, not a verdict. That is why the synteny module\u2019s independent orthology '
          + 'call sits beside it on every gene page.'),
        el('li', {}, el('strong', {}, 'No homolog found can also mean no genome searched. '),
          'A homolog is only detectable in species that were searched. The panel lacks P. alba, '
          + 'one parent of this hybrid, so most PS19 genes (no homolog outside this assembly) are '
          + 'likely inherited from one parent species, not born in the hybrid. Adding Idesia '
          + 'removed 66 HAP1 and 48 HAP2 genes from the pool for the same reason.'),
        el('li', {}, el('strong', {}, 'PS18 and PS19 are reported separately. '),
          'Short or poorly supported gene models are a known source of false young genes, and '
          + 'the youngest stratum is where they would collect, so the split is shown rather than '
          + 'pooled.'))),
  );
}

export default {
  id: 'phylostrat',
  label: 'Phylostrat',
  geneCard,
  overview,
  filters: [
    { id: 'ps_lsg', label: 'LSG pool (PS ≥ 18)', options: ['yes', 'no'],
      test: (g, v) => (g.ps && g.ps.rank >= LSG_MIN ? 'yes' : 'no') === v },
    { id: 'ps_rank', label: 'Phylostratum ≥',
      test: (g, v) => g.ps && g.ps.rank != null && g.ps.rank >= parseInt(v, 10) },
  ],
};
