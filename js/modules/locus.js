// Module: locus. Now only the RNA-seq gene-model evidence card; the locus drawing itself is
// the shared two-haplotype browser in the "both haplotypes" card (mirror.js).
import { el } from '../core/dom.js';

// ---- RNA-seq gene-model evidence -----------------------------------------------------------
// Two readouts from 94 public RNA-seq libraries, both FLAGS and never filters.
//
// WHY THE READ COUNT IS NOT OPTIONAL. Intron detection tracks DEPTH: under 20 reads the
// per-intron detection rate is below 27%, over 5,000 reads it is 98.8%. Lineage-specific genes
// are shallow, so "0 of 4 introns seen" for an LSG with 0 reads is a statement about its
// expression and nothing else. Shown without the depth it would read as evidence against the
// gene model, which the data does not support.
//
// NO THRESHOLD IS APPLIED to the gap-coverage discriminator. An intron is spliced OUT, so a
// genuine split shows a junction with LOW coverage inside the gap, while high gap coverage is
// untidy readthrough with a defensible annotation. The owning analysis never published a
// cutoff, so both counts are shown and the reader judges.
function geneModelCard(g) {
  const gm = g.gm;
  if (!gm) return null;
  const kv = el('dl', { class: 'kv' });
  const add = (k, ...v) => { kv.append(el('dt', {}, k), el('dd', {}, ...v)); };
  const N = (v) => (v == null ? '–' : v.toLocaleString());

  const i = gm.in;
  if (i) {
    const shallow = i.rd < 20;
    add('Introns observed',
      el('div', {},
        el('div', {},
          el('strong', {}, `${i.seen} of ${i.n}`),
          ' annotated intron', i.n === 1 ? '' : 's',
          ' seen as spliced reads'),
        el('div', { class: 'muted', style: 'font-size:12px;margin-top:3px' },
          `${N(i.rd)} read${i.rd === 1 ? '' : 's'} on this gene across 94 libraries.`),
        shallow
          ? el('div', { class: 'warn', style: 'margin-top:6px' },
              el('strong', {}, 'Too shallow to conclude anything. '),
              'Below 20 reads, fewer than 27% of annotated introns are detected even for genes '
              + 'whose models are otherwise well supported (20 reads is a display breakpoint, not '
              + 'a calibrated threshold). ',
              i.seen < i.n
                ? 'The unseen introns here are expected at this depth and are not evidence '
                  + 'against the gene model.'
                : '')
          : el('div', { class: 'muted', style: 'font-size:12px;margin-top:3px' },
              i.seen === i.n
                ? 'Every annotated intron is supported by spliced reads.'
                : 'Detection rises steeply with depth, 98.8% of introns are seen above 5,000 '
                  + 'reads, so read the missing ones against this gene\u2019s coverage.')));
  }

  if (gm.sp) {
    const sp = gm.sp;
    const lsgArm = sp.arm === 'A_lsg_conserved';
    add('Possible split model',
      el('div', {},
        el('div', {},
          el('span', { class: 'badge' }, 'flag'),
          ' spliced reads join an annotated exon boundary of this gene to one of ',
          el('a', { href: `#/gene/${sp.p}`, class: 'mono' }, sp.p), '.'),
        el('table', { class: 'data', style: 'margin-top:7px' },
          el('tbody', {},
            el('tr', {}, el('td', {}, 'Junctions at exon boundaries'),
              el('td', { class: 'num mono' }, N(sp.je))),
            el('tr', {}, el('td', {}, 'All spliced reads joining the pair'),
              el('td', { class: 'num mono' }, N(sp.jr))),
            el('tr', {}, el('td', {}, 'Reads inside the intergenic gap'),
              el('td', { class: 'num mono' }, N(sp.gr))),
            el('tr', {}, el('td', {}, 'Libraries with the junction'),
              el('td', { class: 'num mono' }, `${N(sp.lib)} of 94`)),
            el('tr', {}, el('td', {}, 'Gap between the two genes'),
              el('td', { class: 'num mono' }, `${N(sp.gap)} bp`)))),
        el('div', { class: 'muted', style: 'font-size:12px;margin-top:6px' },
          el('strong', {}, 'How to read the gap count. '),
          'An intron is spliced out, so if these really are one gene the gap should be nearly '
          + 'empty of reads. Heavy coverage inside the gap points instead at readthrough '
          + 'transcription or an unannotated UTR, where the annotation is defensible and only '
          + 'the transcription is untidy. No threshold is applied here; both numbers are shown '
          + 'because the published analysis did not set one.'),
        el('div', { class: 'muted', style: 'font-size:12px;margin-top:4px' },
          el('strong', {}, 'A flag, not a verdict, and not a lineage-specific problem. '),
          '1,422 of the 1,435 flagged pairs are two conserved genes; only 13 involve a '
          + 'lineage-specific gene. At matched gap and depth, lineage-specific pairs are flagged '
          + 'LESS often than conserved ones (0.96% compared with 3.61%). ',
          lsgArm ? 'This pair is one of the 13.' : 'This pair is a conserved-with-conserved one.')));
  } else if (gm.nev === 0 && gm.np > 0) {
    add('Possible split model',
      el('div', { class: 'muted' },
        el('strong', {}, 'Not testable for this gene. '),
        `Neither of its ${gm.np === 1 ? 'adjacent pair' : `${gm.np} adjacent pairs`} had enough `
        + 'reads to look for a joining junction, so the absence of a flag here is not evidence '
        + 'that the two genes are separate. 9,079 genes are in this position.'));
  }

  if (!kv.children.length) return null;
  // The last sentence used to read "An 827-library expansion is in progress and will supersede
  // this panel". The enlarged panel is 1,755 libraries on ms2-lsg PR #20, unmerged and held for
  // Chen's review (2026-09-21), so neither the size nor "will supersede" was true.
  return el('div', { class: 'card', 'data-card': 'gene-model' },
    el('h2', {}, 'Gene model, from RNA-seq',
      el('span', { class: 'tag' }, '94 public libraries')), kv,
    el('p', { class: 'muted', style: 'font-size:11px;margin:10px 0 0' },
      'Spliced-read evidence only, the strong readout. Read-pair spanning and gap '
      + 'occupancy alone are not used to call anything, because readthrough and unspliced '
      + 'pre-mRNA both produce them. Each haplotype was aligned to its own genome, so HAP1 and '
      + 'HAP2 are not independent replicates.'));
}

// The single-haplotype canvas strip that lived here is gone: the "both haplotypes" card
// (mirror.js) now uses the shared locus browser, which shows the same exon structure when
// zoomed in and pans. This module keeps the RNA-seq gene-model evidence card.
function geneCard(g) {
  return geneModelCard(g);
}

export default { id: 'locus', label: 'Locus', geneCard };
