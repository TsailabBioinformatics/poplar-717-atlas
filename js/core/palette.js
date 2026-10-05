// Shared canvas palette. Extracted from locus.js so every canvas-drawn view (the locus
// card, the genome map) reads the SAME resolved colours and can never drift apart.
//
// Canvas fillStyle does NOT resolve CSS custom properties -- assigning 'var(--c-tandem)'
// is silently invalid and leaves the previous colour (a real bug, caught once already by
// screenshot, not by any assertion). So the tokens are read off the document and cached,
// re-read only when the theme flips.
const CLASS_VAR = {
  T: '--c-tandem', A: '--c-priv', S: '--c-core', D: '--c-disp', C: '--c-nonsyn',
};
let PALETTE = null;
export function palette() {
  if (PALETTE) return PALETTE;
  const cs = getComputedStyle(document.documentElement);
  const g = (v, fb) => (cs.getPropertyValue(v).trim() || fb);
  PALETTE = {
    T: g('--c-tandem', '#e87ba4'), A: g('--c-priv', '#eda100'),
    S: g('--c-core', '#2a78d6'), D: g('--c-disp', '#eb6834'),
    C: g('--c-nonsyn', '#1baf7a'),
    ink: g('--ink', '#111'), ink3: g('--ink-3', '#888'),
    line: g('--line-2', '#ccc'), accent: g('--accent', '#2a78d6'),
    panel: g('--panel', '#fff'),
  };
  return PALETTE;
}
matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => { PALETTE = null; });

/** A gene's letters are concatenated (e.g. "ADS"); colour by the first structural one. */
export const fillFor = (c) => {
  const p = palette();
  if (!c) return p.ink3;
  for (const k of ['T', 'A', 'S', 'D', 'C']) if (c.includes(k)) return p[k];
  return p.ink3;
};

export { CLASS_VAR };

// pLDDT bands for the 3D model. The AlphaFold DB convention -- the bands and the colours --
// so a reader who has seen one AlphaFold model reads this one without learning a legend.
// Fixed hex rather than theme tokens because the convention IS the colour. Lower bounds are
// inclusive, which makes "low" + "very low" exactly the card's "residues below 70".
export const PLDDT_BANDS = [
  { min: 90, label: 'very high', range: '≥ 90', color: '#0053D6' },
  { min: 70, label: 'confident', range: '70 to 90', color: '#65CBF3' },
  { min: 50, label: 'low', range: '50 to 70', color: '#FFDB13' },
  { min: -Infinity, label: 'very low', range: '< 50', color: '#FF7D45' },
];
export const plddtBand = (v) => PLDDT_BANDS.findIndex((b) => v >= b.min);
