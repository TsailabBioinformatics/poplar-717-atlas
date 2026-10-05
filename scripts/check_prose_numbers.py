#!/usr/bin/env python3
"""Recompute the numbers the About page states in prose, and diff them against the page.

WHY. `js/modules/about.js` is the page a reader is told to consult before citing anything, and
several of its figures are typed into sentences rather than read from the manifest: the
per-haplotype gene counts, the HAP1/HAP2 Pfam gap, the count of syntelog pairs whose omega is
undefined. `check_aggregates.py` recomputes every manifest aggregate and diffs it, but prose is
outside its reach -- so exactly the page that exists to keep the atlas honest is the one place a
number can rot unobserved. This repo has already shipped a withdrawn figure on a rendered page
for two days (0.3.11, the 1,725 LSG pool), which is the same failure with a different surface.

WHAT IT DOES. Recomputes each claim from the committed shards and asserts the formatted string
appears in about.js. It does NOT rewrite the prose: a sentence is a claim with context around
it, and a script that silently edited the number could make the surrounding words false. When
this fails, a human reads the sentence and decides what it should now say.
"""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
ABOUT = os.path.join(ROOT, "js", "modules", "about.js")

fails = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}{'  — ' + detail if detail else ''}")
    if not ok:
        fails.append(name)


PROTEIN = os.path.join(ROOT, "js", "modules", "protein.js")
DEEPEST = "PXD025636"


def protein_claims():
    """The Protein page's age-gradient sentence, per haplotype, recomputed from the shards.

    Added 2026-09-21. That sentence said "0.12% of the lineage-specific pool", pooled over both
    haplotypes: the published Tier 2 figure (2 of 1,611), one of whose hits the owning analysis
    withdrew on 2026-09-09. The shards never carried the withdrawn hit -- only the sentence did,
    and nothing compared the two. Same failure as the About page's numbers, on another page.
    """
    classes = [("ancient", lambda r: r is not None and r <= 2),
               ("conserved", lambda r: r is not None and 3 <= r <= 12),
               ("rosid", lambda r: r is not None and 13 <= r <= 15),
               ("Salicaceae", lambda r: r is not None and 16 <= r <= 17)]
    n, det = {}, {}
    pool = {h: [0, 0] for h in ("hap1", "hap2")}            # [detected in DEEPEST, genes]
    union = {(h, k): [0, 0] for h in ("hap1", "hap2") for k in (18, 19)}  # any dataset
    for hap in ("hap1", "hap2"):
        for dp, _, fs in os.walk(os.path.join(DATA, "genes", hap)):
            for fn in fs:
                if not fn.endswith(".json"):
                    continue
                for rec in json.load(open(os.path.join(dp, fn))).values():
                    ps = rec.get("ps") or {}
                    # the owning analysis ranks by the Idesia rerun, which is ps.rank again since
                    # C5; a leftover ps.rr would mean the C4 arrangement and is refused
                    assert "rr" not in ps, "ps.rr is the C4 arrangement; run patch_lsg_pool.py"
                    r = ps.get("rank")
                    ms = (rec.get("prot") or {}).get("ms") or {}
                    hit = DEEPEST in (ms.get("ds") or [])
                    for name, f in classes:
                        if f(r):
                            n[(hap, name)] = n.get((hap, name), 0) + 1
                            det[(hap, name)] = det.get((hap, name), 0) + hit
                    if r is not None and r >= 18:
                        pool[hap][0] += hit
                        pool[hap][1] += 1
                        union[(hap, r)][0] += (ms.get("n") or 0) > 0
                        union[(hap, r)][1] += 1
    src = open(PROTEIN).read()
    pc = lambda h, c: f"{100 * det[(h, c)] / n[(h, c)]:.1f}%"
    wants = [f"{pc('hap1', 'ancient')} of HAP1 and {pc('hap2', 'ancient')} of HAP2 ancient"]
    for c in ("conserved", "rosid", "Salicaceae"):
        wants.append(f"{pc('hap1', c)[:-1]}% and {pc('hap2', c)[:-1]}% of {c}")
    wants.append(f"{pool['hap1'][0]} of {pool['hap1'][1]:,} HAP1 and "
                 f"{pool['hap2'][0]} of {pool['hap2'][1]:,} HAP2 lineage-specific pool genes")
    wants.append(f"{union[('hap1', 18)][0]} of {union[('hap1', 18)][1]} HAP1 and "
                 f"{union[('hap2', 18)][0]} of {union[('hap2', 18)][1]} HAP2 PS18")
    wants.append(f"{union[('hap1', 19)][0]} of {union[('hap1', 19)][1]} HAP1 and "
                 f"{union[('hap2', 19)][0]} of {union[('hap2', 19)][1]} HAP2 PS19")
    # Join the JS string concatenation ('...' + '...') so a claim can span a source line break.
    flat = re.sub(r"'\s*\n\s*\+\s*'", "", src)
    for w in wants:
        ok = w in flat
        check(f"protein.js states: {w}", ok, "" if ok else "not found -- recomputed from the shards")


def main():
    n_hap = {"hap1": 0, "hap2": 0}
    pfam = {"hap1": 0, "hap2": 0}
    ks0 = 0
    for hap in ("hap1", "hap2"):
        for dp, _, fs in os.walk(os.path.join(DATA, "genes", hap)):
            for fn in sorted(fs):
                if not fn.endswith(".json"):
                    continue
                for gid, rec in json.load(open(os.path.join(dp, fn))).items():
                    n_hap[hap] += 1
                    if ((rec.get("ann") or {}).get("pf") or []):
                        pfam[hap] += 1
                    if (rec.get("allele") or {}).get("ks0"):
                        ks0 += 1

    src = open(ABOUT).read()
    fmt = lambda n: f"{n:,}"

    claims = [
        ("HAP1 gene count", fmt(n_hap["hap1"])),
        ("HAP2 gene count", fmt(n_hap["hap2"])),
        ("HAP1 genes with a Pfam domain", fmt(pfam["hap1"])),
        ("HAP2 genes with a Pfam domain", fmt(pfam["hap2"])),
    ]
    print(f"\nrecomputed from the shards: HAP1 {n_hap['hap1']:,} / HAP2 {n_hap['hap2']:,}, "
          f"Pfam {pfam['hap1']:,} / {pfam['hap2']:,}, allele pairs with Ks 0: {ks0:,}\n")

    for label, s in claims:
        # the page writes these with fmt(), so they appear with thousands separators, but a
        # raw literal is also accepted -- what matters is that the true number is on the page
        bare = s.replace(",", "")
        present = (s in src) or (bare in src)
        check(f"about.js states the right {label}", present,
              f"expected {s}" if not present else s)

    # The page counts PAIRS; the flag lives on both members, so the gene count must be exactly
    # twice the pair count. Assert that first: an odd total would mean one side of some pair
    # carries the flag and the other does not, which is a data defect the page's number would
    # quietly average over. (This is how the first draft of this check got it wrong -- it
    # counted genes, reported 1,432 against the page's 716, and the page was right.)
    check("the Ks-0 allele flag is symmetric across each pair", ks0 % 2 == 0,
          f"{ks0:,} genes carry it, which is odd -- some pair is flagged on one side only"
          if ks0 % 2 else f"{ks0:,} genes = {ks0 // 2:,} pairs")
    pairs = ks0 // 2
    present = (fmt(pairs) in src) or (str(pairs) in src)
    check("about.js states the right count of syntelog PAIRS with Ks 0", present,
          f"recomputed {pairs:,} pairs" + ("" if present else " -- not found in the page"))

    protein_claims()

    print()
    if fails:
        print(f"{len(fails)} prose number(s) disagree with the data.")
        print("Read the sentence and decide what it should say -- do not just swap the digits,")
        print("because the words around a number can stop being true when the number moves.")
        return 1
    print("every number the About page states in prose matches the committed data")
    return 0


if __name__ == "__main__":
    sys.exit(main())
