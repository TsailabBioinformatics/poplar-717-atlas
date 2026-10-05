#!/usr/bin/env python3
"""TODO 29: per-gene Ka/Ks against Populus trichocarpa.

New to the atlas. The only selection-adjacent number here before was HAP1<->HAP2 dS on the
allele card, which compares two haplotypes of one clone; this compares 717 to another
species, which is the quantity a visitor asking "is this gene under purifying selection"
expects.

NG86 IS THE PER-GENE VALUE, and that was a measurement rather than a default. On these
alignments NG86 and YN00 omega correlate at r = 0.07 to 0.36 per gene while their class
medians agree to within 0.06; median dS is about 0.04, close enough that YN00's estimate is
unstable gene by gene. NG86 is the one with external support -- r = 0.937 against the
atlas's existing, independently derived allele omega, where YN00 gives 0.083 on the same
series. `yn00_omega` still ships, but only so the card can SHOW that spread instead of
asserting a precision the data does not have. See the provenance sidecar.

ABSENCE IS INFORMATIVE HERE, so it is stored rather than left blank: a gene with no
reciprocal-best trichocarpa ortholog gets `ks.no = 1`. 1 of 681 PS19 genes has any pair at
all, because a hybrid-only gene has nothing to measure against.

The source is gated at generation against the owning analysis's published per-class table.
"""
import csv
import gzip
import io
import json
import os
import statistics as st
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
SRC = os.path.join(ROOT, "sources", "kaks_ptri.tsv.gz")

# medians the source builder printed and the sidecar records; re-asserted from the shards
EXPECTED_MEDIAN = {"PS1_2": 0.2351, "PS3_12": 0.3272, "PS13_15": 0.4748,
                   "PS16_17": 0.5662, "LSG_PS18_populus": 0.7204}


def num(v):
    try:
        f = float(v)
        return f if f == f else None
    except (TypeError, ValueError):
        return None


def gene_class(rk):
    if rk is None:
        return None
    if rk <= 2:  return "PS1_2"
    if rk <= 12: return "PS3_12"
    if rk <= 15: return "PS13_15"
    if rk <= 17: return "PS16_17"
    if rk == 18: return "LSG_PS18_populus"
    return "LSG_PS19_hybrid"


def shard_paths(root):
    for dp, _, fs in os.walk(root):
        for fn in sorted(fs):
            if fn.endswith(".json"):
                yield os.path.join(dp, fn)


def main():
    with gzip.open(SRC, "rt") as fh:
        rows = {r["gene_id"]: r for r in csv.DictReader(io.StringIO(fh.read()), delimiter="\t")}
    print(f"source: {len(rows)} genes with a trichocarpa pair")

    patched, no_pair, agree, by_cls = 0, 0, Counter(), defaultdict(list)
    for hap in ("hap1", "hap2"):
        for p in shard_paths(os.path.join(DATA, "genes", hap)):
            shard = json.load(open(p))
            for gid, rec in shard.items():
                r = rows.get(gid)
                if r is None:
                    rec["ks"] = {"no": 1}
                    no_pair += 1
                    continue
                w = num(r["ng86_omega"])
                ks = {"id": r["ptri_id"], "n": int(float(r["n_codons"])),
                      "pid": round(float(r["pct_identity"]), 2)}
                for key, col in (("ka", "ng86_dn"), ("kss", "ng86_ds"), ("w", "ng86_omega")):
                    v = num(r[col])
                    if v is not None:
                        ks[key] = round(v, 4)
                v2 = num(r["yn00_omega"])
                if v2 is not None:
                    ks["w2"] = round(v2, 4)
                # does the RBH ortholog agree with the atlas's own syntenic 1:1 call?
                pt = rec.get("ptri") or {}
                if pt.get("one2one") == 1:
                    same = pt.get("id") == r["ptri_id"]
                    ks["ag"] = 1 if same else 0
                    agree["same" if same else "different"] += 1
                rec["ks"] = ks
                patched += 1
                c = gene_class((rec.get("ps") or {}).get("rank"))
                if c and w is not None:
                    by_cls[c].append(w)
            with open(p, "w") as fh:
                json.dump(shard, fh, separators=(",", ":"), sort_keys=True)

    print(f"patched {patched} genes; {no_pair} marked as having no trichocarpa ortholog")
    print(f"  RBH ortholog vs the atlas's syntenic 1:1 call: {dict(agree)} "
          f"({100 * agree['same'] / max(sum(agree.values()), 1):.1f}% agree)")
    assert patched == len(rows), f"{len(rows)} source rows but {patched} patched"

    bad = []
    for c, want in EXPECTED_MEDIAN.items():
        got = st.median(by_cls[c]) if by_cls[c] else None
        print(f"  {c:18s} n={len(by_cls[c]):6d} median omega "
              f"{got:.4f}" if got is not None else f"  {c:18s} n=0")
        if got is None or abs(got - want) > 5e-4:
            bad.append((c, want, got))
    assert not bad, (
        f"the per-gene medians read back off the shards do not reproduce the source's own "
        f"per-class table: {bad}. Either the join dropped genes or the class map moved.")
    print("  medians reproduce the source's per-class table")

    man_path = os.path.join(DATA, "meta", "manifest.json")
    man = json.load(open(man_path))
    man["kaks"] = {
        "target": "Populus trichocarpa, reciprocal-best-hit ortholog",
        "n_genes": patched,
        "n_without_an_ortholog": no_pair,
        "estimator": "NG86 (Nei-Gojobori 1986) for the per-gene value",
        "why_not_yn00": "NG86 and YN00 omega correlate at r = 0.07 to 0.36 PER GENE on these "
                        "alignments while their class medians agree to within 0.06 -- median "
                        "dS is about 0.04, too close for YN00 to be stable gene by gene. NG86 "
                        "is the one with external support: r = 0.937 against the atlas's "
                        "independently derived allele omega, against YN00's 0.083. yn00_omega "
                        "ships only to show the spread, never as a per-gene answer.",
        "median_omega_by_class": EXPECTED_MEDIAN,
        "reading": "Monotone with youth and well below the neutral 1. Young genes are "
                   "constrained, just less constrained -- not neutrally evolving, which is "
                   "what a spurious open reading frame would look like.",
        "selection_bias": "Pairing rate falls with youth (87.6% ancient, 39.1% PS18, ~0% "
                          "PS19), so the young-gene omega is measured on the conserved end of "
                          "its stratum and if anything understates the relaxation.",
        "ps19": "1 of 681 PS19 genes has any trichocarpa pair. A structural absence, not "
                "missing data -- a hybrid-only gene has no ortholog to measure against.",
        "ortholog_agreement_with_syntenic_call": dict(agree),
    }
    with open(man_path, "w") as fh:
        json.dump(man, fh, indent=1, sort_keys=True)
    print("manifest: kaks block written")


if __name__ == "__main__":
    main()
