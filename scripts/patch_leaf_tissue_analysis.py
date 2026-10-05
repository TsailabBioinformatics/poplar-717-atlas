#!/usr/bin/env python3
"""One-time correction: re-group the WT-control tissue baseline (and tau) by
`tissue_analysis` instead of `tissue`, which silently pooled leaf_young + leaf_old
(46 samples, all study 1) into one 'leaf' bucket. See docs/DATA_AUDIT.md.

This is NOT a replacement for `build_expression.py` -- that script is fixed at the
source in the same commit, for the next full Sapelo2 rebuild. This script exists
because the raw sources (BED files, the 219MB TPM matrix on Sapelo2) are not needed
to fix this: the already-shipped per-sample TPM vectors in data/expr/buckets/ are
untouched by this bug (it is purely an AGGREGATION bug), so re-aggregating them
locally, from data already committed to this repo, reproduces exactly what a full
rebuild would produce for this one layer. Verified by full re-derivation of one gene
against the manual calculation before this was trusted to run genome-wide.

Touches: data/expr/samples.json (adds tissue_analysis), data/expr/studies.json
(tissues list), every gene shard's "x" field (t/n/tau/top; "b" and "max" untouched,
since neither depends on tissue grouping), data/meta/manifest.json.
"""
import glob, json, os, sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")

TAU_TISSUES_EXCLUDE = {"callus"}   # identical to build_expression.py


def tau(values):
    mx = max(values)
    if mx <= 0:
        return None
    n = len(values)
    if n <= 1:
        return None
    return sum(1 - (v / mx) for v in values) / (n - 1)


def main():
    # -- 1. join tissue_analysis onto samples.json, in place ----------------------------
    sp = os.path.join(DATA, "expr", "samples.json")
    samples = json.load(open(sp))
    ta_map = {}
    with open(os.path.join(ROOT, "sources", "sample_tissue_analysis.tsv")) as fh:
        next(fh)
        for line in fh:
            label, ta = line.rstrip("\n").split("\t")
            ta_map[label] = ta
    missing = [s["label"] for s in samples if s["label"] not in ta_map]
    assert not missing, f"{len(missing)} atlas samples have no tissue_analysis mapping: {missing[:5]}"
    for s in samples:
        s["tissue_analysis"] = ta_map[s["label"]]
    # every original 'tissue' value must map to a CONSISTENT set of tissue_analysis values
    # across all its samples counted correctly -- not a hard requirement, just print it,
    # since a coarse tissue legitimately maps to >1 fine tissue (that is the whole point).
    by_tissue = defaultdict(set)
    for s in samples:
        by_tissue[s["tissue"]].add(s["tissue_analysis"])
    refined = {t: sorted(v) for t, v in by_tissue.items() if len(v) > 1}
    print("tissue values that refine into >1 tissue_analysis value:", refined)
    with open(sp, "w") as fh:
        json.dump(samples, fh, separators=(",", ":"))
    print(f"patched samples.json: {len(samples)} rows, added tissue_analysis")

    # -- 2. studies.json: recompute the per-study 'tissues' list on tissue_analysis -----
    stp = os.path.join(DATA, "expr", "studies.json")
    studies = json.load(open(stp))
    by_study = defaultdict(set)
    for s in samples:
        by_study[str(s["study"])].add(s["tissue_analysis"])
    changed_studies = 0
    for sid, st in studies.items():
        new_tissues = sorted(by_study[sid])
        if new_tissues != st["tissues"]:
            changed_studies += 1
        st["tissues"] = new_tissues
    with open(stp, "w") as fh:
        json.dump(studies, fh, separators=(",", ":"), sort_keys=True)
    print(f"patched studies.json: {changed_studies} of {len(studies)} studies' tissue list changed")

    # -- 3. WT-control index, on tissue_analysis -----------------------------------------
    wt_ctrl_idx = defaultdict(list)
    for i, s in enumerate(samples):
        if s["genotype"] == "WT" and s["ctrl"] == "control":
            wt_ctrl_idx[s["tissue_analysis"]].append(i)
    print("WT-control samples per tissue_analysis:",
          {t: len(v) for t, v in sorted(wt_ctrl_idx.items())})

    # -- 4. recompute every gene's x.t/x.n/x.tau/x.top from the SHIPPED bucket vectors --
    bucket_cache = {}

    def load_bucket(bkey):
        if bkey not in bucket_cache:
            hap, bnum = bkey.split("/")
            bucket_cache[bkey] = json.load(open(os.path.join(DATA, "expr", "buckets", hap, f"{bnum}.json")))
        return bucket_cache[bkey]

    patched = 0
    for hap in ("hap1", "hap2"):
        for dp, _, fs in os.walk(os.path.join(DATA, "genes", hap)):
            for fn in fs:
                if not fn.endswith(".json"):
                    continue
                p = os.path.join(dp, fn)
                shard = json.load(open(p))
                changed = False
                for gid, rec in shard.items():
                    x = rec.get("x")
                    if not x:
                        continue
                    vals = load_bucket(x["b"])[gid]
                    tmed, tn = {}, {}
                    for t, idxs in wt_ctrl_idx.items():
                        tv = sorted(vals[k] for k in idxs)
                        tmed[t] = round(tv[len(tv) // 2], 2)
                        tn[t] = len(idxs)
                    tau_vals = [tmed[t] for t in tmed if t not in TAU_TISSUES_EXCLUDE]
                    tv = tau(tau_vals) if tau_vals else None
                    top = max(tmed, key=tmed.get) if tmed else None
                    new_x = {"b": x["b"], "t": tmed, "n": tn,
                             "tau": round(tv, 3) if tv is not None else None,
                             "top": top, "max": x["max"]}
                    if new_x != x:
                        rec["x"] = new_x
                        changed = True
                        patched += 1
                if changed:
                    with open(p, "w") as fh:
                        json.dump(shard, fh, separators=(",", ":"), sort_keys=True)
    print(f"patched {patched} gene shard 'x' records (tissue_analysis-grouped baseline + tau)")

    # -- 5. manifest.json: tissue_wt_control_n, data_version, a correction note ---------
    mp = os.path.join(DATA, "meta", "manifest.json")
    man = json.load(open(mp))
    old_tw = man["expression"]["tissue_wt_control_n"]
    man["expression"]["tissue_wt_control_n"] = {t: len(v) for t, v in sorted(wt_ctrl_idx.items())}
    man.setdefault("corrections", []).append({
        "date": "2026-09-05",
        "what": "expression tissue baseline + tau regrouped by tissue_analysis, not tissue",
        "why": "tissue pooled leaf_young + leaf_old (46 samples, all study 1) into one 'leaf' "
               "bucket; the WT-control subset of that pool was 11 + 11 of 77",
        "before_tissue_wt_control_n": old_tw,
        "after_tissue_wt_control_n": man["expression"]["tissue_wt_control_n"],
        "detail": "docs/DATA_AUDIT.md, CHANGELOG.md",
    })
    prev = man["data_version"]
    parts = [int(p) for p in prev.split(".")]
    parts[-1] += 1
    man["data_version"] = ".".join(map(str, parts))
    with open(mp, "w") as fh:
        json.dump(man, fh, indent=1, sort_keys=True)
    print(f"manifest: data_version {prev} -> {man['data_version']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
