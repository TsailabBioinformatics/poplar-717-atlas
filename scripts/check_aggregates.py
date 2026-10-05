#!/usr/bin/env python3
"""TODO 25's helper: recompute every shard-derivable manifest aggregate and diff it.

WHY THIS EXISTS, and why it is separate from validate_data.py. An aggregate is stale-by-default
the moment a patcher touches the fields it summarises, because the stage that WROTE it may not
be runnable here -- build_data.py wipes data/genes, and build_locus.py needs the OneDrive tree.
A stale aggregate renders fine, validates fine and is internally consistent, so nothing catches
it. That has now happened three times:

  1. the facet index, after the leaf-tissue patch (54,705 wrong tau)   -- v0.3.11
  2. data/locus/overview.json, after the pool correction (104 bins)    -- v0.6.0
  3. manifest.summary.ps_rank, after the pool correction               -- found BY this script

The third is the reason it exists rather than remaining a TODO: the home page computes
"lineage-specific candidates" by summing ps_rank for rank >= 18, so it was showing the
WITHDRAWN 1,725 while every other page said 1,611.

This recomputes from the shards -- the shards are the artifact, everything else is a summary of
them -- and prints a diff. Run it after ANY patcher.
"""
import json, os, sys
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")


def shards():
    for hap in ("hap1", "hap2"):
        for dp, _, fs in os.walk(os.path.join(DATA, "genes", hap)):
            for fn in sorted(fs):
                if fn.endswith(".json"):
                    for gid, rec in json.load(open(os.path.join(dp, fn))).items():
                        yield hap, gid, rec


def main():
    man = json.load(open(os.path.join(DATA, "meta", "manifest.json")))
    n = 0
    by_chr, ps_rank, syn_cls = Counter(), Counter(), Counter()
    changed = unstable = n_allele = 0
    for hap, gid, rec in shards():
        n += 1
        by_chr[f"{hap}/{rec['chr']}"] += 1
        ps = rec.get("ps") or {}
        ps_rank[str(ps["rank"]) if ps.get("rank") is not None else "None"] += 1
        syn = rec.get("syn") or {}
        if syn.get("cls"):
            syn_cls[syn["cls"]] += 1
        if syn.get("changed"):
            changed += 1
        if syn.get("stable") is False:
            unstable += 1
        if rec.get("allele"):
            n_allele += 1

    want = {
        "n_genes": (man.get("n_genes"), n),
        "summary.by_chr": (man.get("summary", {}).get("by_chr"), dict(by_chr)),
        "summary.ps_rank": (man.get("summary", {}).get("ps_rank"), dict(ps_rank)),
        "summary.synteny_class": (man.get("summary", {}).get("synteny_class"), dict(syn_cls)),
        "summary.changed": (man.get("summary", {}).get("changed"), changed),
        "summary.unstable": (man.get("summary", {}).get("unstable"), unstable),
        # genes carrying an allele partner; the home page shows half of it as "allele pairs".
        # Stayed at the old map's 51,862 for a release after the allele switch (2026-10-03).
        "annotation.n_allele": (man.get("annotation", {}).get("n_allele"), n_allele),
    }
    stale = []
    for key, (shipped, recomputed) in want.items():
        if shipped == recomputed:
            print(f"  OK    {key}")
            continue
        stale.append(key)
        print(f"  STALE {key}")
        if isinstance(shipped, dict) and isinstance(recomputed, dict):
            for k in sorted(set(shipped) | set(recomputed), key=str):
                a, b = shipped.get(k), recomputed.get(k)
                if a != b:
                    print(f"          {k}: manifest {a} -> shards {b}")
        else:
            print(f"          manifest {shipped} -> shards {recomputed}")

    # the number the home page actually renders, spelled out, because that is what a reader sees
    sh = man.get("summary", {}).get("ps_rank") or {}
    home_shipped = sum(v for k, v in sh.items() if k != "None" and int(k) >= 18)
    home_true = sum(v for k, v in ps_rank.items() if k != "None" and int(k) >= 18)
    print(f"\nhome page 'lineage-specific candidates (PS >= 18)': renders {home_shipped}, "
          f"shards say {home_true}"
          + ("" if home_shipped == home_true else "   <<< the home page contradicts every other page"))
    lp = (man.get("lsg_pool") or {}).get("pool", {})
    print(f"manifest lsg_pool.pool (the genEra rerun that added Idesia): HAP1 {lp.get('hap1')} / HAP2 "
          f"{lp.get('hap2')}")

    if stale:
        print(f"\n{len(stale)} STALE aggregate(s): {stale}")
        return 1
    print("\nevery shard-derivable manifest aggregate matches the shards")
    return 0


if __name__ == "__main__":
    sys.exit(main())
