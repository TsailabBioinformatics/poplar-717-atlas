#!/usr/bin/env python3
"""The gene's founder family (`ps.fnd`) from the same genEra run as its rank and the LSG pool.

build_phylostrat.py writes `ps.fnd` from the frozen 2026-04 run. Since C5 the pool and `ps.rank`
come from the rerun that added Idesia, and joining the two runs made a family look younger than
its own member on 419 genes, which cannot happen inside one run. This replaces `ps.fnd` with the
rerun's own founder families (sources/founder_families_rerun.tsv.gz, from ms2-lsg; see its
provenance record), so the two fields are one run again.

WRITES, per gene: ps.fnd = {rank, fam} (the family's origin rank and size), or no ps.fnd where the
rerun has none: 39 genes whose family has no age row and 6 with no family. The card says "not
available" for those; nothing is filled in.

GATES, enforced rather than trusted:
  1. The source's gene_rank equals the ps.rank this atlas serves on every gene that has one, so
     the rank on the card and the family beside it are the same run. Run patch_lsg_pool.py first.
  2. Afterwards, no gene carries a family younger than itself (fnd.rank <= ps.rank).
  3. The genes left without ps.fnd are exactly the genes with no rank.
"""
import csv
import gzip
import io
import json
import os
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
SRC = os.path.join(ROOT, "sources", "founder_families_rerun.tsv.gz")


def rk(v):
    return None if v in ("", "nan", "NA", None) else int(float(v))


def shard_paths(root):
    for dp, _, fs in os.walk(root):
        for fn in sorted(fs):
            if fn.endswith(".json"):
                yield os.path.join(dp, fn)


def main():
    with gzip.open(SRC, "rt") as fh:
        src = {r["gene_id"]: r for r in csv.DictReader(io.StringIO(fh.read()), delimiter="\t")}
    assert {"family_id", "gene"}.isdisjoint(next(iter(src.values()))), "source carries dropped columns"

    paths = [p for h in ("hap1", "hap2") for p in shard_paths(os.path.join(DATA, "genes", h))]
    off, young, missing, unranked, n = [], [], set(), set(), Counter()
    changed = 0
    for p in paths:
        shard = json.load(open(p))
        for gid, rec in shard.items():
            ps = rec.setdefault("ps", {})
            if ps.get("rank") is None:
                unranked.add(gid)
            r = src.get(gid)
            if r is None or rk(r["gene_rank"]) is None:
                missing.add(gid)
                changed += "fnd" in ps
                ps.pop("fnd", None)
                continue
            if rk(r["gene_rank"]) != ps.get("rank"):
                off.append(gid)
            new = {"rank": rk(r["family_rank"]), "fam": int(r["family_size"])}
            changed += ps.get("fnd") != new
            ps["fnd"] = new
            if new["rank"] > (ps.get("rank") if ps.get("rank") is not None else 99):
                young.append(gid)
            n["origin older than gene" if new["rank"] < ps["rank"] else "same stratum"] += 1
        with open(p, "w") as fh:
            json.dump(shard, fh, separators=(",", ":"), sort_keys=True)

    assert not off, (f"GATE 1 FAILED: {len(off)} genes whose family table rank is not the served "
                     f"ps.rank, e.g. {off[:5]}. Run patch_lsg_pool.py first.")
    assert not young, f"GATE 2 FAILED: {len(young)} genes in a family younger than themselves: {young[:5]}"
    assert missing == unranked, (f"GATE 3 FAILED: genes without a shown family {len(missing)} vs unranked "
                                 f"{len(unranked)}; differ on {sorted(missing ^ unranked)[:5]}")
    no_row = sum(1 for g in missing if g not in src)
    print(f"ps.fnd from the pool's run on {sum(n.values())} genes ({dict(n)}); "
          f"{changed} records changed; not shown on {len(missing)}, all unranked: "
          f"{no_row} with no family assignment, {len(missing) - no_row} with a family but no gene age")


if __name__ == "__main__":
    main()
