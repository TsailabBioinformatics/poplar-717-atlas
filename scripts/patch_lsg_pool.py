#!/usr/bin/env python3
"""The LSG pool is the genEra rerun that added Idesia; the frozen call is kept as provenance.

RULING (Chen, 2026-09-21 evening, reconfirmed 2026-09-26 to two sessions independently; ms2-lsg
docs/CORRECTIONS.md C5, commit 8eb65b0): "make the new canon with the idesia included run as the
canon, and only AbSENSE and synteny check are the QC steps". The pool is 776 HAP1 + 835 HAP2,
rank >= 18 in the rerun that added Idesia polycarpa (ranks 18/19: 433/343 and 497/338).

HISTORY. Named for what it patches rather than for which call is the pool, because that has
flipped three times: this file was patch_pool_idesia.py (2026-09-10, C2: the rerun is the pool,
frozen call under `ps.frz`), then patch_pool_frozen.py (2026-09-21 morning, C4: the frozen call
is the pool, rerun under `ps.rr`). C5 reverses C4, so the 2026-09-10 arrangement is back: `ps` is
the rerun's call and the frozen call sits under `ps.frz` wherever its rank differs.

The frozen call stays on every gene it differs on, as provenance only: "frozen 2026-04 genEra
call, superseded". Nothing reads it for a filter, a count or a default view. Dropping it would
hide why 114 genes left the pool, and a reader holding an older figure could not tell.

WHAT IT WRITES, per gene:
  ps.rank, ps.name, ps.lsg, ps.tr   the rerun's call (sources/lsg_pool_idesia.tsv.gz)
  ps.frz = {rank, name, lsg}        the frozen call, only where its rank differs (766 genes)
  ps.rr                             removed (it was the rerun, now primary)
and the manifest aggregates that summarise those ranks, recomputed from the same counts.

GATES, enforced rather than trusted:
  1. The frozen source must reproduce the frozen rank the atlas is serving for every gene, in
     any of the states a shard can be in: freshly built or patched by patch_pool_frozen.py
     (ps.rank is frozen), or already patched by this one (ps.frz holds it where it differs).
  2. The rerun's frozen columns equal the frozen source, gene by gene.
  3. The rerun never moves a gene INTO the pool (upstream's premise; asserted per row).
  4. No candidate-label column in either source.
  5. The pool is exactly 776 HAP1 + 835 HAP2, split PS18/PS19 433/343 and 497/338, and is_lsg is
     rank >= 18 in the rerun on every gene. check_ms2_canon.py then compares the membership with
     ms2-lsg's own flag table, gene by gene.

AFTER THIS RUNS, patch_lsg_qc.py must run (it flags the pool genes this defines), then
patch_proteomics_v2.py, patch_map_overlays.py and build_facets.py (all read ps.lsg / ps.rank).
"""
import csv
import gzip
import io
import json
import os
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FROZEN = os.path.join(ROOT, "sources", "lsg_pool_frozen.tsv.gz")
RERUN = os.path.join(ROOT, "sources", "lsg_pool_idesia.tsv.gz")

FORBIDDEN_COLUMNS = {"canonical_label", "group"}
WANT_POOL = {"HAP1": 776, "HAP2": 835}
WANT_STRATA = {"HAP1_PS18": 433, "HAP1_PS19": 343, "HAP2_PS18": 497, "HAP2_PS19": 338}
WANT_FROZEN_POOL = {"HAP1": 842, "HAP2": 883}


def tf(v):
    return str(v).strip().lower() in ("true", "1")


def iv(v):
    v = (v or "").strip()
    return None if v in ("", "NA", "nan", "None") else int(float(v))


def fv(v):
    v = (v or "").strip()
    return None if v in ("", "NA", "nan", "None") else float(v)


def read(path):
    with gzip.open(path, "rt") as fh:
        rows = list(csv.DictReader(io.StringIO(fh.read()), delimiter="\t"))
    leaked = FORBIDDEN_COLUMNS & set(rows[0])
    assert not leaked, f"{os.path.basename(path)} carries embargoed column(s) {leaked}"
    return {r["gene_id"]: r for r in rows}


def shard_paths(root):
    for dp, _, fs in os.walk(root):
        for fn in sorted(fs):
            if fn.endswith(".json"):
                yield os.path.join(dp, fn)


def main():
    frozen, rerun = read(FROZEN), read(RERUN)
    assert set(frozen) == set(rerun), "the frozen and rerun sources cover different genes"
    for g, r in rerun.items():
        assert not tf(r["entered_pool_under_idesia"]), f"{g} claims to ENTER the pool in the rerun"
        assert iv(r["frozen_rank"]) == iv(frozen[g]["rank"]), f"{g}: the two sources disagree"
        assert tf(r["frozen_is_lsg"]) == tf(frozen[g]["is_lsg"]), f"{g}: frozen membership differs"
    print(f"sources: {len(frozen)} genes, both free of {sorted(FORBIDDEN_COLUMNS)}")

    paths = list(shard_paths(os.path.join(DATA, "genes", "hap1"))) + \
            list(shard_paths(os.path.join(DATA, "genes", "hap2")))

    # ---- gate 1: the frozen rank the atlas serves, in whichever state the shard is in ------
    live = {}
    for p in paths:
        for gid, rec in json.load(open(p)).items():
            ps = rec.get("ps") or {}
            live[gid] = ps["frz"]["rank"] if ps.get("frz") else ps.get("rank")
    assert set(live) == set(frozen), (
        f"gene sets differ: {len(set(live) - set(frozen))} only in atlas, "
        f"{len(set(frozen) - set(live))} only in source")
    # Already patched by this script: a gene without ps.frz serves the rerun's rank, which equals
    # the frozen one by construction. Before this script: ps.rank is the frozen rank.
    bad = [(g, live[g], iv(frozen[g]["rank"])) for g in live if live[g] != iv(frozen[g]["rank"])]
    assert not bad, (f"GATE FAILED: {len(bad)} genes where the frozen source is not the frozen rank "
                     f"this atlas serves, e.g. {bad[:5]}. Refusing to write.")
    print(f"GATE PASS: the frozen source reproduces the served frozen rank for all {len(live)} genes")

    # ---- patch ---------------------------------------------------------------------------
    pool, strata, moved, transitions = Counter(), Counter(), 0, Counter()
    for p in paths:
        shard = json.load(open(p))
        for gid, rec in shard.items():
            f, r = frozen[gid], rerun[gid]
            ps = rec.setdefault("ps", {})
            rank, fz_rank = iv(r["rank"]), iv(f["rank"])
            ps.pop("rr", None)
            ps["rank"] = rank
            ps["name"] = r["phylostratum"] or None
            ps["lsg"] = tf(r["is_lsg"])
            tr = fv(r["tax_rep"])
            # an unplaced gene carries no representativeness, as build_data.py always wrote it
            if tr is not None and rank is not None:
                ps["tr"] = tr
            else:
                ps.pop("tr", None)
            if fz_rank != rank:
                ps["frz"] = {"rank": fz_rank, "name": f["phylostratum"] or None,
                             "lsg": tf(f["is_lsg"])}
                moved += 1
                transitions[f"{fz_rank}->{rank}"] += 1
            else:
                ps.pop("frz", None)
            if ps["lsg"]:
                pool[r["hap"]] += 1
                strata[f"{r['hap']}_PS{rank}"] += 1
        with open(p, "w") as fh:
            json.dump(shard, fh, separators=(",", ":"), sort_keys=True)

    off = [g for g, r in rerun.items() if tf(r["is_lsg"]) != ((iv(r["rank"]) or -1) >= 18)]
    assert not off, f"{len(off)} genes where is_lsg disagrees with rank >= 18, e.g. {off[:5]}"
    assert dict(pool) == WANT_POOL, f"pool is {dict(pool)}, expected {WANT_POOL}"
    assert dict(strata) == WANT_STRATA, f"PS18/PS19 split is {dict(strata)}, expected {WANT_STRATA}"
    fz_pool = Counter(f["hap"] for f in frozen.values() if tf(f["is_lsg"]))
    fz_strata = Counter(f"{f['hap']}_PS{iv(f['rank'])}" for f in frozen.values() if tf(f["is_lsg"]))
    left = Counter(r["hap"] for g, r in rerun.items()
                   if tf(frozen[g]["is_lsg"]) and not tf(r["is_lsg"]))
    assert dict(fz_pool) == WANT_FROZEN_POOL, f"frozen pool is {dict(fz_pool)}"
    print(f"patched {len(rerun)} genes; the frozen call gives a different rank for {moved}")
    print(f"  pool (the rerun): {dict(pool)}  PS18/PS19: {dict(sorted(strata.items()))}")
    print(f"  frozen call, superseded: {dict(fz_pool)}; left the pool: {dict(left)}")
    print(f"  largest rank changes, frozen -> rerun: {dict(transitions.most_common(5))}")

    # ---- manifest: every aggregate over these ranks, recomputed from the same counts --------
    man_path = os.path.join(DATA, "meta", "manifest.json")
    man = json.load(open(man_path))
    by_rank, names = {"hap1": Counter(), "hap2": Counter()}, {}
    for g, r in rerun.items():
        hap = "hap1" if r["hap"] == "HAP1" else "hap2"
        rk = iv(r["rank"])
        by_rank[hap]["NA" if rk is None else str(rk)] += 1
        if rk is not None and r["phylostratum"]:
            names[rk] = r["phylostratum"]
    ps_man = man.setdefault("phylostrat", {})
    ps_man["by_rank"] = {h: dict(sorted(c.items())) for h, c in by_rank.items()}
    ps_man["ladder"] = {h: [{"rank": rk, "name": names[rk], "n": by_rank[h][str(rk)]}
                            for rk in sorted(names)] for h in ("hap1", "hap2")}
    ps_man["note"] = (
        "LSG pool = rank >= 18 in the genEra rerun that added Idesia polycarpa and searched the "
        "full comparison panel (PS18 Populus + PS19 hybrid-only): 776 HAP1 + 835 HAP2 (Chen, "
        "2026-09-21 evening; ms2-lsg CORRECTIONS C5). Ranks here are that call. The frozen "
        "2026-04 call is superseded and kept per gene under ps.frz wherever its rank differs, as "
        "provenance only.")
    # summary.ps_rank is read by the home page; "NA" there is spelled "None" (build_data.py's key).
    merged = Counter()
    for hcnt in by_rank.values():
        for k, v in hcnt.items():
            merged["None" if k == "NA" else k] += v
    man.setdefault("summary", {})["ps_rank"] = dict(sorted(merged.items()))
    man["lsg_pool"] = {
        "pool": {"hap1": pool["HAP1"], "hap2": pool["HAP2"], "total": sum(pool.values()),
                 "ps18_ps19": dict(sorted(strata.items())),
                 "call": "the genEra rerun that added Idesia polycarpa and was the first run to "
                         "search the full comparison panel"},
        "ruling": "Chen, 2026-09-21 evening, reconfirmed 2026-09-26: the pool is the rerun that "
                  "added Idesia, 776/835, and only abSENSE and the synteny check are QC steps; a "
                  "flag never changes the pool. ms2-lsg docs/CORRECTIONS.md C5, commit 8eb65b0.",
        "frozen": {"what": "the frozen 2026-04 genEra call (an NCBI nr search; its configured "
                           "comparison panel was never searched), superseded as the pool",
                   "at_rank_18_or_above": {"hap1": fz_pool["HAP1"], "hap2": fz_pool["HAP2"]},
                   "ps18_ps19": dict(sorted(fz_strata.items())),
                   "left_the_pool": {"hap1": left["HAP1"], "hap2": left["HAP2"]},
                   "genes_with_a_different_rank": moved,
                   "note": "Provenance only, never a filter, count or default view. Numbers "
                           "published against this call are correct statements about it and "
                           "are not retracted. A rank change is not always due to Idesia: the "
                           "rerun was also the first to search the comparison panel at all."},
        "never": "PS18 and PS19 are never pooled without their split.",
    }
    with open(man_path, "w") as fh:
        json.dump(man, fh, indent=1, sort_keys=True)
    print("manifest: phylostrat, summary.ps_rank and lsg_pool rewritten")
    print("\nNEXT: patch_lsg_qc.py, patch_proteomics_v2.py, patch_map_overlays.py, build_facets.py.")


if __name__ == "__main__":
    main()
