#!/usr/bin/env python3
"""Attach the route C combined Ks + topology WGD call to each gene, as `dup.cc`.

`dup.s` / `dup.a` (the canonical v9 letters) are NOT touched. Both calls ship, because they
disagree on 344 HAP1 and 413 HAP2 genes and the atlas should not be the place that decides
which is right -- MS1's own record still lists v9 as canonical and marks the combined call
"nothing applied, for the WGD track owner".

  dup.cc = {"s": 0|1, "a": 0|1, "r": 0|1}   present only where the gene is in an anchor pair
  dup.ccd = 1                               set when cc disagrees with the v9 letters

`ccd` exists so the page can lead with the disagreement instead of asking a reader to compare
two rows themselves, and so Browse can offer "genes whose WGD call is disputed" without
recomputing it client-side.
"""
import csv, gzip, io, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
SRC = os.path.join(ROOT, "sources", "wgd_combined_call.tsv.gz")


def shard_paths(root):
    for dp, _, fs in os.walk(root):
        for fn in sorted(fs):
            if fn.endswith(".json"):
                yield os.path.join(dp, fn)


def main():
    with gzip.open(SRC, "rt") as fh:
        cc = {}
        for r in csv.DictReader(io.StringIO(fh.read()), delimiter="\t"):
            cc[r["gene_id"]] = (int(r["cc_salicoid"]), int(r["cc_gamma"]), int(r["cc_review"]))
    print(f"{len(cc):,} genes with a combined call")

    n_set, n_dis, moved = 0, 0, {"into_salicoid": 0, "out_of_salicoid": 0,
                                 "into_gamma": 0, "out_of_gamma": 0}
    for hap in ("hap1", "hap2"):
        for p in shard_paths(os.path.join(DATA, "genes", hap)):
            shard = json.load(open(p))
            changed = False
            for gid, rec in shard.items():
                d = rec.get("dup")
                if not d:
                    continue
                v = cc.get(gid)
                if not v:
                    if d.pop("cc", None) is not None or d.pop("ccd", None) is not None:
                        changed = True
                    continue
                s, a, rev = v
                d["cc"] = {"s": s, "a": a, "r": rev}
                n_set += 1
                v9s, v9a = int(bool(d.get("s"))), int(bool(d.get("a")))
                if s != v9s or a != v9a:
                    d["ccd"] = 1
                    n_dis += 1
                    if s and not v9s: moved["into_salicoid"] += 1
                    if v9s and not s: moved["out_of_salicoid"] += 1
                    if a and not v9a: moved["into_gamma"] += 1
                    if v9a and not a: moved["out_of_gamma"] += 1
                else:
                    d.pop("ccd", None)
                changed = True
            if changed:
                with open(p, "w") as fh:
                    json.dump(shard, fh, separators=(",", ":"), sort_keys=True)

    print(f"{n_set:,} genes carry a combined call; {n_dis:,} disagree with the v9 letters")
    print(f"   {moved}")

    man_path = os.path.join(DATA, "meta", "manifest.json")
    man = json.load(open(man_path))
    man["wgd_combined_call"] = {
        "genes_with_a_call": n_set,
        "genes_disagreeing_with_v9": n_dis,
        "movement": moved,
        "rule": "topology where grape and segment Ks independently agree; v8 otherwise; "
                "split-evidence pairs set aside as `review`.",
        "canonical": "master_duplication_table_v9 (v8 rule 5c + dS resolution) remains the "
                     "atlas's called class. The combined call is shown beside it as a "
                     "proposed correction, NOT applied upstream as of 2026-09-16.",
    }
    json.dump(man, open(man_path, "w"), indent=1)


if __name__ == "__main__":
    main()
