#!/usr/bin/env python3
"""Write every WGD anchor partner onto the gene record as `wgd`: [{id, e, ks, ss, cp, m, c, r}].

e = event (S salicoid, A gamma), ks = dS, ss = synonymous sites, cp = chromosome pair,
m = multiplicon, c = dS conflict resolved upstream, r = expression correlation (absent = n/c).
All pairs, no cap, no filter -- see make_wgd_pairs_source.py. Run build_facets.py afterwards.
"""
import csv, gzip, json, os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")


def num(s):
    return None if s in ("", "nan", None) else float(s)


def main():
    part = defaultdict(list)
    with gzip.open(os.path.join(ROOT, "sources", "wgd_pairs.tsv.gz"), "rt") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for x in rows:
        for g, o in ((x["gene_a"], x["gene_b"]), (x["gene_b"], x["gene_a"])):
            rec = {"id": o, "e": x["event"], "ks": num(x["dS"]), "ss": num(x["syn_sites"]),
                   "cp": x["chrpair"], "m": x["multiplicon"]}
            if x["conflict_resolved"] == "True":
                rec["c"] = 1
            if x["r"] != "":
                rec["r"] = float(x["r"])
            part[g].append(rec)
    n_gene = n_link = 0
    for hap in ("hap1", "hap2"):
        for dp, _, fs in os.walk(os.path.join(DATA, "genes", hap)):
            for fn in sorted(fs):
                if not fn.endswith(".json"):
                    continue
                p = os.path.join(dp, fn)
                shard = json.load(open(p))
                for gid, rec in shard.items():
                    ps = sorted(part.get(gid, []), key=lambda q: (q["e"] != "S", q["ks"] if q["ks"] is not None else 9))
                    if ps:
                        rec["wgd"] = ps; n_gene += 1; n_link += len(ps)
                    else:
                        rec.pop("wgd", None)
                with open(p, "w") as fh:
                    json.dump(shard, fh, separators=(",", ":"), sort_keys=True)
    assert n_link == 2 * len(rows), f"{n_link} links for {len(rows)} pairs"
    print(f"GATE PASS: {len(rows):,} pairs -> {n_gene:,} genes, {n_link:,} links (2 per pair)")


if __name__ == "__main__":
    main()
