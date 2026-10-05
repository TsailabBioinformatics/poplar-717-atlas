#!/usr/bin/env python3
"""Name the duplication partners on the gene page, with their Ks and expression correlation.

The card previously showed a duplication CLASS ("ADS"), an array id and a member count -- so a
reader could see that a gene was in a 7-member tandem array and never learn which seven. It
also showed paralog dS/dN/omega but no expression relationship at all.

  para[i].r   expression correlation with that paralog
  td.mem      the array's OTHER members, each with its correlation

WHAT IS DELIBERATELY NOT SHOWN PER MEMBER: a Ks. `td.ks` is the ARRAY MEDIAN, present on 12,082
of 13,663 genes, and rendering it beside an individual member would read as that pair's
divergence. Paralog pairs do have real per-pair dS and it is already shown there.

Absence of `r` means the pair is not comparable on the ComBat matrix (one partner is outside
its 59,599 genes), NOT that the two are uncorrelated. The card must say which.
"""
import csv, gzip, io, json, os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
SRC = os.path.join(ROOT, "sources", "duplicate_coexpr.tsv.gz")
MAX_MEMBERS = 24          # a 56-member array would otherwise dominate the page


def shard_paths(root):
    for dp, _, fs in os.walk(root):
        for fn in sorted(fs):
            if fn.endswith(".json"):
                yield os.path.join(dp, fn)


def main():
    with gzip.open(SRC, "rt") as fh:
        rr = {}
        for r in csv.DictReader(io.StringIO(fh.read()), delimiter="\t"):
            rr[(r["gene_a"], r["gene_b"])] = float(r["r"])
    def corr(a, b):
        return rr.get((a, b), rr.get((b, a)))
    arrays = json.load(open(os.path.join(DATA, "index", "arrays.json")))
    member_of = {}
    for aid, members in arrays.items():
        for g in members:
            member_of[g] = members
    print(f"{len(rr):,} correlations, {len(arrays):,} arrays")

    n_para_r, n_para, n_mem, n_mem_r, truncated = 0, 0, 0, 0, 0
    for hap in ("hap1", "hap2"):
        for p in shard_paths(os.path.join(DATA, "genes", hap)):
            shard = json.load(open(p))
            changed = False
            for gid, rec in shard.items():
                for q in (rec.get("para") or []):
                    n_para += 1
                    c = corr(gid, q["id"])
                    if c is None:
                        q.pop("r", None)
                    else:
                        q["r"] = c
                        n_para_r += 1
                    changed = True
                members = member_of.get(gid)
                td = rec.get("td")
                if members and td:
                    others = [m for m in members if m != gid]
                    if len(others) > MAX_MEMBERS:
                        td["more"] = len(others) - MAX_MEMBERS
                        others = others[:MAX_MEMBERS]
                        truncated += 1
                    else:
                        td.pop("more", None)
                    mem = []
                    for m in others:
                        c = corr(gid, m)
                        mem.append([m] if c is None else [m, c])
                        n_mem += 1
                        if c is not None:
                            n_mem_r += 1
                    td["mem"] = mem
                    changed = True
            if changed:
                with open(p, "w") as fh:
                    json.dump(shard, fh, separators=(",", ":"), sort_keys=True)

    print(f"paralog entries: {n_para:,}, of which {n_para_r:,} got a correlation "
          f"({100*n_para_r/max(n_para,1):.1f}%)")
    print(f"array member entries: {n_mem:,}, of which {n_mem_r:,} got one "
          f"({100*n_mem_r/max(n_mem,1):.1f}%); {truncated} genes had their list truncated")

    man_path = os.path.join(DATA, "meta", "manifest.json")
    man = json.load(open(man_path))
    man["duplicate_coexpr"] = {
        "pairs_with_a_correlation": len(rr),
        "basis": "Pearson r on log2(ComBat-corrected TPM + 1) over 652 samples -- the same "
                 "matrix the neighbour co-expression (nb) uses, because a 652-sample "
                 "correlation spans 28 studies and is therefore a cross-study comparison.",
        "coverage": "A pair needs BOTH partners in the ComBat matrix (59,599 of 63,960 genes), "
                    "so 7,707 of 32,333 duplicate pairs have no value. Absence means 'not "
                    "comparable', not 'uncorrelated'.",
        "distribution": {"median": 0.310, "p10": -0.101, "p90": 0.803},
        "not_a_divergence_measure": "A descriptive correlation only. No model, no threshold, "
                                    "and no claim about duplicate fate. Nothing upstream "
                                    "computes duplicate-pair expression divergence.",
        "per_member_ks_deliberately_absent": "td.ks is the ARRAY MEDIAN, not a pair value, so "
                                             "it is never rendered beside an individual member.",
        "max_members_listed": MAX_MEMBERS,
    }
    with open(man_path, "w") as fh:
        json.dump(man, fh, indent=1, sort_keys=True)
    print("manifest: duplicate_coexpr block written")


if __name__ == "__main__":
    main()
