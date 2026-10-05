#!/usr/bin/env python3
"""Attach the per-PAIR tandem dS to each gene's array members, as `td.pks`.

A new key rather than a third slot in `td.mem`: mem entries are `[id]` or `[id, r]`, so a
positional dS would have to be written as `[id, null, dS]` wherever the correlation is absent
-- a null that means "not comparable" sitting in a list read positionally by two modules. A
map keyed by partner id cannot be read in the wrong position, and a partner simply missing
from it is unambiguous.

`td.ks` (the array median) is left exactly where it is and keeps its own meaning. The two are
not interchangeable and the page must not present them as the same number.

ABSENCE IS INFORMATIVE AND MUST SURVIVE TO THE PAGE. Only ~half of member entries get a value,
because `td.mem` is all-vs-all within an array while the source holds detected edges that also
passed alignment-coverage QC. `td.pks_n` records how many of THIS gene's members were scored,
so the card can say "3 of 7 pairs scored" instead of rendering four silent blanks.
"""
import csv, gzip, io, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
SRC = os.path.join(ROOT, "sources", "td_pair_ks.tsv.gz")


def shard_paths(root):
    for dp, _, fs in os.walk(root):
        for fn in sorted(fs):
            if fn.endswith(".json"):
                yield os.path.join(dp, fn)


def main():
    with gzip.open(SRC, "rt") as fh:
        pairs = {}
        for r in csv.DictReader(io.StringIO(fh.read()), delimiter="\t"):
            pairs[(r["gene_a"], r["gene_b"])] = float(r["dS"])
    print(f"{len(pairs):,} scored tandem pairs")

    n_gene, n_mem, n_hit, n_cleared = 0, 0, 0, 0
    for hap in ("hap1", "hap2"):
        for p in shard_paths(os.path.join(DATA, "genes", hap)):
            shard = json.load(open(p))
            changed = False
            for gid, rec in shard.items():
                td = rec.get("td")
                if not td:
                    continue
                members = [m[0] for m in (td.get("mem") or [])]
                if not members:
                    # a rerun after the member list shrank must not leave a stale map behind
                    if td.pop("pks", None) is not None or td.pop("pks_n", None) is not None:
                        n_cleared += 1
                        changed = True
                    continue
                n_gene += 1
                pks = {}
                for m in members:
                    n_mem += 1
                    ds = pairs.get((gid, m), pairs.get((m, gid)))
                    if ds is not None:
                        pks[m] = ds
                        n_hit += 1
                if pks:
                    td["pks"] = pks
                    td["pks_n"] = [len(pks), len(members)]
                else:
                    td.pop("pks", None)
                    td["pks_n"] = [0, len(members)]
                changed = True
            if changed:
                with open(p, "w") as fh:
                    json.dump(shard, fh, separators=(",", ":"), sort_keys=True)

    print(f"{n_gene:,} genes with array members; {n_hit:,} of {n_mem:,} member entries "
          f"got a per-pair dS ({100 * n_hit / max(n_mem, 1):.1f}%)")

    man_path = os.path.join(DATA, "meta", "manifest.json")
    man = json.load(open(man_path))
    man["td_pair_ks"] = {
        "scored_pairs": len(pairs),
        "member_entries_with_ds": n_hit,
        "member_entries_total": n_mem,
        "basis": "Per-pair dS from the adopted TD union rebuild (union_pairs_final.csv), "
                 "retained edges passing the source's own Ka/Ks coverage QC, dS > 0.",
        "not_the_array_median": "`td.ks` remains the ARRAY MEDIAN and is a different quantity. "
                                "A member without a `pks` entry was never scored as a pair; it "
                                "does not mean the two copies are identical.",
    }
    json.dump(man, open(man_path, "w"), indent=1)


if __name__ == "__main__":
    main()
