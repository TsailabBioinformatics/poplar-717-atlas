#!/usr/bin/env python3
"""Write each gene's place in the ancestral eudicot karyotype onto its record as `aek`.

  aek = {c:  ancestral chromosome, 1-7 (AEK1 to AEK7)
         b:  1 = inside an AEK collinear block, 0 = chromosome inherited from neighbouring anchors
         g:  gamma copy, 1-3      s: salicoid copy, 1-2
         sp: share of the 10 nearest collinear anchors agreeing with that label
         x:  1 where the AEK block and the P. trichocarpa naming disagree (so no g/s/sp)
         an: [{id: AEK gene, sl: [[gamma copy, salicoid copy, [genes]], ...]}]   direct anchors only}

An absent key is an absent value: no g/s/sp means no copy was assigned, and no `aek` at all means
the gene was not placed (or, for the 131 HAP1 unplaced-scaffold genes, not assessed). `sl` lists
only the OCCUPIED slots of that AEK gene, genes of the same haplotype only; every anchor of one AEK
gene carries the identical table, so a gene page needs no second fetch. Exploratory layer, see
make_karyotype_source.py.

Reads only the committed sources, so it runs on any machine and is a `repo` stage in build_all.sh.
Removes `aek` from every gene the source does not place, so a changed source leaves nothing stale.
Run build_facets.py afterwards: the facet index carries an `aek` column.
"""
import csv, gzip, json, os, sys
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
SRC = os.path.join(ROOT, "sources", "karyotype_route_c.tsv.gz")
SRC_ANCHORS = os.path.join(ROOT, "sources", "karyotype_route_c_anchors.tsv.gz")


def rows(path):
    with gzip.open(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    rec, assessed = {}, Counter()
    for x in rows(SRC):
        assessed[x["hap"]] += 1
        if x["aek_chromosome"] == "":
            continue
        r = {"c": int(x["aek_chromosome"]), "b": int(x["inside_aek_block"])}
        if x["gamma_copy"] != "":
            r.update(g=int(x["gamma_copy"]), s=int(x["salicoid_copy"]), sp=float(x["label_support"]))
        if x["conflict"] == "1":
            assert "g" not in r, f"{x['gene_id']}: a conflict gene cannot carry a copy"
            r["x"] = 1
        rec[x["gene_id"]] = (x["hap"].lower(), r)

    occupied = defaultdict(lambda: defaultdict(list))     # (hap, AEK gene) -> (gamma, salicoid) -> genes
    anchor_of = defaultdict(set)                          # gene -> {(hap, AEK gene)}
    for a in rows(SRC_ANCHORS):
        gid, hap, slot = a["gene_id"], a["hap"].lower(), (int(a["gamma_copy"]), int(a["salicoid_copy"]))
        h, r = rec[gid]
        assert h == hap and (r["c"], r.get("g"), r.get("s")) == (int(a["aek_chromosome"]), *slot), \
            f"{gid}: anchor slot {a} disagrees with the gene's own label {r}"
        occupied[(hap, a["aek_gene"])][slot].append(gid)
        anchor_of[gid].add((hap, a["aek_gene"]))
    table = {k: [[g, s, sorted(v[(g, s)])] for g, s in sorted(v)] for k, v in occupied.items()}
    for gid, keys in anchor_of.items():
        rec[gid][1]["an"] = [{"id": e, "sl": table[(h, e)]} for h, e in sorted(keys)]

    n = defaultdict(Counter)
    for hap in ("hap1", "hap2"):
        for dp, _, fs in os.walk(os.path.join(DATA, "genes", hap)):
            for fn in sorted(fs):
                if not fn.endswith(".json"):
                    continue
                p = os.path.join(dp, fn)
                shard = json.load(open(p))
                for gid, g in shard.items():
                    if gid not in rec:
                        g.pop("aek", None)
                        continue
                    h, r = rec[gid]
                    assert h == hap, f"{gid}: source haplotype {h} vs shard {hap}"
                    g["aek"] = r
                    c = n[hap]
                    c["placed"] += 1
                    c["copy_assigned"] += "g" in r
                    c["inside_aek_block"] += r["b"] == 1
                    c["inherited_from_neighbouring_anchors"] += r["b"] == 0
                    c["conflict_no_copy"] += r.get("x", 0)
                    c["direct_anchors"] += "an" in r
                with open(p, "w") as fh:
                    json.dump(shard, fh, separators=(",", ":"), sort_keys=True)

    written = sum(c["placed"] for c in n.values())
    assert written == len(rec), f"{written} records written vs {len(rec)} placed genes in the source"
    assert sum(c["direct_anchors"] for c in n.values()) == len(anchor_of), "an anchor gene was not written"

    man_path = os.path.join(DATA, "meta", "manifest.json")
    man = json.load(open(man_path))
    man["karyotype"] = {
        "status": "exploratory; the gene-level WGD events remain the duplication call",
        "reference": "Wang et al. 2022, BMC Biology 20:216, doi:10.1186/s12915-022-01420-1",
        "source": "sources/karyotype_route_c.tsv.gz and sources/karyotype_route_c_anchors.tsv.gz (route C)",
        "copy_numbers": "gamma copy 1-3 and salicoid copy 1-2 are completeness ranks within one ancestral "
                        "chromosome from Wang et al.'s karyotype projection, not subgenomes",
        **{hap: {"assessed": assessed[hap.upper()],
                 "not_placed": assessed[hap.upper()] - n[hap]["placed"],
                 "ancestral_genes_with_an_anchor": sum(1 for h, _ in occupied if h == hap),
                 **dict(sorted(n[hap].items()))}
           for hap in ("hap1", "hap2")},
    }
    with open(man_path, "w") as fh:
        json.dump(man, fh, indent=1, sort_keys=True)
    for hap in ("hap1", "hap2"):
        print(f"GATE PASS {hap.upper()}: {json.dumps(man['karyotype'][hap], sort_keys=True)}")


if __name__ == "__main__":
    sys.exit(main())
