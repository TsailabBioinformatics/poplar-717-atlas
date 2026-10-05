#!/usr/bin/env python3
"""TODO 30: chromatin at the P. trichocarpa ortholog's promoter (GEO GSE128434).

This is a property OF THE ORTHOLOG. The measurement was made in P. trichocarpa leaf and
scored against P. trichocarpa gene starts; no chromatin data exists for 717 or for either
parent species. It attaches to a 717 gene only through the atlas's reciprocal 1:1 ortholog,
and the card says so in words. Reading it as this gene's own promoter would assert promoter
synteny that a reciprocal-ortholog call does not establish.

Attached ONLY where `ptri.one2one == 1`. An ambiguous ortholog is not a licence to pick a
top hit, and this script asserts it never attaches to one.

`ex` (expressed in leaf) is stored first and gates the rest on the card, because the
chromatin is leaf only: a gene silent in leaf correctly has no promoter mark there, so an
absent mark is uninterpretable without it. Same discipline as the proteomics detectability
column, and the source analysis's first run failed its own gate for precisely this reason.
"""
import csv
import gzip
import io
import json
import os
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
SRC = os.path.join(ROOT, "sources", "chromatin_ptri.tsv.gz")


def shard_paths(root):
    for dp, _, fs in os.walk(root):
        for fn in sorted(fs):
            if fn.endswith(".json"):
                yield os.path.join(dp, fn)


def main():
    with gzip.open(SRC, "rt") as fh:
        rows = {r["locus"]: r for r in csv.DictReader(io.StringIO(fh.read()), delimiter="\t")}
    assert "rank" not in next(iter(rows.values())), (
        "the P. trichocarpa `rank` column is present in the shipped source. It would be read "
        "as the 717 rank on a gene page -- make_chromatin_source.py must drop it.")
    print(f"source: {len(rows)} P. trichocarpa loci")

    attached, no_ortho, ambiguous_skipped, stats = 0, 0, 0, Counter()
    used = set()
    for hap in ("hap1", "hap2"):
        for p in shard_paths(os.path.join(DATA, "genes", hap)):
            shard = json.load(open(p))
            changed = False
            for gid, rec in shard.items():
                pt = rec.get("ptri") or {}
                r = rows.get(pt["id"]) if pt.get("one2one") == 1 and pt.get("id") else None
                if r is None:
                    if rec.pop("atac", None) is not None:
                        changed = True
                    if pt and pt.get("one2one") != 1:
                        ambiguous_skipped += 1
                    else:
                        no_ortho += 1
                    continue
                ex = r["expressed"] == "1"
                a = {"ex": 1 if ex else 0}
                if ex:
                    a["rpk"] = round(float(r["rpk"]), 2)
                    a["acr"] = 1 if r["acr"] == "1" else 0
                    a["k4"] = 1 if r["h3k4me3"] == "1" else 0
                    a["k36"] = 1 if r["h3k36me3_body"] == "1" else 0
                    a["own"] = 1 if r["own_promoter"] == "1" else 0
                    for k in ("acr", "k4", "k36", "own"):
                        stats[k] += a[k]
                    stats["expressed"] += 1
                else:
                    stats["silent_in_leaf"] += 1
                rec["atac"] = a
                used.add(pt["id"])
                attached += 1
                changed = True
            if changed:
                with open(p, "w") as fh:
                    json.dump(shard, fh, separators=(",", ":"), sort_keys=True)

    print(f"attached to {attached} genes via {len(used)} distinct P. trichocarpa loci")
    print(f"  no 1:1 ortholog: {no_ortho}   ortholog ambiguous, deliberately skipped: "
          f"{ambiguous_skipped}")
    print(f"  of the attached: {stats['expressed']} leaf-expressed, "
          f"{stats['silent_in_leaf']} silent in leaf")
    print(f"  among leaf-expressed: acr {stats['acr']}, H3K4me3 {stats['k4']}, "
          f"own promoter {stats['own']}, H3K36me3 body {stats['k36']}")

    # a phased diploid should give at most two atlas genes per trichocarpa locus
    per_locus = Counter()
    for hap in ("hap1", "hap2"):
        for p in shard_paths(os.path.join(DATA, "genes", hap)):
            for gid, rec in json.load(open(p)).items():
                pt = rec.get("ptri") or {}
                if rec.get("atac") and pt.get("id"):
                    per_locus[pt["id"]] += 1
    worst = per_locus.most_common(1)
    assert not worst or worst[0][1] <= 2, (
        f"{worst[0][0]} carries chromatin for {worst[0][1]} atlas genes; a reciprocal 1:1 "
        "ortholog in a phased diploid should reach at most one gene per haplotype")
    print(f"  loci reaching both haplotypes: {sum(1 for v in per_locus.values() if v == 2)}")

    man_path = os.path.join(DATA, "meta", "manifest.json")
    man = json.load(open(man_path))
    man["chromatin"] = {
        "source": "GEO GSE128434 (Lu, Marand, Schmitz et al., Nature Plants 2019), "
                  "P. trichocarpa leaf: ATAC-seq, H3K4me3, H3K36me3, and matched RNA-seq",
        "scope": "A property OF THE P. TRICHOCARPA ORTHOLOG, not of the 717 gene. No "
                 "chromatin data exists for 717 or for either parent species. Attached only "
                 "where the ortholog call is reciprocal 1:1.",
        "n_genes": attached,
        "n_loci": len(used),
        "genes_without_a_1to1_ortholog": no_ortho,
        "genes_whose_ortholog_is_ambiguous": ambiguous_skipped,
        "leaf_expressed": stats["expressed"],
        "silent_in_leaf": stats["silent_in_leaf"],
        "why_expression_gates_it": "The chromatin is leaf only. A gene silent in leaf "
                                   "correctly carries no promoter mark there, so an absent "
                                   "mark says nothing about the gene unless you know it was "
                                   "expressed. The source analysis's first run failed its own "
                                   "gate for exactly this reason.",
        "the_standing_limit": "ATAC and H3K4me3 establish THAT a locus has a promoter, never "
                              "WHERE transcription starts. Only a 5'-end method does that, "
                              "and none exists for any Populus.",
        "accessible_region_definition": "genic or within 2 kb of the start, replicates merged "
                                        "-- the source paper's own classification. A 500 bp "
                                        "window fails, because 55% of conserved gene starts "
                                        "have their nearest accessible region more than 2 kb "
                                        "away, which is that paper's headline finding about "
                                        "distal elements in plants.",
    }
    with open(man_path, "w") as fh:
        json.dump(man, fh, indent=1, sort_keys=True)
    print("manifest: chromatin block written")


if __name__ == "__main__":
    main()
