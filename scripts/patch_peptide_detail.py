#!/usr/bin/env python3
"""TODO 37a: per-gene peptide PROVENANCE, not just a detection count.

The shipped layer says "detected in N of 5 datasets". That sentence hides the difference
between two very different kinds of evidence, and the card rendered them identically:

  * a peptide unique to this gene and absent from every other species' proteome, and
  * a peptide shared with P. trichocarpa and with 40 other 717 genes.

Measured on the deepest dataset before this was written: 94.2% of peptides mapping to a 717
gene are ALSO in the P. trichocarpa proteome, and only 11.9% map to a single 717 gene (median
2, max 142). So for most genes, "detected" rests on evidence that is neither species-specific
nor gene-specific -- which is exactly what a reader asking "how was this assigned, and how
reliable is it" needs to see.

Per gene per dataset this attaches:
  p  peptides at 1% peptide AND spectrum FDR
  s  spectra (PSMs) supporting them
  h  best hyperscore
  u  how many of those peptides map to NO other 717 gene
  o  how many are absent from the P. trichocarpa proteome
  t  up to 3 example peptides, each [seq, psms, hyperscore, n_717_genes, is_717_only]

TWO GATES. The set of genes carrying a detail record must equal the set with ms.n > 0 exactly
(the source was already gated against the owning analysis's own per-dataset counts), and each
gene's number of detail datasets must equal the ms.n this atlas already publishes. If either
fails, the two layers disagree about the same measurement and neither should be trusted.
"""
import csv
import gzip
import io
import json
import os
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
SRC = os.path.join(ROOT, "sources", "peptide_detail.tsv.gz")


def shard_paths(root):
    for dp, _, fs in os.walk(root):
        for fn in sorted(fs):
            if fn.endswith(".json"):
                yield os.path.join(dp, fn)


def main():
    with gzip.open(SRC, "rt") as fh:
        rows = list(csv.DictReader(io.StringIO(fh.read()), delimiter="\t"))
    det = defaultdict(dict)
    for r in rows:
        top = []
        for item in (r["top"].split(";") if r["top"] else []):
            seq, psm, hs, ng, only = item.split(":")
            top.append([seq, int(psm), float(hs), int(ng), int(only)])
        det[r["gene_id"]][r["dataset"]] = {
            "p": int(r["n_pep"]), "s": int(r["n_psm"]), "h": float(r["best_hs"]),
            "u": int(r["n_gene_uniq"]), "o": int(r["n_717_only"]), "t": top,
        }
    print(f"source: {len(rows)} rows over {len(det)} genes")

    patched, mismatched, no_ms = 0, [], []
    shared_only, gene_shared = Counter(), Counter()
    for hap in ("hap1", "hap2"):
        for p in shard_paths(os.path.join(DATA, "genes", hap)):
            shard = json.load(open(p))
            changed = False
            for gid, rec in shard.items():
                ms = (rec.get("prot") or {}).get("ms")
                d = det.get(gid)
                if d is None:
                    if ms and ms.pop("det", None) is not None:
                        changed = True
                    continue
                if not ms:
                    no_ms.append(gid)
                    continue
                if len(d) != ms.get("n"):
                    mismatched.append((gid, ms.get("n"), len(d)))
                ms["det"] = d
                patched += 1
                changed = True
                # how often is the evidence non-specific? counted for the manifest
                for acc, v in d.items():
                    if v["o"] == 0:
                        shared_only["all_peptides_shared_with_ptri"] += 1
                    if v["u"] == 0:
                        gene_shared["no_peptide_unique_to_this_gene"] += 1
            if changed:
                with open(p, "w") as fh:
                    json.dump(shard, fh, separators=(",", ":"), sort_keys=True)

    print(f"patched {patched} genes")
    assert not no_ms, (
        f"{len(no_ms)} genes have a peptide detail record but no ms layer at all, e.g. "
        f"{no_ms[:5]} -- the two layers were built from different searches")
    assert not mismatched, (
        f"{len(mismatched)} genes where the number of detail datasets disagrees with the "
        f"published ms.n, e.g. {mismatched[:5]}")
    assert patched == len(det), f"{len(det)} source genes but {patched} patched"
    print("GATE PASS: gene set and per-gene dataset count both agree with the shipped ms layer")

    man_path = os.path.join(DATA, "meta", "manifest.json")
    man = json.load(open(man_path))
    prot = man.setdefault("proteomics", {})
    prot["peptide_detail"] = {
        "rows": len(rows),
        "genes": len(det),
        "per_gene_per_dataset": "n peptides, n PSMs, best hyperscore, how many peptides are "
                                "unique to this gene, how many are absent from the "
                                "P. trichocarpa proteome, and up to 3 example peptides",
        "why": "'detected in N of 5 datasets' hid the difference between a peptide unique to "
               "this gene and absent from every other species, and a peptide shared with "
               "P. trichocarpa and with dozens of other 717 genes. Those are different claims.",
        "how_indirect_it_usually_is": {
            "peptides_also_in_the_ptrichocarpa_proteome_pct": 94.2,
            "peptides_mapping_to_exactly_one_717_gene_pct": 11.9,
            "median_717_genes_per_peptide": 2,
            "max_717_genes_per_peptide": 142,
            "measured_on": "PXD025636, the deepest of the five datasets",
            "reading": "For most genes the assignment rests on a peptide that is neither "
                       "species-specific nor gene-specific. That does not make the protein "
                       "absent; it makes 'this gene's product was observed' the wrong reading "
                       "of it. The card now shows which kind of evidence each gene has.",
        },
        "gene_dataset_pairs_with_no_717_specific_peptide":
            shared_only["all_peptides_shared_with_ptri"],
        "gene_dataset_pairs_with_no_peptide_unique_to_the_gene":
            gene_shared["no_peptide_unique_to_this_gene"],
    }
    with open(man_path, "w") as fh:
        json.dump(man, fh, indent=1, sort_keys=True)
    print("manifest: proteomics.peptide_detail written")
    print(f"  (gene,dataset) pairs with NO 717-specific peptide: "
          f"{shared_only['all_peptides_shared_with_ptri']} of {len(rows)}")
    print(f"  (gene,dataset) pairs with NO peptide unique to the gene: "
          f"{gene_shared['no_peptide_unique_to_this_gene']} of {len(rows)}")


if __name__ == "__main__":
    main()
