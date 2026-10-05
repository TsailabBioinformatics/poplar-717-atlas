#!/usr/bin/env python3
"""Fold the CRE promoter layer into every gene's shard, under key "cre".

Source: sources/cre_promoter_hits.tsv.gz, delivered by the ms2-lsg agent (see
sources/cre_promoter_hits.provenance.json for the full derivation). Reviewed and
verified before this script was written: 63,960 rows match the atlas gene set
exactly (0 missing either direction), family-sum equals cre_total on every row,
tata_core rate over the LSG pool reproduces the dissertation chapter's own
promoter_ness.tsv to six decimals on both haplotypes.

Deliberate choices, made HERE rather than left to the UI to paper over:

  - cre_LFY is dropped. It is identically zero for all 63,960 genes (one 19bp
    matrix that never clears the length-biased threshold anywhere in the
    genome) -- the source's own provenance says this should not be rendered as
    a lane, so it is not carried at all.
  - The two family keys that read as cis-element names in the source data
    (ABRE_bZIP, DREB_ERF) are NOT the columns' actual meaning: they are
    TF-FAMILY groupings by JASPAR matrix name, not confirmed ABRE/DRE cis-
    elements. Renamed at ingestion to bZIP and AP2_ERF -- the family names,
    not the element names -- so nothing downstream can accidentally display
    "ABRE" as if a specific cis-regulatory element were identified. The
    original tsv/provenance keep their own names for traceability back to the
    chapter's s08_promoter_cre.py, which this ingestion does not touch.
  - matrices_per_family (how many JASPAR matrices define each family, wildly
    uneven: 1 for LFY, 73 for AP2/ERF) is genome-constant, so it goes in the
    manifest once rather than being repeated 63,960 times.

Per gene, "cre" holds only what the card needs, ~11 small numbers:
  {gc, tata: {core: bool, pos: int|null}, ov: int, nf: float, fam: {...}, tot: int}
"tot" excludes LFY (dropped above), so it is the sum of the KEPT family columns,
not the source file's cre_total. Documented in the manifest, not just here.
"""
import csv, gzip, json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
SRC = os.path.join(ROOT, "sources", "cre_promoter_hits.tsv.gz")
PROV = os.path.join(ROOT, "sources", "cre_promoter_hits.provenance.json")

# source column -> (shard key, provenance key, display label). Order is display order.
# ABRE_bZIP / DREB_ERF renamed per the docstring above; everything else already
# names a TF family, not a cis-element, so it is carried through unchanged.
FAMILY_MAP = [
    ("cre_ABRE_bZIP", "ABRE_bZIP", "bZIP", "bZIP"),
    ("cre_DREB_ERF", "DREB_ERF", "AP2_ERF", "AP2/ERF"),
    ("cre_WRKY", "WRKY", "WRKY", "WRKY"),
    ("cre_MYB", "MYB", "MYB", "MYB"),
    ("cre_NAC", "NAC", "NAC", "NAC"),
    ("cre_MADS_box", "MADS_box", "MADS", "MADS-box"),
    ("cre_SPL", "SPL", "SPL", "SPL"),
    ("cre_TCP", "TCP", "TCP", "TCP"),
    ("cre_ARF", "ARF", "ARF", "ARF"),
    ("cre_other", "other", "other", "other"),
]
DROPPED_FAMILY = "cre_LFY"   # identically zero genome-wide; see docstring


def shard_paths(root):
    out = []
    for dp, _, fs in os.walk(root):
        out += [os.path.join(dp, f) for f in fs if f.endswith(".json")]
    return out


def main():
    if not os.path.exists(SRC):
        sys.exit(f"CRE source not found: {SRC}\nExpected the ms2-lsg deliverable to be present.")

    prov = json.load(open(PROV))
    print(f"provenance: {prov['rows_written']} rows, JASPAR {prov['motifs']['n_matrices']} "
          f"matrices, window {prov['promoter_window']['shipped_in_this_file_bp']} bp")

    cre = {}
    with gzip.open(SRC, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            # Absent since 2026-09-06: the upstream contract now drops all-zero
            # family columns at source (Chen's call), so cre_LFY is no longer
            # shipped. Older copies still carry it as all-zero, so accept both
            # and keep asserting it is zero whenever it IS present.
            assert row.pop(DROPPED_FAMILY, "0") == "0", f"{row['gene_id']}: {DROPPED_FAMILY} not zero"
            fam = {key: int(row[col]) for col, _prov_key, key, _label in FAMILY_MAP}
            tot = sum(fam.values())
            cre[row["gene_id"]] = {
                "gc": round(float(row["gc"]), 3),
                "ov": int(row["n_upstream_overlap"]),
                "nf": round(float(row["n_frac"]), 3),
                "tata": {"core": row["tata_core"] == "True",
                        "pos": int(row["tata_pos"]) if row["tata_pos"] else None},
                "fam": fam,
                "tot": tot,
            }
    print(f"parsed {len(cre):,} genes")

    patched = 0
    for dp, _, fs in os.walk(os.path.join(DATA, "genes")):
        for fn in fs:
            if not fn.endswith(".json"):
                continue
            p = os.path.join(dp, fn)
            shard = json.load(open(p))
            changed = False
            for gid, rec in shard.items():
                if gid in cre:
                    rec["cre"] = cre[gid]
                    changed = True
                    patched += 1
            if changed:
                with open(p, "w") as fh:
                    json.dump(shard, fh, separators=(",", ":"), sort_keys=True)
    print(f"patched {patched:,} gene records with cre")
    assert patched == len(cre), f"expected to patch {len(cre)}, patched {patched}"

    man_path = os.path.join(DATA, "meta", "manifest.json")
    man = json.load(open(man_path))
    man.setdefault("modules", [])
    if "promoter" not in man["modules"]:
        man["modules"].append("promoter")
    man["cre"] = {
        "n_genes": patched,
        "tier": "public",
        "source": "ms2-lsg pipeline/s08d_promoter_cre_genomewide.py, reviewed 2026-09-06",
        "jaspar": {"database": prov["motifs"]["database"], "n_matrices": prov["motifs"]["n_matrices"],
                  "sha256_16": prov["motifs"]["sha256"][:16]},
        "hit_definition": prov["hit_definition"]["rule"],
        "hit_threshold_frac": prov["hit_definition"]["HIT_FRAC"],
        "promoter_window_bp": prov["promoter_window"]["shipped_in_this_file_bp"],
        "promoter_scanned_bp": prov["promoter_window"]["scanned_bp"],
        "family_labels": {key: label for _col, _prov_key, key, label in FAMILY_MAP},
        # matrix counts per KEPT family, from the provenance's matrices_per_family, for
        # normalizing the card's bars -- a family with 73 matrices will out-count a family
        # with 17 on raw hits alone, which is the source's own documented caution.
        "matrices_per_family": {
            key: prov["family_mapping"]["matrices_per_family"][prov_key]
            for _col, prov_key, key, _label in FAMILY_MAP
        },
        "dropped_family_lfy": "identically zero genome-wide (n_matrices=1, length-biased "
                              "threshold), see build_cre.py docstring",
        "not_included": "per-matrix enrichment statistics and LSG-vs-control comparisons; "
                        "those are population-level chapter results, not per-gene facts",
        "caveat_family_size_uneven": True,
    }
    man["sources"][SRC.replace(ROOT + "/", "")] = {
        "sha256_16": "see sources/cre_promoter_hits.provenance.json", "rows": prov["rows_written"]}
    with open(man_path, "w") as fh:
        json.dump(man, fh, indent=1, sort_keys=True)
    print("manifest updated: cre block + matrices_per_family (family-size caveat carried through)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
