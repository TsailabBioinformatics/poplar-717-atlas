#!/usr/bin/env python3
"""TODO 28: replace the Tier-1 proteomics string match with the Tier-2 / aspen re-search.

WHAT CHANGES AND WHY. The atlas shipped `prot.ms` as three datasets of Tier-1 evidence:
a 717 protein's tryptic peptides string-matched against peptide tables a third party
produced by searching a DIFFERENT species' database. A negative there is uninformative by
construction -- a peptide that exists only in 717 could never have been in that database.

This ships the re-search instead: the same class of raw spectra searched against a database
that contains the 717 proteins, across five public Populus datasets, 1% peptide AND spectrum
FDR. Now "not detected" is a measurement.

AND IT SHIPS THE DENOMINATOR. `detectable` says whether the gene has any fully-tryptic
peptide in the searchable window at all. Without it, "not detected" and "could not have been
detected" render identically, which is the same defect as a control that never ran.

WHAT IS LOST, stated rather than left implicit: Tier 1 covered PXD023826 and PXD064017,
which are not among the five re-searched here (one is DIA and needs a different engine, the
other a Waters vendor format). This layer is deeper and interpretable where Tier 1 was
neither, but it is not a strict superset of datasets. The provenance sidecar says so too.

The source table is gated at generation: make_proteomics_v2_source.py refuses to write
unless its per-gene roll-up reproduces the owning analysis's published per-class table in
all 84 cells. This script re-asserts the two totals it can check cheaply.

AFTER THIS RUNS, build_facets.py must run if the facet index ever carries a prot field.
It does not today, but TODO 24's rule is the default, not the exception.
"""
import csv
import gzip
import io
import json
import os
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
SRC = os.path.join(ROOT, "sources", "proteomics_ms_v2.tsv.gz")
PROV = os.path.join(ROOT, "sources", "proteomics_ms_v2.provenance.json")

# what the owning analysis published, per class, pooled over haplotypes; re-asserted here
PUBLISHED_UNION_POOL_DETECTED = 8   # LSG-pool genes with any hit across all five datasets

# The same union, per haplotype and stratum, as the owning analysis states it (ms2-lsg
# analyses/proteomics_translation_evidence/RUN_LOG.md, "The headline stands, and is now
# multi-dataset"): (detected, genes in the stratum), Idesia-corrected ranks.
PUBLISHED_UNION_POOL = {"hap1": {18: (2, 433), 19: (1, 343)},
                        "hap2": {18: (5, 497), 19: (0, 338)}}

# The deepest dataset. The per-class gradient in the manifest is read off it, per haplotype.
GRADIENT_DATASET = "PXD025636"


def rerun_rank(rec):
    """The gene's rank in the genEra rerun that added Idesia, which is the ranking the owning
    analysis's published tables use. Since C5 (Chen, 2026-09-21 evening, reconfirmed 2026-09-26)
    that rerun is the pool again and `ps.rank` IS its rank; the frozen call sits under `ps.frz`.
    From 2026-09-21 morning to 2026-09-26 the atlas had it the other way round (`ps.rr`), which is
    why this is a function rather than a field read."""
    ps = rec.get("ps") or {}
    assert "rr" not in ps, "ps.rr is the C4 arrangement; run patch_lsg_pool.py first"
    return ps.get("rank")


def gene_class(rank):
    """Copied from make_proteomics_v2_source.py, which copies the owning analysis."""
    if rank is None:
        return "unassigned"
    if rank <= 2:
        return "PS1_2_ancient"
    if rank <= 12:
        return "PS3_12_conserved"
    if rank <= 15:
        return "PS13_15_rosids"
    if rank <= 17:
        return "PS16_17_salicaceae"
    return "LSG_PS18_populus" if rank == 18 else "LSG_PS19_hybrid"


def age_gradient(data_dir):
    """Per-haplotype detection rate by age class in GRADIENT_DATASET, from the shards as shipped,
    plus the five-dataset pool union per stratum.

    Replaces a hand-typed, haplotype-pooled dict ({"LSG_pool": 0.12, ...}) that quoted the
    PUBLISHED Tier 2 figures, including a pool hit the owning analysis withdrew on 2026-09-09
    after the expanded-database re-search this layer actually ships. Computing it from the
    shards means the manifest cannot disagree with the genes it summarises. Asserts the pool
    union against the owner's published numbers.
    """
    n, det, union = Counter(), Counter(), Counter()
    for hap in ("hap1", "hap2"):
        for p in shard_paths(os.path.join(data_dir, "genes", hap)):
            for rec in json.load(open(p)).values():
                rank = rerun_rank(rec)
                c = gene_class(rank)
                ms = (rec.get("prot") or {}).get("ms") or {}
                n[(hap, c)] += 1
                det[(hap, c)] += GRADIENT_DATASET in (ms.get("ds") or [])
                if rank is not None and rank >= 18:
                    union[(hap, rank, "n")] += 1
                    union[(hap, rank, "det")] += (ms.get("n") or 0) > 0
    out = {"dataset": GRADIENT_DATASET,
           "note": "percent of genes with a peptide in this dataset, per haplotype, never pooled; strata are those of the genEra rerun that added Idesia, as the owning analysis reports them",
           "pool_union_all_datasets": {}}
    for hap in ("hap1", "hap2"):
        out[hap] = {c: round(100 * det[(hap, c)] / n[(hap, c)], 3)
                    for c in sorted({k[1] for k in n if k[0] == hap}) if n[(hap, c)]}
        out["pool_union_all_datasets"][hap] = {}
        for rank, (want_det, want_n) in PUBLISHED_UNION_POOL[hap].items():
            got = (union[(hap, rank, "det")], union[(hap, rank, "n")])
            assert got == (want_det, want_n), (
                f"{hap} PS{rank}: shards give {got[0]} of {got[1]} detected across the five "
                f"datasets; the owning analysis published {want_det} of {want_n}")
            out["pool_union_all_datasets"][hap][f"PS{rank}"] = {"detected": got[0], "genes": got[1]}
    return out


def iv(v):
    v = (v or "").strip()
    return None if v == "" else int(float(v))


def shard_paths(root):
    for dp, _, fs in os.walk(root):
        for fn in sorted(fs):
            if fn.endswith(".json"):
                yield os.path.join(dp, fn)


def main():
    with gzip.open(SRC, "rt") as fh:
        rows = list(csv.DictReader(io.StringIO(fh.read()), delimiter="\t"))
    prov = json.load(open(PROV))
    n_ds = len(prov["datasets"])
    print(f"source: {len(rows)} rows, {n_ds} datasets")

    ms = {}
    for r in rows:
        n = iv(r["n_datasets_detected"]) or 0
        rec = {"n": n, "of": n_ds}
        if n:
            rec["pep"] = iv(r["n_peptides_total"])
            rec["psm"] = iv(r["n_psms_total"])
            if r["best_hyperscore"]:
                rec["hs"] = float(r["best_hyperscore"])
            rec["sp"] = iv(r["any_717_specific_peptide"])
            rec["ds"] = r["datasets_detected"].split(",")
        tw = iv(r["n_tryptic_in_window"])
        if tw is not None:
            rec["tw"] = tw
            rec["uq"] = iv(r["n_unique_tryptic"])
            rec["abl"] = iv(r["detectable"])
        ms[r["gene_id"]] = rec

    patched, dropped_v1, detected, undetectable = 0, 0, Counter(), Counter()
    pool_detected = 0
    for hap in ("hap1", "hap2"):
        for p in shard_paths(os.path.join(DATA, "genes", hap)):
            shard = json.load(open(p))
            for gid, rec in shard.items():
                prot = rec.setdefault("prot", {})
                if isinstance(prot.get("ms"), dict) and any(
                        k.isdigit() for k in prot["ms"]):
                    dropped_v1 += 1     # the Tier-1 {"023826": [n, det]} shape
                m = ms.get(gid)
                if m is None:
                    prot.pop("ms", None)
                    continue
                prot["ms"] = m
                patched += 1
                if m["n"]:
                    detected[hap] += 1
                    if (rerun_rank(rec) or -1) >= 18:   # the owner counts rerun strata
                        pool_detected += 1
                if m.get("abl") == 0:
                    undetectable[hap] += 1
            with open(p, "w") as fh:
                json.dump(shard, fh, separators=(",", ":"), sort_keys=True)

    print(f"patched {patched} genes; removed {dropped_v1} Tier-1 records")
    print(f"  detected in >= 1 dataset: {dict(detected)} (total {sum(detected.values())})")
    print(f"  undetectable by construction: {dict(undetectable)}")
    print(f"  LSG-pool genes with any hit: {pool_detected}")
    assert patched == len(ms), f"{len(ms)} source rows but {patched} patched"
    assert pool_detected == PUBLISHED_UNION_POOL_DETECTED, (
        f"the owning analysis reports {PUBLISHED_UNION_POOL_DETECTED} pool genes with a hit "
        f"across the union of datasets; this build says {pool_detected}. Either the pool "
        "moved under this layer or the roll-up drifted -- stop and find out which.")

    man_path = os.path.join(DATA, "meta", "manifest.json")
    man = json.load(open(man_path))
    man["proteomics"] = {
        "version": 2,
        # 135,843, not the original Tier 2's 98,659: all five datasets were searched under the
        # expanded database (Tier 2 set + P. tremula Potra02), per the owning RUN_LOG.
        "method": "raw-spectra re-search against a database CONTAINING the 717 proteins "
                  "(135,843 targets: P. trichocarpa, HAP1, HAP2 and P. tremula Potra02), "
                  "Sage 0.14.7, 1% peptide AND spectrum FDR",
        "datasets": {k: v["species"] for k, v in prov["datasets"].items()},
        "n_datasets": n_ds,
        "genes_with_a_hit": dict(detected),
        "undetectable_by_construction": dict(undetectable),
        "lsg_pool_genes_with_a_hit": pool_detected,
        "supersedes": "Tier-1 string matching against three third-party peptide tables built "
                      "from another species' database (PXD023826 / PXD025636 / PXD064017). A "
                      "Tier-1 negative was uninformative by construction. PXD023826 and "
                      "PXD064017 are not covered by this layer -- the replacement is deeper "
                      "and interpretable, not a strict superset of datasets.",
        "engine_caveat": "A Sage 1% FDR is not an MSFragger 1% FDR. Report the engine with "
                         "the number.",
        "union_fdr_caveat": "1% per dataset ACCUMULATES over five datasets. Read a "
                            "single-dataset singleton against the pooled false-peptide "
                            "expectation, not one dataset's.",
        "cross_species_caveat": "Three of the five datasets are other Populus species, whose "
                                "samples contain no 717 protein. A hit there says the gene "
                                "FAMILY is translated somewhere in Populus, not that 717's "
                                "copy is.",
        "detectability_note": "`abl` says whether the gene has any fully-tryptic peptide in "
                              "the searchable window. Read it before reading a zero: a "
                              "protein with none is invisible however abundant it is, and "
                              "that is a composition property on which young genes differ.",
        "age_gradient_pct": age_gradient(DATA),
    }
    with open(man_path, "w") as fh:
        json.dump(man, fh, indent=1, sort_keys=True)
    print("manifest: proteomics block written")


if __name__ == "__main__":
    main()
