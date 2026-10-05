#!/usr/bin/env python3
"""
Checks that every source the data build reads is present, BEFORE anything is written.

Why this exists: all build inputs live on the UGA Sapelo2 cluster (/scratch, /work). On any
other machine the build cannot run at all. Without this check the scripts would half-succeed --
build_data.py would write shards and the patchers would silently add nothing -- leaving a
plausible-looking data/ directory missing eleven layers. A loud stop is much safer.

Exit 0 = every source present, safe to build.  Exit 1 = do not build here.
"""
import os, sys

SOURCES = {
    "gene models (BED)": [
        "<scratch>/ms2_genespace/workdir/bed/Ptrxalhap1.bed",
        "<scratch>/ms2_genespace/workdir/bed/Ptrxalhap2.bed"],
    "synteny (step 7 diff)": [
        "<scratch>/ms2_genespace/step7_diff_20260901/diff_all_genes.tsv"],
    "synteny (raw pangenes)": [
        "<scratch>/ms2_genespace/tier3_input_v2_20260901_runA/pangenes_long_RAW.tsv.gz"],
    "phylostratigraphy (genEra)": [
        "<scratch>/genera/output/hap1/genera_717hap1_v5.1/80863_gene_ages.tsv",
        "<scratch>/genera/output/hap2/genera_717hap2_v5.1/80863_gene_ages.tsv"],
    "expression (652-sample v3)": [
        "<scratch>/ms2_genespace/expression/COLUMN_LEVEL_sample_metadata.csv",
        "<scratch>/ms2_genespace/expression/combined_tpm_k21_rounded_20260814_study_sorted.csv"],
    "functional annotation (Phytozome v5.1)": [
        "<scratch>/genome_resources/pangenome/input/717Official/annotation/"
        "PtremulaxPopulusalbaHAP1_717_v5.1.annotation_info.txt"],
    "allele Ka/Ks": ["<scratch>/synt_ks/syntelog_ks.csv"],
    "duplication class (ms1-dup)": [
        "<scratch>/ms1-dup/ms1-figures/data/dupclass_membership_v8.csv"],
    "tandem arrays (curated union)": [
        "<scratch>/ms1-dup/data/09_td_pairs/td_union_20260831/union_genes_final.csv"],
    "perturbation responsiveness": [
        "<scratch>/ms1-dup/data/54_v8_freeze_20260825/perturbation_per_gene_v8.csv"],
    "protein layers (ms2-lsg)": [
        "<scratch>/ms2-lsg/analyses/strict_control_audit/results/predicted_disorder_genomewide.tsv",
        "<scratch>/ms2-lsg/analyses/esmfold_secondary_structure/results/secondary_structure_summary.tsv",
        "<scratch>/ms2-lsg/analyses/proteomics_translation_evidence/results/tier1_pxd064017_by_gene.tsv"],
    "reference annotation (/work)": [
        "/work/cjtlab/717Official/annotation/PtremulaxPopulusalbaHAP1_717_v5.1.gene_exons.gff3",
        "/work/cjtlab/717Official/annotation/PtremulaxPopulusalbaHAP1_717_v5.1.repeatmasked_assembly_v5.0.gff3"],
    "cross-species synteny (own GENESPACE run)": [
        "<scratch>/ms2_genespace/workdir_local_47798208/syntenicHits/"
        "Ptrichocarpa_vs_Ptrxalhap1.synHits.txt.gz"],
}

def main():
    missing = {}
    for label, paths in SOURCES.items():
        gone = [p for p in paths if not os.path.exists(p)]
        status = "OK  " if not gone else "MISS"
        print(f"  {status} {label}")
        if gone:
            missing[label] = gone
    if not missing:
        print("\nAll sources present — safe to run scripts/build_all.sh")
        return 0
    print(f"\n{len(missing)} SOURCE GROUP(S) MISSING. Do not build here.\n")
    for label, paths in missing.items():
        print(f"  {label}:")
        for p in paths:
            print(f"    {p}")
    print("\nThese live on the UGA Sapelo2 cluster. The data build only runs there.")
    print("The committed data/ directory IS the artifact — a clone can serve and edit the")
    print("site, run scripts/validate_data.py and scripts/check.mjs, but cannot rebuild data/.")
    return 1

if __name__ == "__main__":
    sys.exit(main())
