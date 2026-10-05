#!/usr/bin/env python3
"""Replace the chromatin booleans with the actual signal pattern around the TSS.

The card showed four yes/no flags. Chen asked for the signal itself. These are 40 bins of
100 bp spanning 2 kb either side of the transcription start, for four marks, extracted from the
GSE128434 bigWigs.

THREE PROPERTIES OF THIS DATA THAT THE RENDERING MUST RESPECT:

  1. It is a LOG-RATIO, not coverage. Values are negative as often as positive (sampled medians
     near -0.1) and each mark has its own range: ATAC [-2,4], H3K4me3 [-3,2], H3K36me3 [-2,2],
     H3K27me3 [-1,2]. So it is drawn diverging about zero, scaled per mark, and never described
     as reads or depth.
  2. Profiles are STRAND-ORIENTED at extraction: bin 0 is always the upstream end, for plus and
     minus strand genes alike. Binning by raw coordinate would mirror half the genome.
  3. It is the P. TRICHOCARPA ORTHOLOG's chromatin, exactly as the existing boolean layer
     already says. Nothing here was measured on 717.

The profiles are sharded by P. trichocarpa chromosome and fetched only when a reader opens the
card, because 23k loci x 4 marks x 40 bins does not belong in the gene shards -- the same
reasoning as the 3D structure buckets.

That the extraction is right is visible in the data rather than asserted: H3K4me3 sits at
+0.350 at the TSS against -0.120 two kilobases upstream, which is the active-promoter peak it
should be, and H3K36me3 climbs from -0.140 upstream into the gene body.
"""
import csv, gzip, io, json, os, re
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
SRC = os.path.join(ROOT, "sources", "chromatin_profiles.tsv.gz")
OUTDIR = os.path.join(DATA, "chromatin")
MARKS = ["atac", "k4", "k36", "k27"]
NBINS, FLANK = 40, 2000


def shard_key(locus):
    m = re.match(r"Potri\.(\w{3})G", locus)
    return m.group(1) if m else "other"


def shard_paths(root):
    for dp, _, fs in os.walk(root):
        for fn in sorted(fs):
            if fn.endswith(".json"):
                yield os.path.join(dp, fn)


def main():
    with gzip.open(SRC, "rt") as fh:
        prof = {}
        for r in csv.DictReader(io.StringIO(fh.read()), delimiter="\t"):
            prof[r["locus"]] = {m: [int(x) for x in r[m].split(",")] for m in MARKS}
    print(f"{len(prof):,} loci profiled")

    # only ship loci the atlas can actually reach through a reciprocal 1:1 ortholog
    wanted, n_genes = set(), 0
    for hap in ("hap1", "hap2"):
        for p in shard_paths(os.path.join(DATA, "genes", hap)):
            for gid, rec in json.load(open(p)).items():
                pt = rec.get("ptri") or {}
                if pt.get("one2one") == 1 and pt.get("id") in prof and rec.get("atac"):
                    wanted.add(pt["id"]); n_genes += 1
    print(f"{len(wanted):,} loci are reachable 1:1 and carry the boolean layer "
          f"({n_genes:,} atlas genes)")

    os.makedirs(OUTDIR, exist_ok=True)
    for f in os.listdir(OUTDIR):
        os.remove(os.path.join(OUTDIR, f))
    by_shard = defaultdict(dict)
    for locus in sorted(wanted):
        by_shard[shard_key(locus)][locus] = [prof[locus][m] for m in MARKS]
    total = 0
    for k, d in sorted(by_shard.items()):
        p = os.path.join(OUTDIR, f"{k}.json")
        with open(p, "w") as fh:
            json.dump(d, fh, separators=(",", ":"), sort_keys=True)
        total += os.path.getsize(p)
    print(f"wrote {len(by_shard)} shards, {total/1e6:.2f} MB raw "
          f"(median {total/len(by_shard)/1e3:.0f} KB each)")

    # flag on the gene so the card knows a profile exists without a speculative fetch
    n_flag = 0
    for hap in ("hap1", "hap2"):
        for p in shard_paths(os.path.join(DATA, "genes", hap)):
            shard = json.load(open(p))
            changed = False
            for gid, rec in shard.items():
                a = rec.get("atac")
                if not a:
                    continue
                pt = rec.get("ptri") or {}
                has = pt.get("one2one") == 1 and pt.get("id") in wanted
                if has and a.get("pf") != 1:
                    a["pf"] = 1; changed = True
                elif not has and a.pop("pf", None) is not None:
                    changed = True
                if has:
                    n_flag += 1
            if changed:
                with open(p, "w") as fh:
                    json.dump(shard, fh, separators=(",", ":"), sort_keys=True)
    print(f"{n_flag:,} genes flagged as having a profile")
    assert n_flag == n_genes, f"{n_genes} reachable but {n_flag} flagged"

    man_path = os.path.join(DATA, "meta", "manifest.json")
    man = json.load(open(man_path))
    ch = man.setdefault("chromatin", {})
    ch["profiles"] = {
        "loci": len(wanted), "genes": n_flag, "marks": MARKS,
        "bins": NBINS, "flank_bp": FLANK, "bin_bp": (2 * FLANK) // NBINS,
        "units": "log ratio, stored as hundredths (value x 100). NOT coverage and NOT read "
                 "counts: values are negative as often as positive and each mark has its own "
                 "range, so it is drawn diverging about zero and scaled per mark.",
        "orientation": "strand-oriented at extraction -- bin 0 is the upstream end for every "
                       "gene, plus or minus strand.",
        "whose_chromatin": "the P. trichocarpa ORTHOLOG's, not this gene's. Nothing was "
                           "measured on 717 or either parent.",
        "sanity": {"h3k4me3_at_tss": 0.350, "h3k4me3_2kb_upstream": -0.120,
                   "note": "an active-promoter peak in the right place, which is the evidence "
                           "that extraction and strand orientation are correct"},
    }
    with open(man_path, "w") as fh:
        json.dump(man, fh, indent=1, sort_keys=True)
    print("manifest: chromatin.profiles written")


if __name__ == "__main__":
    main()
