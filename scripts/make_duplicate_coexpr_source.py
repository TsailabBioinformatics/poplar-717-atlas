#!/usr/bin/env python3
"""Expression correlation between a gene and each of its duplication partners.

Chen asked the duplication card to name WHO shares the relationship and show their Ks and
expression divergence. Ks already ships per paralog; the co-expression does not exist anywhere
upstream -- everything found was NEIGHBOUR co-expression (neighbor_coexpr_652, shipped as `nb`)
or SYNTELOG/allele divergence (ase_deepening_prott5). So this computes it.

BASIS: Pearson r on log2(ComBat-corrected TPM + 1) across all 652 samples.

Two reasons for ComBat rather than the raw TPM the atlas ships as primary, and the second is
the real one:
  1. it matches the shipped `nb` field, which s06/s36 compute as mean PCC on
     combat_log2_652.parquet -- otherwise one gene page would carry two co-expression numbers
     on different matrices;
  2. a correlation across 652 samples from 28 studies IS a cross-study comparison, and the
     atlas's own About page says raw TPM is not suitable for that ("cross-study comparison
     needs batch handling that this build does not apply").
The cost is coverage: 12,440 of 13,995 duplicate-partner genes have a ComBat vector, so pairs
where either partner lacks one get NO value rather than a value on a different basis.

This is a DESCRIPTIVE correlation, deliberately not called a divergence measure. No model, no
threshold, no claim about what a given r implies for duplicate fate.
"""
import csv, gzip, io, json, os
import numpy as np

from expr_range_reader import load_combat_vectors

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "sources", "duplicate_coexpr.tsv.gz")


def shard_paths(root):
    for dp, _, fs in os.walk(root):
        for fn in sorted(fs):
            if fn.endswith(".json"):
                yield os.path.join(dp, fn)


def main():
    # gene -> combat bucket pointer, and the pair list
    ptr, pairs = {}, set()
    for hap in ("hap1", "hap2"):
        for p in shard_paths(os.path.join(DATA, "genes", hap)):
            for gid, r in json.load(open(p)).items():
                x = r.get("x") or {}
                if x.get("cb"):
                    ptr[gid] = {"hap": hap, "range": x.get("co"), "bucket": x["cb"]}
                for q in (r.get("para") or []):
                    pairs.add(tuple(sorted((gid, q["id"]))))
    arrays = json.load(open(os.path.join(DATA, "index", "arrays.json")))
    n_para = len(pairs)
    for members in arrays.values():
        m = sorted(members)
        for i in range(len(m)):
            for j in range(i + 1, len(m)):
                pairs.add((m[i], m[j]))
    pairs = sorted(pairs)
    print(f"{len(pairs):,} distinct duplicate pairs ({n_para:,} paralog, "
          f"{len(pairs)-n_para:,} tandem co-member)")
    print(f"{len(ptr):,} genes have a ComBat vector")

    need = {g for pr in pairs for g in pr if g in ptr}
    vectors, storage = load_combat_vectors(DATA, {g: ptr[g] for g in need})
    vec = {g: np.log2(np.asarray(values, dtype=np.float64) + 1.0)
           for g, values in vectors.items()}
    print(f"loaded {len(vec):,} vectors ({storage['range']:,} packed ranges; "
          f"{storage['legacy']:,} from {storage['legacy_buckets']:,} legacy buckets)")

    # centre once; Pearson is then a normalised dot product
    cent = {}
    for g, v in vec.items():
        d = v - v.mean()
        n = float(np.sqrt((d * d).sum()))
        cent[g] = (d, n)

    rows, skipped_missing, skipped_flat = [], 0, 0
    for a, b in pairs:
        if a not in cent or b not in cent:
            skipped_missing += 1
            continue
        da, na = cent[a]
        db, nb_ = cent[b]
        if na == 0 or nb_ == 0:          # a gene flat across all 652 samples has no correlation
            skipped_flat += 1
            continue
        rows.append({"gene_a": a, "gene_b": b,
                     "r": round(float((da * db).sum() / (na * nb_)), 4)})
    print(f"computed {len(rows):,} correlations; {skipped_missing:,} pairs lack a ComBat vector "
          f"on one side, {skipped_flat:,} have a flat partner")

    # GATES ------------------------------------------------------------------------------
    # 1. the implementation itself, against numpy on a sample
    rng = np.random.RandomState(0)
    idx = rng.choice(len(rows), min(300, len(rows)), replace=False)
    worst = 0.0
    for i in idx:
        a, b, r = rows[i]["gene_a"], rows[i]["gene_b"], rows[i]["r"]
        worst = max(worst, abs(np.corrcoef(vec[a], vec[b])[0, 1] - r))
    assert worst < 1e-3, f"Pearson implementation disagrees with numpy by {worst}"
    print(f"GATE PASS: matches numpy.corrcoef on 300 sampled pairs (max delta {worst:.2e})")
    # 2. a gene against itself must be exactly 1
    g0 = rows[0]["gene_a"]
    d0, n0 = cent[g0]
    assert abs(float((d0 * d0).sum() / (n0 * n0)) - 1.0) < 1e-9, "self-correlation is not 1"
    print("GATE PASS: self-correlation is 1")
    # 3. r must lie in [-1, 1]
    bad = [r for r in rows if not (-1.0001 <= r["r"] <= 1.0001)]
    assert not bad, f"{len(bad)} correlations outside [-1,1]"
    print("GATE PASS: every r inside [-1, 1]")

    rr = np.array([r["r"] for r in rows])
    print(f"  distribution: median {np.median(rr):+.3f}  "
          f"p10 {np.percentile(rr,10):+.3f}  p90 {np.percentile(rr,90):+.3f}")

    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=["gene_a", "gene_b", "r"], delimiter="\t",
                       lineterminator="\n")
    w.writeheader(); w.writerows(rows)
    with open(OUT, "wb") as fh:
        with gzip.GzipFile(fileobj=fh, mode="wb", mtime=0) as gz:
            gz.write(buf.getvalue().encode())
    print(f"wrote {OUT}: {len(rows):,} pairs, {os.path.getsize(OUT)/1e6:.2f} MB gz")


if __name__ == "__main__":
    main()
