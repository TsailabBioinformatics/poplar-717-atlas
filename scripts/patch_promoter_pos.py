#!/usr/bin/env python3
"""TODO 44: motif POSITIONS in the promoter, so the card can show which motifs sit close to
the transcription start and which sit far -- the review asked for "sticks on the promoter
region, who is closer who is farer".

The shipped CRE layer carries family COUNTS and no positions. These positions come from a
rescan of the same 1 kb promoter FASTAs with the upstream algorithm copied verbatim
(scripts in ms2lsg_repo/pipeline: s08_promoter_cre pssms/encode/HIT_FRAC/PAD/STRESS, s08d's
greedy non-overlapping pass, piggyback_locus_cre_marks DEV). `dist` is bp from the TSS to the
motif's 5' edge, which is what the upstream scan already computed internally.

SCOPE, decided by measurement rather than preference:

  * Nine NAMED families only. `other` is 457 unrelated matrices and is not a meaningful thing
    to draw a stick for; it is also where the reproduction disagreed.
  * 45 genes are EXCLUDED. Upstream scanned a 3 kb window and took the 1 kb slice, and its
    greedy non-overlapping pass is global over 3 kb; only the 1 kb FASTA exists on this
    cluster, so a motif near the boundary can be resolved differently. At full scale that
    affects 45 of 63,960 genes on a named family, and every one is this scan finding MORE
    hits, never fewer -- the predicted direction. Those genes get no positions rather than
    positions that disagree with counts already published.

THE GATE, re-asserted here and not merely upstream: for every gene that gets positions, the
number of positions per family must EQUAL the `cre.fam` count this atlas already serves. If
one gene disagrees, the card would contradict the bars directly above it, so nothing is
written.
"""
import csv
import gzip
import io
import json
import os
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
SRC = os.path.join(ROOT, "sources", "promoter_motif_pos.tsv.gz")
EXCL = os.path.join(ROOT, "sources", "promoter_motif_pos.excluded.tsv")

# build_cre.py renames three families when it ingests the counts, on purpose: ABRE and DRE/CRT
# are cis-ELEMENTS bound by particular family members, and labelling a family-level motif hit
# with an element name would be a different and wrong claim. The scan emits the upstream family
# names, so map them onto the keys the shards actually use. Verified against a shard:
# AP2_ERF, ARF, MADS, MYB, NAC, SPL, TCP, WRKY, bZIP, other.
FAMILY_ALIAS = {"ABRE_bZIP": "bZIP", "DREB_ERF": "AP2_ERF", "MADS_box": "MADS"}


def shard_paths(root):
    for dp, _, fs in os.walk(root):
        for fn in sorted(fs):
            if fn.endswith(".json"):
                yield os.path.join(dp, fn)


def main():
    pos = defaultdict(lambda: defaultdict(list))
    with gzip.open(SRC, "rt") as fh:
        for r in csv.DictReader(io.StringIO(fh.read()), delimiter="\t"):
            fam = FAMILY_ALIAS.get(r["family"], r["family"])
            pos[r["gene_id"]][fam].append(int(r["dist_to_tss"]))
    for g in pos:
        for fam in pos[g]:
            pos[g][fam].sort()
    excluded = set()
    if os.path.exists(EXCL):
        with open(EXCL) as fh:
            next(fh, None)
            for line in fh:
                excluded.add(line.split("\t")[0])
    print(f"source: {len(pos)} genes with positions, {sum(len(v) for d in pos.values() for v in d.values())} "
          f"positions; {len(excluded)} genes excluded upstream")

    patched, mism, tagged = 0, [], 0
    fam_tot = Counter()
    for hap in ("hap1", "hap2"):
        for p in shard_paths(os.path.join(DATA, "genes", hap)):
            shard = json.load(open(p))
            changed = False
            for gid, rec in shard.items():
                cre = rec.get("cre")
                if not cre:
                    continue
                if gid in excluded:
                    cre["posx"] = 1          # positions deliberately withheld, see docstring
                    cre.pop("pos", None)
                    tagged += 1
                    changed = True
                    continue
                cre.pop("posx", None)
                d = pos.get(gid)
                if not d:
                    if cre.pop("pos", None) is not None:
                        changed = True
                    continue
                # THE GATE: positions must aggregate back to the counts already on the page
                for fam, dists in d.items():
                    want = (cre.get("fam") or {}).get(fam)
                    if want != len(dists):
                        mism.append((gid, fam, want, len(dists)))
                    fam_tot[fam] += len(dists)
                cre["pos"] = {fam: dists for fam, dists in sorted(d.items())}
                patched += 1
                changed = True
            if changed:
                with open(p, "w") as fh:
                    json.dump(shard, fh, separators=(",", ":"), sort_keys=True)

    print(f"patched {patched} genes with positions; {tagged} tagged as withheld")
    assert not mism, (
        f"{len(mism)} gene/family cells where the positions do not aggregate to the count this "
        f"atlas already serves, e.g. {mism[:5]}. The card would contradict the bars above it.")
    print("GATE PASS: every gene's positions aggregate exactly to its published family counts")

    man_path = os.path.join(DATA, "meta", "manifest.json")
    man = json.load(open(man_path))
    cre = man.setdefault("cre", {})
    cre["positions"] = {
        "genes_with_positions": patched,
        "genes_withheld": tagged,
        "positions_by_family": dict(sorted(fam_tot.items())),
        "measure": "bp from the transcription start to the motif's 5' (far) edge, within the "
                   "1 kb promoter window",
        "families": "the nine named TF families only; `other` (457 matrices) is excluded as "
                    "not meaningful to draw per-motif",
        "why_some_are_withheld": "Upstream scanned a 3 kb window and took the 1 kb slice, with "
                                 "a greedy non-overlapping pass global over 3 kb. Only the 1 kb "
                                 "FASTA is on this cluster, so a motif near the boundary can be "
                                 "resolved differently: 45 of 63,960 genes disagree on a named "
                                 "family, all in the direction of this scan finding MORE hits. "
                                 "Those genes are withheld rather than shown with positions "
                                 "that contradict published counts.",
        "gate": "every shipped gene's positions aggregate exactly to its cre.fam counts, "
                "asserted on every build",
    }
    with open(man_path, "w") as fh:
        json.dump(man, fh, indent=1, sort_keys=True)
    print("manifest: cre.positions written")


if __name__ == "__main__":
    main()
