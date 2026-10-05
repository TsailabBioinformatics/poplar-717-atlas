#!/usr/bin/env python3
"""Ship sources/copy_search.json as data/copies/<locus>.json, plus a `copies` manifest section.

The search itself is make_copy_search_source.py (Sapelo2, run by hand). This stage only checks
the committed source against the gene shards and publishes it, so it runs on any machine.

  data/copies/<slug>.json   the locus as the search wrote it: regions, reference, rows, pieces
  manifest.copies.genes     gene id -> locus slug, for every ANNOTATED gene drawn in a panel,
                            so a gene page knows synchronously whether it has one
  manifest.copies.loci      per locus: name, and the counts the gene page's summary strip
                            states (so the strip computes nothing)

Gene shards are NOT touched. Refuses to write if the source disagrees with them: a panel row for
an annotated gene must carry exactly that gene's coordinates, or the panel and the gene page
would show two different positions for one gene.
"""
import hashlib
import json
import os
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
SRC = os.path.join(ROOT, "sources", "copy_search.json")
OUT = os.path.join(DATA, "copies")


def shard_record(gid):
    """(hap, record) from the gene's shard; record is None if absent. Imports the shard rule
    from build_data.py rather than restating it."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from build_data import shard_of
    hap, key = shard_of(gid)
    path = os.path.join(DATA, "genes", hap, f"{key}.json")
    return hap, (json.load(open(path)).get(gid) if os.path.exists(path) else None)


def main():
    raw = open(SRC, "rb").read()
    src = json.loads(raw)
    fails, genes, loci = [], {}, {}

    for slug, L in sorted(src["loci"].items()):
        n_got, n_arr = L["controls"]["array_members_recovered"]
        if n_got != n_arr:
            fails.append(f"{slug}: control failed at source ({n_got}/{n_arr} members recovered)")
        drawn = {r["gene"] for r in L["rows"] if r["gene"]}
        if not set(L["array"]) <= drawn:
            fails.append(f"{slug}: array members missing from the panel: {sorted(set(L['array']) - drawn)}")
        for r in L["rows"]:
            if not r["gene"]:
                continue
            hap, rec = shard_record(r["gene"])
            if rec is None:
                fails.append(f"{slug}: {r['gene']} is not in the gene shards")
                continue
            if hap != r["hap"] or (rec["chr"], rec["start"], rec["end"]) != (r["chr"], r["start"], r["end"]):
                fails.append(f"{slug}: {r['gene']} panel says {r['hap']} {r['chr']}:{r['start']}-{r['end']}, "
                             f"shard says {hap} {rec['chr']}:{rec['start']}-{rec['end']}")
            genes[r["gene"]] = slug
        for h, reg in L["regions"].items():
            for f in reg["flanks"]:
                if shard_record(f)[1] is None:
                    fails.append(f"{slug}: region flank {f} is not in the gene shards")

        summary = {}
        for h in sorted({r["hap"] for r in L["rows"]}):
            rows = [r for r in L["rows"] if r["hap"] == h]
            summary[h] = {
                "annotated": sum(1 for r in rows if r["gene"]),
                "unannotated": sum(1 for r in rows if not r["gene"]),
                "pieces": sum(len(r["pieces"]) for r in rows),
                "max_cov_pct": max(r["cov_pct"] for r in rows),
            }
        loci[slug] = {"name": L["name"], "desc": L["desc"], "hap": L["hap"],
                      "reference": L["reference"]["id"], "summary": summary}

    if fails:
        print("REFUSING TO WRITE -- the source disagrees with the gene shards:")
        for f in fails:
            print("  " + f)
        sys.exit(1)

    if os.path.isdir(OUT):
        shutil.rmtree(OUT)          # data/copies is owned by this stage alone
    os.makedirs(OUT)
    for slug, L in src["loci"].items():
        with open(os.path.join(OUT, f"{slug}.json"), "w") as fh:
            json.dump({"slug": slug, **L}, fh, separators=(",", ":"), sort_keys=True)

    man_path = os.path.join(DATA, "meta", "manifest.json")
    man = json.load(open(man_path))
    man["copies"] = {
        "genes": dict(sorted(genes.items())),
        "loci": loci,
        "source": "sources/copy_search.json",
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "method": src["method"],
        "note": "Copies found in the genome SEQUENCE, so a panel can show pieces the annotation "
                "does not call. Only the loci listed here have been searched; a gene without a "
                "panel has not been searched, which is not the same as having no copies.",
    }
    json.dump(man, open(man_path, "w"), indent=1)
    print(f"copies: {len(loci)} locus, {len(genes)} annotated genes carry a panel -> data/copies/")
    for slug, l in loci.items():
        print(f"  {slug}: " + "; ".join(f"{h} {s['annotated']} annotated + {s['unannotated']} not, "
                                         f"{s['pieces']} pieces, longest {s['max_cov_pct']}%"
                                         for h, s in l["summary"].items()))


if __name__ == "__main__":
    main()
