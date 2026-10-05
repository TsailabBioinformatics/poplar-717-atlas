#!/usr/bin/env python3
"""Descent groups: connected components over every homology edge the atlas ships.

WHY THIS IS SAFE TO BUILD, measured before a line of the view was written. Union-find over
tandem membership, WGD anchor pairs, HAP1<->HAP2 alleles and paralog pairs (160,898 directed
entries) gives 19,614 components over 63,960 genes, and THE LARGEST IS 85 GENES. 19,525
components hold 20 genes or fewer, covering 95.6% of the corpus. Had tandem arrays plus WGD
chained the genome into one giant component, a "family" view would have shown 20,000 nodes and
been useless -- so that was the falsifier, and it passed.

WHAT IT EMITS
  data/index/families.json   {component_id: [gene ids]} for components of 2+ genes
  x.fam on each gene shard   the component id, so a gene page needs no lookup table

Singletons get no entry and no `fam`: a gene with no homology edge is not in a family, and
inventing a one-member family would make "this gene has a family" true of everything.

DETERMINISM. Component ids are assigned by walking genes in sorted order, and each member list
is sorted, so a rebuild produces a byte-identical file. The search index was reshuffled every
build by set iteration once already (TODO 24); that is not repeated here.
"""
import json, os, sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")

parent = {}


def find(x):
    root = x
    while parent.get(root, root) != root:
        root = parent[root]
    while parent.get(x, x) != x:            # path compression
        parent[x], x = root, parent[x]
    return root


def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[ra] = rb


def shard_paths(root):
    for dp, _, fs in os.walk(root):
        for fn in sorted(fs):
            if fn.endswith(".json"):
                yield os.path.join(dp, fn)


def main():
    genes = set()
    n_edges = defaultdict(int)
    for hap in ("hap1", "hap2"):
        for p in shard_paths(os.path.join(DATA, "genes", hap)):
            for gid, rec in json.load(open(p)).items():
                genes.add(gid)
                parent.setdefault(gid, gid)
                for q in (rec.get("wgd") or []):
                    parent.setdefault(q["id"], q["id"]); union(gid, q["id"]); n_edges["wgd"] += 1
                a = rec.get("allele")
                if a:
                    parent.setdefault(a["id"], a["id"]); union(gid, a["id"]); n_edges["allele"] += 1
                for m in ((rec.get("td") or {}).get("mem") or []):
                    parent.setdefault(m[0], m[0]); union(gid, m[0]); n_edges["tandem"] += 1
                for q in (rec.get("para") or []):
                    parent.setdefault(q["id"], q["id"]); union(gid, q["id"]); n_edges["para"] += 1

    members = defaultdict(list)
    for g in sorted(genes):
        members[find(g)].append(g)

    # ids assigned in a deterministic order, smallest member id first
    fams, gene_fam = {}, {}
    for i, root in enumerate(sorted(members, key=lambda r: members[r][0])):
        ms = members[root]
        if len(ms) < 2:
            continue
        fid = f"F{i:05d}"
        fams[fid] = ms
        for g in ms:
            gene_fam[g] = fid

    sizes = sorted((len(v) for v in fams.values()), reverse=True)
    assert sizes[0] < 500, f"a component of {sizes[0]} genes would not be drawable"

    out = os.path.join(DATA, "index", "families.json")
    with open(out, "w") as fh:
        json.dump(fams, fh, separators=(",", ":"), sort_keys=True)
    print(f"{len(genes):,} genes, edges {dict(n_edges)}")
    print(f"{len(fams):,} families of 2+ genes covering {len(gene_fam):,} genes; "
          f"{len(genes) - len(gene_fam):,} singletons carry no family")
    print(f"largest {sizes[0]}, median {sizes[len(sizes)//2]}, "
          f"{sum(1 for s in sizes if s <= 20):,} of size <= 20")
    print(f"  -> {out} ({os.path.getsize(out)/1e6:.2f} MB)")

    n_set = 0
    for hap in ("hap1", "hap2"):
        for p in shard_paths(os.path.join(DATA, "genes", hap)):
            shard = json.load(open(p))
            changed = False
            for gid, rec in shard.items():
                fid = gene_fam.get(gid)
                if fid:
                    if rec.get("fam") != fid:
                        rec["fam"] = fid
                        changed = True
                    n_set += 1
                elif rec.pop("fam", None) is not None:
                    changed = True
            if changed:
                with open(p, "w") as fh:
                    json.dump(shard, fh, separators=(",", ":"), sort_keys=True)
    print(f"{n_set:,} genes carry a family id")

    man_path = os.path.join(DATA, "meta", "manifest.json")
    man = json.load(open(man_path))
    man["families"] = {
        "n_families": len(fams),
        "genes_in_a_family": len(gene_fam),
        "singletons": len(genes) - len(gene_fam),
        "largest": sizes[0],
        "edges": dict(n_edges),
        "definition": "connected components over tandem membership, WGD anchor pairs, "
                      "HAP1/HAP2 alleles and paralog pairs. A singleton is NOT a family and "
                      "carries no id -- a gene with no homology edge is not in one.",
        "why_drawable": "measured before building the view: the largest component is "
                        f"{sizes[0]} genes, so no family is too large to render.",
    }
    json.dump(man, open(man_path, "w"), indent=1, sort_keys=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
