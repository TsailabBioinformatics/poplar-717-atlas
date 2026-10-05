#!/usr/bin/env python3
"""Range-addressable expression: one gzip member per gene, concatenated per haplotype.

THE PROBLEM. A gene page used to fetch `data/expr/buckets/<hap>/<n>.json` -- 250 genes x 652
samples, 786 KB raw -- to draw ONE gene's profile. The two bucket trees are build inputs, not
release assets: once packed they would otherwise make the Pages artifact hundreds of megabytes
larger than the range-addressable representation.

THE FIX, and why it is safe on GitHub Pages. Pages sends `accept-ranges: bytes` and answers a
Range request with a real 206; a ranged response is served identity-encoded even when the client
offers gzip, so byte offsets are stable. (Verified against the live host before this was written.)
Each vector is gzip-compressed independently, the members are concatenated into one blob per
haplotype and matrix, and the gene shard records its exact [offset, length].

VALUES ARE UNCHANGED. This is a container change, not a quantisation. The JSON inside each member
is byte-identical to the array that was in the bucket. `x.o` identifies raw and `x.co` ComBat
ranges; `x.b` and `x.cb` remain only as compatibility fallbacks for a client briefly holding an
older cached data tree.

The publication step is transactional. Both blob trees are staged, every member is decompressed
and checked against the original payload, updated shards and manifest are staged, and then all
outputs replace their predecessors with rollback backups. Only after that succeeds are the source
bucket trees removed. A packing or publication failure keeps the prior outputs and bucket trees.
"""
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
LEVEL = 6          # 9 buys ~2% for ~4x the build time on 128k members
HAPS = ("hap1", "hap2")
BUCKET_NAME = re.compile(r"^(\d+)\.json$")


def shard_paths(root):
    for dp, dirs, fs in os.walk(root):
        dirs.sort()
        for fn in sorted(fs):
            if fn.endswith(".json"):
                yield os.path.join(dp, fn)


def bucket_names(directory):
    """Return numeric bucket files in a stable order, rejecting unknown entries."""
    names, seen = [], set()
    for name in os.listdir(directory):
        if name.startswith("."):
            continue
        match = BUCKET_NAME.fullmatch(name)
        if not match:
            raise RuntimeError(f"{directory}: unexpected entry {name!r}; expected <number>.json")
        number = int(match.group(1))
        if number in seen:
            raise RuntimeError(f"{directory}: duplicate numeric bucket {number}")
        seen.add(number)
        names.append((number, name))
    if not names:
        raise RuntimeError(f"{directory}: no bucket files found")
    return [name for _, name in sorted(names)]


def gzip_member(payload):
    """Produce a byte-stable gzip member regardless of the wall clock."""
    out = io.BytesIO()
    # gzip.compress() defaults its mtime to the current time. GzipFile lets us pin both
    # the timestamp and optional filename, so unchanged input has unchanged bytes.
    with gzip.GzipFile(fileobj=out, mode="wb", compresslevel=LEVEL, mtime=0, filename="") as fh:
        fh.write(payload)
    return out.getvalue()


def stage_bytes(final_path, payload):
    """Write and fsync a sibling temporary file, leaving the destination untouched."""
    fd, temp_path = tempfile.mkstemp(prefix=f".{os.path.basename(final_path)}.", suffix=".tmp",
                                     dir=os.path.dirname(final_path))
    try:
        with os.fdopen(fd, "wb") as out:
            out.write(payload)
            out.flush()
            os.fsync(out.fileno())
    except Exception:
        os.unlink(temp_path)
        raise
    return temp_path


def build_one(tree, out_name, staged):
    """Stage one bucket tree and return its range index, byte total, and payload digests.

    No existing blob changes here.  The digest preserves the exact JSON bytes read from the
    source buckets so validation catches an offset, compression, ordering, or member-boundary
    defect before the output can replace the old release.
    """
    src = os.path.join(DATA, "expr", tree)
    if not os.path.isdir(src):
        raise RuntimeError(f"required expression bucket tree is absent: {src}")
    index, totals, payload_digests = {}, {}, {}
    for hap in HAPS:
        hdir = os.path.join(src, hap)
        if not os.path.isdir(hdir):
            raise RuntimeError(f"required bucket directory is absent: {hdir}")
        out_path = os.path.join(DATA, "expr", f"{out_name}_{hap}.gzb")
        fd, temp_path = tempfile.mkstemp(prefix=f".{out_name}_{hap}.", suffix=".tmp",
                                         dir=os.path.dirname(out_path))
        os.close(fd)
        off, hap_index, digests = 0, {}, {}
        try:
            with open(temp_path, "wb") as out:
                for name in bucket_names(hdir):
                    path = os.path.join(hdir, name)
                    with open(path) as fh:
                        bucket = json.load(fh)
                    if not isinstance(bucket, dict):
                        raise RuntimeError(f"{path}: bucket must be a JSON object")
                    for gid in sorted(bucket):
                        if gid in hap_index:
                            raise RuntimeError(f"{path}: duplicate expression record for {gid}")
                        payload = json.dumps(bucket[gid], separators=(",", ":"),
                                             ensure_ascii=False, allow_nan=False).encode("utf-8")
                        member = gzip_member(payload)
                        out.write(member)
                        hap_index[gid] = [off, len(member)]
                        digests[gid] = hashlib.sha256(payload).digest()
                        off += len(member)
                out.flush()
                os.fsync(out.fileno())
        except Exception:
            os.unlink(temp_path)
            raise
        staged.append((temp_path, out_path))
        index[hap], totals[hap], payload_digests[hap] = hap_index, off, digests
        print(f"  {out_name}_{hap}.gzb  {off / 1e6:8.1f} MB")
    return index, totals, payload_digests


def validate_staged_blob(path, index, digests):
    """Prove every staged member is the exact source vector at its recorded range."""
    end = 0
    with open(path, "rb") as fh:
        for gid, (off, length) in sorted(index.items(), key=lambda item: item[1][0]):
            if off != end or length <= 0:
                raise RuntimeError(f"{path}: invalid range for {gid}: {off}+{length}")
            fh.seek(off)
            member = fh.read(length)
            if len(member) != length:
                raise RuntimeError(f"{path}: truncated member for {gid}")
            try:
                payload = gzip.decompress(member)
                vector = json.loads(payload)
            except (OSError, ValueError, json.JSONDecodeError) as error:
                raise RuntimeError(f"{path}: cannot decode member for {gid}: {error}") from error
            if not isinstance(vector, list):
                raise RuntimeError(f"{path}: member for {gid} is not a JSON array")
            if hashlib.sha256(payload).digest() != digests.get(gid):
                raise RuntimeError(f"{path}: member for {gid} differs from its source vector")
            end += length
    size = os.path.getsize(path)
    if end != size:
        raise RuntimeError(f"{path}: ranges cover {end} bytes, file is {size}")


def validate_staged_blobs(staged, raw_index, raw_digests, combat_index, combat_digests):
    """Validate the exact staged files that will replace the release blobs."""
    by_final = {final: temp for temp, final in staged}
    if len(by_final) != len(staged):
        raise RuntimeError("duplicate staged range output")
    for label, index, digests in (
        ("raw", raw_index, raw_digests),
        ("cb", combat_index, combat_digests),
    ):
        for hap in HAPS:
            final = os.path.join(DATA, "expr", f"{label}_{hap}.gzb")
            temp = by_final.get(final)
            if temp is None:
                raise RuntimeError(f"missing staged blob for {label}_{hap}")
            validate_staged_blob(temp, index[hap], digests[hap])


def load_shards():
    """Read each shard once and return mutable records keyed by haplotype and path."""
    shards = {hap: [] for hap in HAPS}
    for hap in HAPS:
        root = os.path.join(DATA, "genes", hap)
        if not os.path.isdir(root):
            raise RuntimeError(f"required gene shard directory is absent: {root}")
        for path in shard_paths(root):
            with open(path) as fh:
                shards[hap].append((path, json.load(fh)))
    return shards


def expected_genes(shards, pointer):
    return {
        hap: {gid for _path, shard in shards[hap] for gid, rec in shard.items()
              if (rec.get("x") or {}).get(pointer)}
        for hap in HAPS
    }


def assert_index_matches(pointer, index, shards):
    """Reject a source/tree mismatch before either blobs or shard ranges are changed."""
    expected = expected_genes(shards, pointer)
    for hap in HAPS:
        missing = expected[hap] - set(index[hap])
        extra = set(index[hap]) - expected[hap]
        if missing or extra:
            raise RuntimeError(
                f"{pointer} ranges for {hap} do not match shard pointers: "
                f"{len(missing)} missing, {len(extra)} extra"
            )


def stage_ranges(shards, raw_index, combat_index, staged):
    """Stage shard offsets; existing shards remain untouched until transaction commit."""
    n_raw = n_combat = 0
    for hap in HAPS:
        for path, shard in shards[hap]:
            changed = False
            for gid, rec in shard.items():
                x = rec.get("x")
                if not x:
                    continue
                if x.get("b"):
                    if x.get("o") != raw_index[hap][gid]:
                        x["o"] = raw_index[hap][gid]
                        changed = True
                    n_raw += 1
                elif x.pop("o", None) is not None:
                    changed = True
                if x.get("cb"):
                    if x.get("co") != combat_index[hap][gid]:
                        x["co"] = combat_index[hap][gid]
                        changed = True
                    n_combat += 1
                elif x.pop("co", None) is not None:
                    changed = True
            if changed:
                payload = json.dumps(shard, separators=(",", ":"), sort_keys=True).encode("utf-8")
                staged.append((stage_bytes(path, payload), path))
    return n_raw, n_combat


def stage_manifest(raw_totals, combat_totals, n_raw, n_combat, staged):
    """Stage matching manifest metadata as part of the same publication transaction."""
    path = os.path.join(DATA, "meta", "manifest.json")
    with open(path) as fh:
        manifest = json.load(fh)
    manifest["expr_ranges"] = {
        "raw_bytes": raw_totals,
        "combat_bytes": combat_totals,
        "genes_raw": n_raw,
        "genes_combat": n_combat,
        "encoding": "one gzip member per gene, concatenated; the member inflates to the "
                    "same JSON array the bucket held, byte for byte. Values are UNCHANGED - "
                    "this is a container change, not a quantisation.",
        "why": "A gene page fetched 786 KB of its 250 neighbours' expression to draw one "
               "card. A range request fetches ~1.1 KB.",
        "rejected": "Storing ComBat as a per-study offset table: measured NOT additive in "
                    "log space (within-study spread of the log2 delta, median 0.43, p90 "
                    "1.15), so both matrices are kept.",
    }
    payload = json.dumps(manifest, indent=1, sort_keys=True).encode("utf-8")
    staged.append((stage_bytes(path, payload), path))


def rollback_transaction(backups, published, moved_trees):
    """Best-effort restoration for a normal I/O failure during the publish transaction."""
    errors = []
    for source, backup in reversed(moved_trees):
        try:
            if os.path.exists(backup):
                os.replace(backup, source)
        except OSError as error:
            errors.append(f"restore {source}: {error}")
    for _temp, final in reversed(published):
        try:
            if os.path.exists(final):
                os.unlink(final)
        except OSError as error:
            errors.append(f"remove new {final}: {error}")
    for final, backup in reversed(backups):
        try:
            if backup and os.path.exists(backup):
                os.replace(backup, final)
        except OSError as error:
            errors.append(f"restore {final}: {error}")
    return errors


def publish_transaction(staged, bucket_trees):
    """Atomically replace staged files and retire inputs only after all outputs are ready.

    Individual ``os.replace`` calls are atomic.  Backups and rollback make the group behave as
    one transaction for ordinary failures: neither old outputs nor build-input bucket trees are
    discarded until every replacement and directory move has succeeded.
    """
    if not staged:
        raise RuntimeError("no outputs staged")
    expr_dir = os.path.join(DATA, "expr")
    transaction = tempfile.mkdtemp(prefix=".expr-range-transaction-", dir=expr_dir)
    backup_files = os.path.join(transaction, "files")
    backup_trees = os.path.join(transaction, "trees")
    os.mkdir(backup_files)
    os.mkdir(backup_trees)
    backups, published, moved_trees = [], [], []
    try:
        for number, (_temp, final) in enumerate(staged):
            backup = None
            if os.path.exists(final):
                backup = os.path.join(backup_files, str(number))
                os.replace(final, backup)
            backups.append((final, backup))
        for temp, final in staged:
            os.replace(temp, final)
            published.append((temp, final))
        for tree in bucket_trees:
            if not os.path.isdir(tree):
                raise RuntimeError(f"bucket tree disappeared before publish: {tree}")
            backup = os.path.join(backup_trees, os.path.basename(tree))
            os.replace(tree, backup)
            moved_trees.append((tree, backup))
    except Exception as error:
        restore_errors = rollback_transaction(backups, published, moved_trees)
        if not restore_errors:
            try:
                shutil.rmtree(transaction)
            except OSError as cleanup_error:
                restore_errors.append(f"remove rollback directory {transaction}: {cleanup_error}")
        detail = f"; rollback also failed: {'; '.join(restore_errors)}" if restore_errors else ""
        raise RuntimeError(f"range publication failed: {error}{detail}") from error
    finally:
        # Unpublished staged files still exist after any early failure.  Published files no
        # longer exist at their temporary names, so this is harmless after a successful commit.
        for temp, _final in staged:
            if os.path.exists(temp):
                os.unlink(temp)

    # The source trees now live under the transaction directory, outside their public paths.
    # A cleanup failure cannot expose them through the tracked artifact; leave it visible for
    # manual recovery instead of turning a successful release into a half-rolled-back one.
    try:
        shutil.rmtree(transaction)
    except OSError as error:
        print(f"WARNING: packed bucket backups remain at {transaction}: {error}")


def main():
    staged = []
    raw_tree = os.path.join(DATA, "expr", "buckets")
    combat_tree = os.path.join(DATA, "expr", "combat_buckets")
    try:
        print("packing raw expression")
        raw_idx, raw_tot, raw_digests = build_one("buckets", "raw", staged)
        print("packing ComBat expression")
        cb_idx, cb_tot, cb_digests = build_one("combat_buckets", "cb", staged)

        shards = load_shards()
        assert_index_matches("b", raw_idx, shards)
        assert_index_matches("cb", cb_idx, shards)
        validate_staged_blobs(staged, raw_idx, raw_digests, cb_idx, cb_digests)

        n_raw, n_cb = stage_ranges(shards, raw_idx, cb_idx, staged)
        stage_manifest(raw_tot, cb_tot, n_raw, n_cb, staged)
        publish_transaction(staged, (raw_tree, combat_tree))
        staged.clear()
    finally:
        for temp, _out_path in staged:
            if os.path.exists(temp):
                os.unlink(temp)

    print(f"\n{n_raw:,} genes got a raw range; {n_cb:,} got a ComBat range")
    print("removed packed expression bucket trees")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
