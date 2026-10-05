"""Read packed expression vectors for source generators, with legacy bucket compatibility.

The public release stores one gzip member per gene in ``cb_<hap>.gzb``.  Source generators may
also be run against an older checkout before the packer stage, so this reader falls back to the
old ``combat_buckets`` tree only when a range pointer or its packed blob is absent.  A malformed
range/blob is a hard error: silently substituting different values would corrupt a derived source.
"""
from __future__ import annotations

import gzip
import json
import math
import os
from collections import defaultdict


def valid_vector(value, label):
    if not isinstance(value, list) or not value:
        raise RuntimeError(f"{label}: expected a non-empty expression array")
    if any(not isinstance(number, (int, float)) or isinstance(number, bool)
           or not math.isfinite(number) for number in value):
        raise RuntimeError(f"{label}: expression array contains a non-finite value")
    return value


def range_vector(handle, offset_length, label):
    if (not isinstance(offset_length, list) or len(offset_length) != 2
            or not all(isinstance(value, int) for value in offset_length)):
        raise RuntimeError(f"{label}: invalid [offset, length] range")
    offset, length = offset_length
    if offset < 0 or length <= 0:
        raise RuntimeError(f"{label}: invalid [offset, length] range")
    handle.seek(offset)
    member = handle.read(length)
    if len(member) != length:
        raise RuntimeError(f"{label}: range is outside the packed blob")
    try:
        return valid_vector(json.loads(gzip.decompress(member)), label)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise RuntimeError(f"{label}: cannot decode packed range: {error}") from error


def load_combat_vectors(data_dir, pointers):
    """Return ``gene_id -> ComBat vector`` from packed ranges or legacy bucket files.

    ``pointers`` maps a gene ID to ``{"hap", "range", "bucket"}``, taken from its shard's
    ``x.co`` and ``x.cb`` fields.  Packed ranges are preferred and read with one file handle per
    haplotype.  Missing ranges/blobs use the legacy bucket only when it is present.
    """
    vectors, legacy = {}, defaultdict(list)
    ranges = defaultdict(list)
    for gene_id, pointer in pointers.items():
        hap = pointer.get("hap")
        bucket = pointer.get("bucket")
        offset_length = pointer.get("range")
        if hap not in {"hap1", "hap2"}:
            raise RuntimeError(f"{gene_id}: invalid haplotype pointer {hap!r}")
        if offset_length is not None:
            ranges[hap].append((gene_id, offset_length, bucket))
        elif bucket:
            legacy[bucket].append(gene_id)
        else:
            raise RuntimeError(f"{gene_id}: has neither x.co nor x.cb")

    range_count = legacy_count = 0
    for hap, entries in ranges.items():
        path = os.path.join(data_dir, "expr", f"cb_{hap}.gzb")
        if not os.path.isfile(path):
            for gene_id, _offset_length, bucket in entries:
                if not bucket:
                    raise RuntimeError(f"{gene_id}: packed blob {path} is absent and no x.cb fallback exists")
                legacy[bucket].append(gene_id)
            continue
        with open(path, "rb") as handle:
            for gene_id, offset_length, _bucket in entries:
                vectors[gene_id] = range_vector(handle, offset_length, f"{path}:{gene_id}")
                range_count += 1

    for bucket, gene_ids in legacy.items():
        path = os.path.join(data_dir, "expr", "combat_buckets", bucket + ".json")
        if not os.path.isfile(path):
            raise RuntimeError(
                f"missing packed ComBat blob or legacy bucket for {bucket}; "
                "run build_expr_ranges.py or restore the matching bucket tree"
            )
        with open(path) as handle:
            values = json.load(handle)
        if not isinstance(values, dict):
            raise RuntimeError(f"{path}: expected a JSON object")
        for gene_id in gene_ids:
            if gene_id not in values:
                raise RuntimeError(f"{path}: missing vector for {gene_id}")
            vectors[gene_id] = valid_vector(values[gene_id], f"{path}:{gene_id}")
            legacy_count += 1

    return vectors, {"range": range_count, "legacy": legacy_count, "legacy_buckets": len(legacy)}
