#!/usr/bin/env python3
"""Focused fixture for transactional range-expression packing.

Run directly with ``python3 scripts/test_build_expr_ranges.py``.  It proves a successful pack
retires the two redundant bucket trees, and an injected publication failure restores every old
output and keeps both input trees intact.
"""
from __future__ import annotations

import gzip
import importlib.util
import json
import os
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "build_expr_ranges.py"
READER = ROOT / "scripts" / "expr_range_reader.py"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_builder():
    return load_module(SCRIPT, "build_expr_ranges_fixture")


def load_reader():
    return load_module(READER, "expr_range_reader_fixture")


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, separators=(",", ":")), encoding="utf-8")


def fixture(root: Path):
    data = root / "data"
    expr = data / "expr"
    expected = {}
    for hap, gene, raw, combat in (
        ("hap1", "PtXaTreH.01G000100", [1.0, 2.5, 3.0], [1.1, 2.6, 3.1]),
        ("hap2", "PtXaAlbH.01G000100", [4.0, 5.5, 6.0], [4.1, 5.6, 6.1]),
    ):
        write_json(expr / "buckets" / hap / "0.json", {gene: raw})
        write_json(expr / "combat_buckets" / hap / "0.json", {gene: combat})
        write_json(data / "genes" / hap / "Chr01" / "0.json", {
            gene: {"id": gene, "x": {"b": f"{hap}/0", "cb": f"{hap}/0"}},
        })
        for prefix in ("raw", "cb"):
            (expr / f"{prefix}_{hap}.gzb").write_bytes(f"old-{prefix}-{hap}".encode())
        expected[hap] = {"gene": gene, "raw": raw, "combat": combat}
    write_json(data / "meta" / "manifest.json", {"data_version": "fixture"})
    return data, expected


def read_member(blob: Path, offset_length):
    offset, length = offset_length
    with blob.open("rb") as fh:
        fh.seek(offset)
        return json.loads(gzip.decompress(fh.read(length)))


def assert_no_temps(data: Path) -> None:
    leftovers = list((data / "expr").glob(".*.tmp")) + list((data / "expr").glob(".expr-range-transaction-*"))
    assert not leftovers, f"leftover transaction files: {leftovers}"


def pointers_from_shards(data: Path):
    pointers = {}
    for hap in ("hap1", "hap2"):
        shard = json.loads((data / "genes" / hap / "Chr01" / "0.json").read_text())
        for gene, record in shard.items():
            x = record["x"]
            pointers[gene] = {"hap": hap, "range": x.get("co"), "bucket": x.get("cb")}
    return pointers


def test_success(module, reader) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        data, expected = fixture(Path(tmp))
        module.DATA = str(data)
        legacy, storage = reader.load_combat_vectors(str(data), pointers_from_shards(data))
        assert storage == {"range": 0, "legacy": 2, "legacy_buckets": 2}
        for hap, values in expected.items():
            assert legacy[values["gene"]] == values["combat"]
        assert module.main() == 0
        assert not (data / "expr" / "buckets").exists()
        assert not (data / "expr" / "combat_buckets").exists()
        assert_no_temps(data)
        manifest = json.loads((data / "meta" / "manifest.json").read_text())
        assert manifest["expr_ranges"]["genes_raw"] == 2
        assert manifest["expr_ranges"]["genes_combat"] == 2
        packed, storage = reader.load_combat_vectors(str(data), pointers_from_shards(data))
        assert storage == {"range": 2, "legacy": 0, "legacy_buckets": 0}
        for hap, values in expected.items():
            assert packed[values["gene"]] == values["combat"]
            shard = json.loads((data / "genes" / hap / "Chr01" / "0.json").read_text())
            record = shard[values["gene"]]
            assert read_member(data / "expr" / f"raw_{hap}.gzb", record["x"]["o"]) == values["raw"]
            assert read_member(data / "expr" / f"cb_{hap}.gzb", record["x"]["co"]) == values["combat"]


def test_publication_failure_rolls_back(module) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        data, _expected = fixture(Path(tmp))
        module.DATA = str(data)
        old_outputs = {
            path: path.read_bytes()
            for path in (data / "expr").glob("*.gzb")
        }
        old_shards = {
            path: path.read_bytes()
            for path in (data / "genes").glob("**/*.json")
        }
        old_manifest = (data / "meta" / "manifest.json").read_bytes()
        raw_hap1 = str(data / "expr" / "raw_hap1.gzb")
        real_replace = module.os.replace

        def fail_first_publish(source, destination):
            if str(destination) == raw_hap1 and str(source).endswith(".tmp"):
                raise OSError("fixture: injected publish failure")
            return real_replace(source, destination)

        module.os.replace = fail_first_publish
        try:
            try:
                module.main()
            except RuntimeError as error:
                assert "range publication failed" in str(error)
            else:
                raise AssertionError("fixture expected publication failure")
        finally:
            module.os.replace = real_replace

        assert (data / "expr" / "buckets").is_dir()
        assert (data / "expr" / "combat_buckets").is_dir()
        assert {path: path.read_bytes() for path in (data / "expr").glob("*.gzb")} == old_outputs
        assert {path: path.read_bytes() for path in (data / "genes").glob("**/*.json")} == old_shards
        assert (data / "meta" / "manifest.json").read_bytes() == old_manifest
        assert_no_temps(data)


def main() -> int:
    module = load_builder()
    reader = load_reader()
    test_success(module, reader)
    test_publication_failure_rolls_back(module)
    print("range-expression transactional fixture: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
