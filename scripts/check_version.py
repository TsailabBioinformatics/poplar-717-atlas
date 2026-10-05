#!/usr/bin/env python3
"""
Assert the four places that carry a version number still agree.

This exists because they did not. `data_version` was READ by the validator, the site footer
and the About page, and documented in CHANGELOG.md as "stamped into manifest.json" -- but no
build script ever wrote it. It survived only because it had been hand-added to the committed
manifest after the last rebuild. The next `build_all.sh` would have dropped it, and nothing
would have failed: the footer hides the chip when the key is absent (`if (f && man.data_version)`)
and the About page falls back to an em dash. A silently missing version on a citable dataset.

Runs anywhere -- it reads only committed repo files, so CI catches drift off-cluster.
"""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
fails = []

def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}{'  — ' + detail if detail else ''}")
    if not ok:
        fails.append(name)

def main():
    version = open(os.path.join(ROOT, "VERSION")).read().strip()
    print(f"\nVERSION = {version}\n")
    check("VERSION is a bare semver", bool(re.fullmatch(r"\d+\.\d+\.\d+", version)), version)

    man = json.load(open(os.path.join(ROOT, "data", "meta", "manifest.json")))
    check("manifest carries data_version", "data_version" in man,
          "absent — build_data.py did not stamp it" if "data_version" not in man else "")
    check("manifest data_version matches VERSION", man.get("data_version") == version,
          f"manifest {man.get('data_version')!r} vs VERSION {version!r}")

    cff = open(os.path.join(ROOT, "CITATION.cff")).read()
    m = re.search(r'^version:\s*"?([^"\n]+)"?\s*$', cff, re.M)
    cff_ver = m.group(1).strip() if m else None
    check("CITATION.cff version matches VERSION", cff_ver == version,
          f"CITATION.cff {cff_ver!r} vs VERSION {version!r}")

    m = re.search(r'^date-released:\s*"?([\d-]+)"?\s*$', cff, re.M)
    cff_date = m.group(1) if m else None

    # Top-most VERSIONED heading, e.g. "## 0.3.10 — 2026-09-07 — Expression scrubber".
    # Releases that changed only the site carry a plain "## title" heading and are skipped
    # here by design -- see the convention note at the top of CHANGELOG.md.
    chg = open(os.path.join(ROOT, "CHANGELOG.md")).read()
    m = re.search(r'^##\s*(\d+\.\d+\.\d+)\s*[—-]\s*(\d{4}-\d{2}-\d{2})\b', chg, re.M)
    chg_ver, chg_date = (m.group(1), m.group(2)) if m else (None, None)
    check("CHANGELOG's newest versioned entry is VERSION", chg_ver == version,
          f"CHANGELOG {chg_ver!r} vs VERSION {version!r}")
    check("CITATION.cff date-released matches that entry", cff_date == chg_date,
          f"CITATION.cff {cff_date!r} vs CHANGELOG {chg_date!r}")

    # index.html carries the build version the shipped JavaScript belongs to; app.js compares it
    # against the manifest at boot and warns a reader holding a cached build. If it drifts from
    # VERSION the warning either never fires or fires forever, so it is held equal here.
    html = open(os.path.join(ROOT, "index.html")).read()
    m = re.search(r'<meta name="atlas-build-version" content="([^"]+)">', html)
    check("index.html build version matches VERSION", bool(m) and m.group(1) == version,
          f"index.html '{m.group(1) if m else 'MISSING'}' vs VERSION '{version}'")

    print(f"\n{len(fails)} failed\n" if fails else "\nall version references agree\n")
    if fails:
        print("To release: edit VERSION, add a CHANGELOG entry, sync CITATION.cff, rebuild "
              "(or restamp the manifest), then re-run this check.")
    return 1 if fails else 0

if __name__ == "__main__":
    sys.exit(main())
