#!/usr/bin/env bash
# Move the pin on ms1-dup's numbers of record to a newer commit.
#
# The atlas checks itself against FACTS_v10.json (scripts/check_facts.py) but deliberately does
# NOT keep a copy: this repo's root is published by GitHub Pages, and that file carries values
# held off the public site. sources/facts_v10.pin.json records only the commit and sha256.
#
#   scripts/refresh_facts.sh [ref]      # ref defaults to origin/main
#
# MS1_DUP_ROOT overrides where ms1-dup is. After moving the pin, run check_facts.py and READ
# WHAT CHANGED: a new pin means the owning analysis moved, which is exactly the event this repo
# has historically failed to notice.
set -euo pipefail

REF="${1:-origin/main}"
MS1="${MS1_DUP_ROOT:-/scratch/$USER/ms1-dup}"
SRC="data/81_master_tables_v10_20260920/tables/FACTS_v10.json"
PIN="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/sources/facts_v10.pin.json"

git -C "$MS1" rev-parse --is-inside-work-tree >/dev/null 2>&1 \
  || { echo "!! no ms1-dup checkout at $MS1 (set MS1_DUP_ROOT)"; exit 1; }
git -C "$MS1" fetch --quiet origin || true
COMMIT="$(git -C "$MS1" rev-parse "$REF^{commit}")"
SHA="$(git -C "$MS1" show "$COMMIT:$SRC" | sha256sum | cut -d' ' -f1)"

python3 - "$PIN" "$COMMIT" "$SHA" <<'PY'
import json, sys
path, commit, sha = sys.argv[1:4]
pin = json.load(open(path))
old = pin["ms1_dup_commit"][:9]
pin["ms1_dup_commit"], pin["sha256"] = commit, sha
with open(path, "w") as fh:
    json.dump(pin, fh, indent=1)
    fh.write("\n")
print(f"pin moved {old} -> {commit[:9]}  (sha256 {sha[:16]})")
PY
echo "now run: python3 scripts/check_facts.py   -- and read what changed before committing"
