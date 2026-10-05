#!/usr/bin/env python3
"""Keep index.html's <link rel="modulepreload"> list equal to js/app.js's static import graph.

Why: the site has no bundler, so the browser discovers modules one import level at a time --
app.js, then registry.js, then each module, then core helpers -- one network round trip per
level. On a 170 ms mobile link that chain cost ~0.7 s before any data was requested. Listing
every module in the head lets the browser fetch them all in parallel from the first byte.

The list is GENERATED, never hand-edited, and `--check` fails when it drifts (a new module
that is not preloaded still works, just slower; a preloaded file that no longer exists is a 404).

  scripts/sync_preload.py          rewrite the block between the markers
  scripts/sync_preload.py --check  exit 1 if the block is stale
"""
import re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INDEX = ROOT / "index.html"
BEGIN, END = "<!-- modulepreload:begin -->", "<!-- modulepreload:end -->"
IMPORT = re.compile(r"""(?:^|\n)\s*import\s+(?:[^'"]*?\s+from\s+)?['"](\.{1,2}/[^'"]+)['"]""")


def graph(entry: Path) -> list[str]:
    seen, stack = set(), [entry]
    while stack:
        f = stack.pop()
        if f in seen or not f.exists():
            continue
        seen.add(f)
        for spec in IMPORT.findall(f.read_text()):
            stack.append((f.parent / spec).resolve())
    return sorted(p.relative_to(ROOT).as_posix() for p in seen if p != entry.resolve())


def block() -> str:
    lines = [f'<link rel="modulepreload" href="{p}">' for p in graph((ROOT / "js/app.js").resolve())]
    return BEGIN + "\n" + "\n".join(lines) + "\n" + END


def main() -> int:
    html = INDEX.read_text()
    want = block()
    if BEGIN in html:
        cur = html[html.index(BEGIN):html.index(END) + len(END)]
    else:
        cur = None
    if "--check" in sys.argv:
        if cur != want:
            print("sync_preload: index.html modulepreload list is stale; run scripts/sync_preload.py")
            return 1
        print(f"sync_preload: {want.count('modulepreload\"')} modules preloaded, list current")
        return 0
    if cur is None:
        html = html.replace('<script type="module" src="js/app.js"></script>',
                            '<script type="module" src="js/app.js"></script>', 1)
        html = html.replace("</head>", want + "\n</head>", 1)
    else:
        html = html.replace(cur, want)
    INDEX.write_text(html)
    print(f"sync_preload: wrote {want.count('modulepreload\"')} modulepreload links")
    return 0


if __name__ == "__main__":
    sys.exit(main())
