#!/usr/bin/env python3
"""Put the keys inside every multi-key \\cite command of main.tex in year order.

Years come from the .bib files next to this script; ties keep the order the
author wrote.  Run from anywhere:

    python paper/sort_citations.py            # rewrites main.tex in place
    python paper/sort_citations.py --check    # exit 1 if anything is out of order
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEX = HERE / "main.tex"
BIBS = sorted(HERE.glob("*.bib"))

ENTRY = re.compile(r"@\w+\{\s*([^,\s]+)\s*,(.*?)\n\}", re.S)
YEAR = re.compile(r"\byear\s*=\s*\{?\"?(\d{4})", re.I)
CITE = re.compile(r"(\\cite[tp]?\*?(?:\[[^\]]*\]){0,2}\{)([^}]*)(\})")


def years() -> dict[str, int]:
    out: dict[str, int] = {}
    for bib in BIBS:
        for m in ENTRY.finditer(bib.read_text(encoding="utf-8")):
            y = YEAR.search(m.group(2))
            if y:
                out[m.group(1)] = int(y.group(1))
    return out


def main() -> int:
    check = "--check" in sys.argv
    yr = years()
    text = TEX.read_text(encoding="utf-8")
    changed = 0
    missing: set[str] = set()

    def fix(m: re.Match) -> str:
        nonlocal changed
        keys = [k.strip() for k in m.group(2).split(",")]
        if len(keys) < 2:
            return m.group(0)
        for k in keys:
            if k not in yr:
                missing.add(k)
        ordered = sorted(keys, key=lambda k: yr.get(k, 9999))  # stable: ties keep author order
        if ordered != keys:
            changed += 1
        return m.group(1) + ", ".join(ordered) + m.group(3)

    new = CITE.sub(fix, text)
    if missing:
        print("no year found for:", ", ".join(sorted(missing)))
    if check:
        print(f"{changed} citation(s) out of chronological order")
        return 1 if changed else 0
    if changed:
        TEX.write_text(new, encoding="utf-8")
    print(f"reordered {changed} citation(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
