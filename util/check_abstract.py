#!/usr/bin/env python3
"""Check the abstract of paper/main.tex against the arXiv 1920-character limit.

Counts the collapsed body between \\begin{abstract} and \\end{abstract},
including TeX commands and squeezing whitespace to single spaces. Exits
nonzero if that body is over the limit, uses a preamble macro or a citation
command, or opens with \\cite.

    python util/check_abstract.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEX = ROOT / "paper" / "main.tex"
LIMIT = 1920

CITE_COMMANDS = {
    "cite",
    "citealp",
    "citealt",
    "citeauthor",
    "citep",
    "citet",
    "citeyear",
    "citeyearpar",
    "nocite",
}

MACRO_DEF = re.compile(r"\\(?:newcommand|renewcommand|providecommand|def)\s*\*?\{?\\([A-Za-z@]+)")
ABSTRACT = re.compile(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", re.S)
COMMAND = re.compile(r"\\([A-Za-z@]+)")


def collapsed_body(tex: str) -> str:
    match = ABSTRACT.search(tex)
    if not match:
        raise SystemExit(f"no abstract in {TEX}")
    return re.sub(r"\s+", " ", match.group(1)).strip()


def defined_macros(tex: str) -> set[str]:
    preamble = tex.split(r"\begin{document}", 1)[0]
    return set(MACRO_DEF.findall(preamble))


def main() -> int:
    tex = TEX.read_text(encoding="utf-8")
    body = collapsed_body(tex)
    used = set(COMMAND.findall(body))
    forbidden = sorted(used & (defined_macros(tex) | CITE_COMMANDS))
    leads = body.startswith("\\cite")
    print(f"abstract characters: {len(body)} / {LIMIT}")
    print("forbidden macros: " + (", ".join(forbidden) if forbidden else "none"))
    print("opens with a citation command: " + ("yes" if leads else "no"))
    if len(body) <= LIMIT and not forbidden and not leads:
        return 0
    if len(body) > LIMIT:
        print(f"abstract is {len(body) - LIMIT} characters over the limit")
    if forbidden:
        print("abstract uses forbidden macros")
    if leads:
        print("abstract opens with a citation command")
    return 1


if __name__ == "__main__":
    sys.exit(main())
