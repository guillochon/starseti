#!/usr/bin/env python3
"""Check the title and abstract of paper/main.tex for arXiv.

Counts the collapsed body between \\begin{abstract} and \\end{abstract},
including TeX commands and squeezing whitespace to single spaces, against the
1920-character limit. Exits nonzero if that body is over the limit, if the
abstract or the \\title uses a preamble macro or a citation command, or if the
abstract opens with \\cite.

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
TITLE = re.compile(r"\\title\s*\{")
COMMAND = re.compile(r"\\([A-Za-z@]+)")


def collapsed_body(tex: str) -> str:
    match = ABSTRACT.search(tex)
    if not match:
        raise SystemExit(f"no abstract in {TEX}")
    return re.sub(r"\s+", " ", match.group(1)).strip()


def title(tex: str) -> str:
    """Body of \\title{...}, with nested braces."""
    match = TITLE.search(tex)
    if not match:
        raise SystemExit(f"no title in {TEX}")
    depth, i = 1, match.end()
    while i < len(tex) and depth:
        depth += {"{": 1, "}": -1}.get(tex[i], 0)
        i += 1
    return re.sub(r"\s+", " ", tex[match.end() : i - 1]).strip()


def defined_macros(tex: str) -> set[str]:
    preamble = tex.split(r"\begin{document}", 1)[0]
    return set(MACRO_DEF.findall(preamble))


def main() -> int:
    tex = TEX.read_text(encoding="utf-8")
    banned = defined_macros(tex) | CITE_COMMANDS
    body = collapsed_body(tex)
    forbidden = sorted(set(COMMAND.findall(body)) & banned)
    in_title = sorted(set(COMMAND.findall(title(tex))) & banned)
    leads = body.startswith("\\cite")
    print(f"abstract characters: {len(body)} / {LIMIT}")
    print("forbidden macros: " + (", ".join(forbidden) if forbidden else "none"))
    print("forbidden macros in the title: " + (", ".join(in_title) if in_title else "none"))
    print("opens with a citation command: " + ("yes" if leads else "no"))
    if len(body) <= LIMIT and not forbidden and not in_title and not leads:
        return 0
    if len(body) > LIMIT:
        print(f"abstract is {len(body) - LIMIT} characters over the limit")
    if forbidden:
        print("abstract uses forbidden macros")
    if in_title:
        print("title uses forbidden macros")
    if leads:
        print("abstract opens with a citation command")
    return 1


if __name__ == "__main__":
    sys.exit(main())
