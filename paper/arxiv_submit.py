#!/usr/bin/env python3
"""
Zip the files arXiv needs to compile the two-column preprint.

main.tex picks its layout from \\papermode, which preprint.tex sets before
\\input{main}. arXiv compiles a single top-level file, so the zip's main.tex is
paper/main.tex with \\def\\papermode{preprint} prepended: the same source,
typeset as preprint.pdf. Every figure named by \\includegraphics keeps its path
(figures/...). The bibliography goes in as main.bbl, copied from the
preprint.bbl of a fresh build in a temporary directory (so a viewer holding
paper/preprint.pdf open cannot block it), and arXiv needs neither BibTeX nor the
.bib files. elsarticle.cls, its .bst and newtx are in arXiv's TeX Live and are
not bundled. The compiled PDF and aux files are omitted.

    python paper/arxiv_submit.py              # build preprint, write starseti-arxiv.zip
    python paper/arxiv_submit.py --out PATH   # choose the zip path
    python paper/arxiv_submit.py --no-build   # reuse the existing preprint.bbl
    python paper/arxiv_submit.py --dry-run    # list files, do not write
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

PAPER = Path(__file__).resolve().parent
ROOT = PAPER.parent
TEX_NAME = "main.tex"
BBL_NAME = "preprint.bbl"
BIBS = ["references.bib", "references-seti.bib"]
DEFAULT_ZIP = ROOT / "starseti-arxiv.zip"

PREAMBLE = (
    "% arXiv build of the two-column preprint layout (see preprint.tex in the source repository).\n"
    "\\def\\papermode{preprint}\n"
)

# Optional [...] before the brace, so \includegraphics[width=...]{...} matches.
GRAPHICS_CMD = re.compile(r"\\includegraphics\s*(?:\[[^\]]*\])?\s*\{([^}]+)\}")


def strip_comments(text):
    """Drop TeX comments. A % preceded by an odd number of backslashes is kept."""
    kept = []
    for line in text.splitlines(keepends=True):
        i = 0
        newline = "\n" if line.endswith("\n") else ""
        body = line[:-1] if newline else line
        while True:
            j = body.find("%", i)
            if j < 0:
                kept.append(body + newline)
                break
            backslashes = 0
            k = j - 1
            while k >= 0 and body[k] == "\\":
                backslashes += 1
                k -= 1
            if backslashes % 2 == 1:
                i = j + 1
                continue
            kept.append(body[:j] + newline)
            break
    return "".join(kept)


def figures(tex_text):
    """Figure paths relative to paper/, in first-use order."""
    ordered = []
    for match in GRAPHICS_CMD.finditer(strip_comments(tex_text)):
        name = match.group(1).strip()
        if name and name not in ordered:
            ordered.append(name)
    return ordered


def build(outdir):
    """Typeset the preprint into outdir, so a PDF viewer holding paper/preprint.pdf cannot block it."""
    cmd = ["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error",
           f"-outdir={outdir}", "preprint.tex"]
    print("Running", " ".join(cmd))
    # BibTeX runs in outdir, so tell it where the .bib files are.
    env = dict(os.environ, BIBINPUTS=f"{PAPER}{os.pathsep}{os.environ.get('BIBINPUTS', '')}")
    result = subprocess.run(cmd, cwd=PAPER, capture_output=True, text=True, env=env)
    if result.returncode != 0:
        print(result.stdout[-3000:], file=sys.stderr)
        print(f"latexmk failed; see {outdir / 'preprint.log'}", file=sys.stderr)
        sys.exit(1)
    return outdir / BBL_NAME


def check_bbl():
    bbl = PAPER / BBL_NAME
    if not bbl.is_file():
        print(f"Missing paper/{BBL_NAME}; run without --no-build", file=sys.stderr)
        sys.exit(1)
    newer = [n for n in [TEX_NAME, *BIBS] if (PAPER / n).stat().st_mtime > bbl.stat().st_mtime]
    if newer:
        print(f"paper/{BBL_NAME} is older than {', '.join(newer)}; run without --no-build",
              file=sys.stderr)
        sys.exit(1)
    return bbl


def collect(bbl):
    """(zip name, bytes) pairs, the .tex first."""
    tex = (PAPER / TEX_NAME).read_text(encoding="utf-8")
    names = figures(tex)
    missing = [name for name in names if not (PAPER / name).is_file()]
    if missing:
        print("Missing from paper/:", file=sys.stderr)
        for name in missing:
            print(f"  {name}", file=sys.stderr)
        sys.exit(1)
    items = [(TEX_NAME, (PREAMBLE + tex).encode("utf-8"))]
    bbl_bytes = bbl.read_bytes()
    if b"\\bibitem" not in bbl_bytes:
        print(f"{bbl} has no \\bibitem entries; BibTeX did not find the references", file=sys.stderr)
        sys.exit(1)
    items.append(("main.bbl", bbl_bytes))
    items += [(name, (PAPER / name).read_bytes()) for name in names]
    return items


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_ZIP,
        help=f"zip path (default: {DEFAULT_ZIP.name} in the repository root)",
    )
    parser.add_argument(
        "--no-build",
        action="store_true",
        help=f"skip latexmk and reuse paper/{BBL_NAME} (refused if it is stale)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="list files and sizes, do not write",
    )
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmp:
        bbl = check_bbl() if args.no_build else build(Path(tmp))
        items = collect(bbl)
    width = max(len(name) for name, _ in items)
    total = 0
    for name, data in items:
        total += len(data)
        print(f"  {name:{width}s} {len(data) / 1024:8.1f} KB")
    print(f"  {'total':{width}s} {total / 1024:8.1f} KB")

    if args.dry_run:
        return

    out = args.out.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(out, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name, data in items:
            zf.writestr(name, data)
    print(f"Wrote {out} ({out.stat().st_size / 1024:.1f} KB)")


if __name__ == "__main__":
    main()
