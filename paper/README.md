# Paper draft

AASTeX 7 (`aastex701`, two-column, no line numbers). Placeholder or provisional
text is set in red with `\colr{...}`.

## Build

```bash
cd paper
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

or `latexmk -pdf main`. MiKTeX and TeX Live both ship `aastex701.cls` and
`aasjournalv7.bst`.

`aasjournalv7.bst` puts a stray space in front of collaboration author names
in the citation labels ("( Gaia Collaboration et al. 2023)") because it always
appends a blank after the (empty) forename. `aasjournalv7collab.bst` is a local
copy with that blank made conditional; `main.tex` uses it. Diff against the
MiKTeX/TeX Live copy to see the two changed lines.

References live in two files: `references.bib` (the idea, physics, surveys,
statistics, software) and `references-seti.bib` (the broader technosignature
context cited in the introduction).

Citations inside one `\citep{}` are kept in chronological order;
`sort_citations.py` reorders every multi-key citation by year and can be re-run
after editing.

## Figures

Every figure is produced by a script in `figures/`; the PDFs are committed, the
scripts are the record of how they were made.

| Script | Output | Needs |
|---|---|---|
| `make_geometry.py` | `geometry.pdf` | matplotlib only |
| `make_completeness.py` | `completeness.pdf` | `data/injection.sqlite` (the `W-<class>-fold` runs) |
| `make_synthetic.py` | `synthetic_examples.pdf` | the `src/` pipeline; runs the three synthetic controls through `stage1_pulse.run_target` (a few minutes) |
| `make_candidates.py` | `candidates.pdf` | `data/cepheid.sqlite` (read-only); TIC 152825245 and TIC 147576267 from the stored O-C, fold, decode and stage-2 products |

Run them from the repository root with the project environment:

```bash
uv run python paper/figures/make_geometry.py
uv run python paper/figures/make_completeness.py
uv run python paper/figures/make_synthetic.py
uv run python paper/figures/make_candidates.py
```

## References

`references.bib` entries were located through the Valency literature index and
written out by hand; arXiv identifiers and DOIs come from that index, while
journal volume and page fields still need to be checked against ADS before
submission (the acknowledgements promise this).
The entries added for the response to review 1 (Gandhi1996, Gandhi1998,
IceCube2021glashow, KM3NeT2025, Smith1971, Edmondson2003) carry volume and
page fields from memory or from Valency's journal-ref strings; Edmondson2003
and Smith1971 are not in Valency at all and must be checked against ADS.
