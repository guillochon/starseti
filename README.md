# Listening to the Standard Candles

Manuscript, bibliography, and figure scripts for the survey-scale search for
phase-coded beacons in Cepheids, RR Lyrae, and beta Cephei stars.

## Build

```bash
cd paper
latexmk -pdf
```

builds both typesettings of `paper/main.tex` in the Acta Astronautica
(`elsarticle`) format: `submit.pdf`, single column for journal submission, and
`preprint.pdf`, two columns. `elsarticle` and `newtx` come with TeX Live and
MiKTeX. `python paper/arxiv_submit.py` zips the two-column version for arXiv.
See `paper/README.md`.

## Figures

```bash
uv sync
uv run python paper/figures/make_geometry.py
uv run python paper/figures/make_synthetic.py
uv run python paper/figures/make_completeness.py
uv run python paper/figures/make_candidates.py
```

`make_geometry.py` needs only matplotlib. `make_synthetic.py` runs the timing
code in `src/` and takes a few minutes. `make_completeness.py` reads
`data/injection.sqlite`, and `make_candidates.py` reads `data/cepheid.sqlite`.
Those stores are products of the search and are not in this snapshot.

The derivations behind the beam powers are in `research/`.

## Checks

`python util/check_abstract.py` enforces the 1920-character abstract limit.
`python util/ads_validate.py --keys` checks that every citation key resolves.
The ADS field check needs `ADS_TOKEN` in `.env`, which is not committed.
