# Why 8,185 TESS RR Lyrae say "no data"

Measured 2026-09-16 against the store and a live MAST sample.

## It is not a search bug

`stage1_summary` holds 10,348 error rows, and 8,185 of them are one cell:
`survey='tess'`, `star_class='rrlyrae'`, `error='no data'`. `load_data` falls
back to `download_tess`, which runs one unfiltered `search_lightcurve` and then
applies the author/cadence preference locally, so "no data" means the search
came back empty or nothing in it could be read.

Ten of the 8,185, drawn at random, were queried live:

```
lk.search_lightcurve("TIC ...", mission="TESS")  ->  0 rows, all ten
```

Four randomly drawn *successful* TESS RR Lyrae, for contrast, returned 1-6 rows
each, and **every one of them was 1800 s** — TASOC, GSFC-ELEANOR-LITE, T16,
TGLC. Those are cycle 1-2 FFI pipelines.

So the failures are real: these stars have no light-curve product at MAST at
all. They are not on the 2-min target list and they are too faint for QLP
(median Tmag 14.6, running to 15.5). Nothing in the download path recovers
them, because there is nothing to download.

## What the TESS RR Lyrae arm is actually made of

Effective cadence inferred from `points_per_cycle` over 4,000 stage-1 rows:

| cadence | stars | dense mode | median O–C points | median baseline |
|---|---|---|---|---|
| ~120 s | 172 | 172 | 64 | 27 d |
| ~200 s | 452 | 452 | 146 | 1509 d |
| ~600 s | 353 | 351 | 179 | 1850 d |
| ~1800 s | 3,023 | 497 | 34 | 27 d |

Three quarters of the arm is a **single 27-day sector at 30-minute cadence**,
and five sixths of that cannot reach dense mode at all
(`DENSE_MIN_PTS_PER_CYCLE = 30`, and P = 0.6 d at 1800 s gives 28). The part of
the arm carrying multi-year baselines and 150+ per-cycle timings is the ~800
stars with 200-600 s products.

## The FFIs have them

The same stars do have FFI coverage. `search_tesscut` on three of the failures:

```
TIC 13813938   2 sectors (54, 92)
TIC 30858850   6 sectors (14, 40, 41, 54, 74, 81)
TIC 456623745  3 sectors (55, 81, 82)
```

Sectors 40+ are 200 s FFIs, which is *better* than what most of the working
arm has.

## Prototype: one star, end to end

TIC 30858850 (Tmag 14.25, P = 0.3398 d), sector 81, 9x9 TESScut cutout:

- download 13 s, 18.7 MB, 11,033 frames at 200 s
- background from the faintest 40 % of pixels, per frame
- aperture placed by WCS, not at the cutout centre — **this is the whole
  game**. The cutout centre is not the star: the field's brightest pixel is a
  neighbour at ~1500 e⁻/s against the target's ~300, and a 9-pixel aperture
  gave a 1.3 % amplitude (an RRab is 20-100 %) and 1286 s timings.

| aperture | amplitude | median timing err | O–C rms |
|---|---|---|---|
| 9x9 box (centred on cutout) | 1.27 % | 1286 s | 795 s |
| 1 px at the WCS position | 7.98 % | 235 s | 459 s |
| **r < 1.0 px** | **8.20 %** | **126 s** | **174 s** |
| r < 1.5 px | 8.83 % | 207 s | 309 s |
| r < 2.0 px | 5.15 % | 491 s | 403 s |
| gaussian sigma = 0.8 px | 7.64 % | 191 s | 236 s |

Through the real `run_target`: dense mode, 147 points per cycle, 66 O–C points
from the one sector, period recovered to 1.4e-5 d of the Gaia value.

8.2 % against an intrinsic ~30 % still means ~4x dilution, so a PRF fit or
difference imaging has more to give, and the TIC `contratio` would say in
advance which stars are hopeless.

## Worth doing?

6,008 of the 8,185 already have a usable ASAS-SN or ZTF result, so this is not
mainly a yield fix — 2,177 stars have nothing anywhere. It is a *regime* fix:
ground-based arms can only give block-averaged timing and a jitter statistic,
while this gives 126 s per-cycle timing over 66+ cycles per sector on multi-year
baselines. And it applies just as much to the 3,023 stars currently stuck on a
single 30-min sector, which is the larger population.

Cost, streaming the cubes and keeping only the light curves: ~6 MB per
star-sector at 5x5, ~13 s each, call it 33,000 star-sectors for the RR Lyrae —
a few hundred GB of transfer over a multi-day background run, and TESScut is
rate-limited. A `contratio` cut would trim it first.
