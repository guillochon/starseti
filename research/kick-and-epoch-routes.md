# Where a real code was being lost, and two routes that read it

*Written 2026-09-21, after the first complete pass over all four surveys.*
The question was not "how do we cut the candidate list" but "is there a real
signal the pipeline cannot see". Most of the answer was upstream of vetting.

## Where the loss was

1. **Amplitude.** Measured thresholds sit 13× or more above the kick ceiling
   of a self-powered transmitter (3e-4) in every class. No route change fixes
   this; only repeats do.
2. **Sparse surveys read no content.** Blocks are 50–135 cycles wide, so a
   block step is the sum of ~100 bits and is Gaussian. A one-bit-per-cycle
   code on a ground-based star arrived as jitter and nothing else. The
   11,000 jitter-only SURVIVORs track the Blazhko incidence (ZTF RRab 21%,
   RRc 5%; ASAS-SN 9.2% vs 1.0%), scale with the number of merged ZTF object
   ids (8% → 21% → 27%) and fall with magnitude: a Blazhko meter.
3. **Slow codes.** The fold searched k ∈ {1, 2, 4, 8} cycles per bit. A code
   held for a night, a week or a month per bit was invisible.
4. **Low duty cycle.** `clip_steps` winsorized any step past 3σ when fewer
   than a tenth were out there. A transmitter that saves for N cycles and
   kicks once — the cheapest beacon there is — had its set bits flattened.
   The >5σ jumps were counted, vetted only for sitting at segment edges, and
   never read.
5. **Prewhitening** may remove a strong code's sidebands (untested).
6. **Stage 2's false-negative rate** is unmeasured on real hosts; the
   boundary-shift test damages exactly a one-bit-per-cycle staircase.

Items 3 and 4 and the sparse half of item 2 are what this note fixes.

## The kick route

Read the *raw* steps (before the clip). A kick is a step more than
`KICK_SIGMA` = 4 robust σ (IQR/1.349) from the median; the side with more of
them is the code's side, the other side counts against it, and the set is
one-sided when the counter-kicks are at most a quarter of the kicks (a bad
cycle is one step down and one back up). With `KICK_MIN` = 4 kicks the
spacing is read two ways:

- **Catalog gaps.** The longest run of consecutive gaps that is a window of
  an integer catalog (primes, counting, Fibonacci, squares, triangular,
  powers of two) in either encoding: kicks *at* the catalog's values as
  cycle numbers (gaps = the catalog's first differences; kicking on every
  prime-numbered cycle gives 1, 2, 2, 4, 2, 4, …), or gaps *equal* to the
  values (kick, 2 cycles, kick, 3, kick, 5, …). A gap between kicks with an
  unobserved cycle between them could hide a kick and ends any run.
  `KICK_MIN_RUN` = 3 to name anything.
- **Prime separations.** The fraction of kick pairs whose separation is
  prime — the Hippke et al. (2015) statistic, asked of the kicks rather than
  of the O–C autocorrelation.

Both are add-one Monte-Carlo p-values over `nshuffle` random sets of as many
kicks drawn from the same observed cycles (`_kick_null_nb`), so the number
of kicks, the gaps in the data, the two encodings and the catalogs' own
trials are all inside the p-value. Route: one-sided and (run ≥ 3 at
p < 1e-3, or prime fraction at p < 1e-3). The clip itself now leaves a
one-sided cluster of same-sized outliers alone (`one_sided_cluster`:
≥ 3 outliers, one side, spread ≤ half their median distance), so the fold
decoder reads a low-duty word as a word; the null loop clips every draw by
the same rule.

Why 4σ and 4 kicks: the Gaussian tail past 4σ is 6e-5 per step, so a
thousand-cycle series produces one chance kick in about one run in sixteen,
and four on one side by chance is negligible; on 12 quiet dense TESS RR Lyrae
the detector never found more than three.

## The epoch fold

`epoch_offsets` gives each epoch's implied timing offset (residual over
template slope) and `cycle_timings` the inverse-variance mean per *observed*
cycle — the series `sparse_timings` averages into blocks, kept at the cycle.
Almost no two adjacent cycles are ever observed, so the staircase cannot be
differenced; instead the detrended levels are folded on (L, k) and the bits
are the cyclic differences of the residue means (`_fold_level_word`,
`_level_fold_cost_tab_nb`). The cost is the description length of the levels
about the staircase the word predicts, residue j at a − s(S_j − Wj/L), fitted
by weighted least squares: n ln(RSS/n) + L ln 2 + 2 ln n, the same currency
as the step model. A repeated one-sided word walks down by dt per set bit
and the ramp of its duty cycle is what the quadratic detrend removed, which
is why the mean-subtracted cumulative sum is the right shape.

The null permutes the levels among the observed cycles and reruns the whole
(L, k) search on each draw. For everything downstream that was written for
steps — squareness against the best free sinusoid, the matched filter, the
per-cycle votes — `pseudo_steps` pairs each point with the nearest
observation of the previous bit's residue: its fold means are the cyclic
differences of the level means, its variance is twice the level noise, and
it carries the original cycle. Routes: the decode gate with
`epoch_seq_square`, or the sequence gate with `epoch_seq_window_square`. The
per-cycle series goes to `stage1_epoch_oc`, so `--recompute-routes` re-reads
the fold from the store.

Noise: a difference of two residue means has variance σ²(1/n_j + 1/n_{j−1}),
the same as folding differenced steps over the same counts, so nothing is
lost relative to the dense path except the √2 a step already pays; what is
gained is that the fold uses every observed cycle instead of the few
consecutive pairs.

A white null on a red clock is the same caveat the rest of stage 1 carries.
On real quiet and Blazhko stars the smooth modulation saturates the null in
a few per cent of stars, and the squareness guard is what keeps them out
(see below).

## The bit-period grid

`bit_cycles_for(n)` extends (1, 2, 4, 8) in powers of two while a 4-bit word
at that period still fits twice into the series. A 40-cycle series searches
what it did; a 3,000-cycle one goes to 256 cycles per bit. The decode nulls
redraw the grid, so the extra look-elsewhere is charged. `SEQ_CODE_VERSION`
is 4; stored words carry 3 and are decoded live by the viewer until
`--reseq` or a stage-1 pass rewrites them.

## Measured

Synthetic (tests/test_kick_epoch.py):

| series | result |
| --- | --- |
| 21-bit counting staircase, unit per-cycle S/N, 15–30% of 3,000 cycles | fold model reads it, names `counting [1..8]`, square, p at the null floor |
| the same with 40 units of noise on 30 of signal | still read and named |
| noise, same sampling | p ≈ 0.1–0.4, not square, unrecognized |
| 17 kicks with gaps 5…20 in 260 cycles (6.5% duty) | all 17 found, `counting as kick gaps [5, 6, 7, …]`, p at the floor, promoted |
| kicks at prime-numbered cycles | `kicks at primes [...]`, p at the floor |
| 12 one-sided kicks at random cycles | one-sided, no label, p ≈ 1 |
| two up, two down glitches | not one-sided; nothing read |
| an 8-bit word held 32 cycles per bit | decoded at k = 32 |

Real (read-only, cached light curves, 1,000 draws):

| sample | new-route promotions | notes |
| --- | ---: | --- |
| 24 quiet ZTF RR Lyrae | 0 | one saturated the epoch null, not square |
| 16 quiet ASAS-SN RR Lyrae | 0 | one saturated, not square |
| 24 ZTF jitter SURVIVORs (Blazhko-like) | 0 | three saturated, none square |
| 12 quiet dense TESS RR Lyrae | 0 | at most three one-sided kicks |

Cost: the epoch fold adds 1–4 s per sparse star at 1,000 draws (rows ≈ 260
(L, k) folds over 200–1,000 observed cycles); the kick null is negligible
because most stars have no kick set to test.

## Not done

- Stage 2 has no test aimed at a kick or an epoch-fold claim. A star promoted
  on either is judged by the existing battery, whose boundary-shift and
  richer-template tests were written for block and per-cycle O–C structure
  and whose `claims` block does not name these routes. The natural vets: the
  kick set should reproduce with the boundaries moved and the richer
  template; the epoch word should decode the same in each half of the
  baseline (the analogue of `fold_stability`).
- The viewer's list and score read the block-series word for a sparse star.
  The epoch word is in the payload (`epoch_code`, `epoch_seq_*`) and the
  detail's `timing` block, but no panel draws the per-cycle series yet, and
  `score.sequence_bonus` does not pay on `epoch_seq_p`.
- Injection has not measured either route's recovery on real hosts. The
  sparse injection machinery can now inject a kick word (`random:<n>@duty`)
  and the epoch fold will read it; that grid is the next measurement.
- The kick catalogs stop at 4,096 cycles per gap and do not include the
  binary expansions; a kick code spelling π in unary is not read.
