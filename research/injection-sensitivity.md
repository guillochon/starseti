# What amplitude of phase code would this search have found?

Measured 2026-09-16 with `starseti inject`, against the light curves in
`data/lc_cache`. Five grids: 180 real TESS hosts, 320 uninjected control runs,
1,980 injections. Every host is a star stage 1 already called quiet, and every
host was also re-run with nothing injected.

The injection is a relabelling of the time stamps, `t' = t − frac·P·bit(cycle)`,
so the flux, the noise, the cadence, the gaps and whatever astrophysics the host
carries are the star's own. Amplitudes are quoted three ways: `dt/P` (the
Learned et al. trigger scaling is 1e-3), `A/err` against the median per-cycle
timing error, and `A/rms` against the host's own O–C scatter, which is the noise
the detector actually faces.

## The headline

**RR Lyrae are the only population this search can hear, and it can hear them at
the right amplitude.** At the Learned scaling itself, dt/P = 1.2e-3, 42% of
injections into TESS RR Lyrae reach triage and 61% have their word decoded
correctly. At dt/P = 2e-3 the decoder is at 90%.

**The Cepheid arm — the one the project is named for — has no sensitivity to any
credible code.** The same dt/P = 1e-3 injection reaches triage 2% of the time.
It would take dt/P ≈ 0.05, a twentieth of a cycle and fifty times the formal
timing error, to be found at all, and by then it is not a subtle beacon.

## The hosts, and why they differ

| class | hosts | O–C points | P (d) | err (s) | rms (s) | rms/err | rms/P |
|---|---:|---:|---:|---:|---:|---:|---:|
| cepheid | 60 | 18 | 3.27 | 285 | 4156 | **24.9** | 0.024 |
| rrlyrae | 80 | 39 | 0.54 | 10 | 20 | **2.2** | 0.000 |
| bcep | 40 | 118 | 0.26 | 129 | 822 | **6.2** | 0.034 |

`rms/err` is the most important number in this table. The per-cycle timing
errors are already inflated by √χ²red inside `_dense_timings_nb`, so this is
genuine between-cycle wander: on a typical TESS Cepheid the O–C moves **25 times
further than the per-cycle precision allows**, and on β Cep 6 times. RR Lyrae
are nearly white (2.2), which is the whole reason they are the population that
still works.

Two consequences beyond this note:

1. `score._snr_factor` compares `oc_rms_seconds` against
   `median_timing_err_seconds`. On this evidence that ratio is measuring the
   red-noise excess, not a signal, and it is ~25 for a *quiet* Cepheid.
2. The stage-1 null is a random permutation, which is white by construction.
   Against a series whose own scatter is 25× its error bars that is the wrong
   null, and it is consistent with the gate firing on almost nothing.

## Completeness at the triage gate

An 8-bit random word, one bit per cycle. "decoded" is ≥ 90% of cycles recovered,
polarity-free; "chance" is the same decode against the *uninjected* run of the
same star.

| class | dt/P | A/err | A/rms | triage | decoded | chance |
|---|---:|---:|---:|---:|---:|---:|
| cepheid | 1.0e-3 | 1.2 | 0.03 | 2% | 2% | 0% |
| cepheid | 1.0e-2 | 11.7 | 0.3 | 6% | 23% | 0% |
| cepheid | 3.0e-2 | 35.0 | 1.0 | 8% | 31% | 0% |
| cepheid | 5.1e-2 | 50.8 | 2.0 | 7% | 40% | 0% |
| rrlyrae | 3.9e-4 | 2.2 | 1.0 | 2% | 8% | 0% |
| rrlyrae | 7.9e-4 | 4.4 | 2.0 | 25% | 29% | 0% |
| rrlyrae | 1.2e-3 | 6.7 | 3.0 | **42%** | 61% | 0% |
| rrlyrae | 2.0e-3 | 11.1 | 5.0 | **57%** | **90%** | 0% |
| bcep | 2.7e-2 | 5.6 | 1.0 | 21% | 15% | 0% |
| bcep | 3.7e-2 | 8.2 | 1.5 | 50% | 32% | 0% |
| bcep | 4.4e-2 | 9.1 | 2.0 | 60% | 56% | 0% |

50% recovery: RR Lyrae at 4.0× the host's own O–C rms (dt/P ≈ 1.6e-3), β Cep at
1.5× (dt/P ≈ 3.7e-2), Cepheids nowhere inside the grid. Nothing reaches 90%.

The chance column is measured, not assumed: **0 of 1,631 uninjected decodes**
reached 90% recall against the word that would have been injected, the highest
being 0.884. And 0 of the 320 uninjected runs reached triage, so the
re-analysis reproduces the store's own "quiet" verdicts exactly.

## Three findings that are not about amplitude

### 1. The gate is blind to long words. The decoder is not.

Same hosts, same amplitudes, two word lengths. `counting` is the 21-bit catalog
word; `random:8` is 8 bits.

| word | dt/P | triage | decoded | catalog match |
|---|---:|---:|---:|---:|
| random:8 (8 bits) | 2.0e-3 | **57%** | 90% | — |
| counting (21 bits) | 3.3e-3 | **5%** | 57% | 52% |

A longer word at a *higher* amplitude is eleven times harder for triage. The
reason is structural: Ljung–Box and max|ACF| are powered by *repetition*, and an
8-bit word repeats five times in a 39-cycle series while a 21-bit word repeats
under twice. The gate is tuned to short, highly redundant codes — and redundancy
is exactly what a message carrying information does not have.

The decoder does not care. It recovers the 21-bit word from 57% of those
injections, and `recognize_sequence` names it as the counting sequence in 52%.

### 2. The pipeline decodes the word and then throws it away

Of the 168 `counting` injections, **34 produced an exact catalog match at
seq_p = 0.002** — a 21-bit word, `bit_recall` 1.0, four orders of magnitude
inside the `SEQ_P_GATE` of 0.01. Only 6 of those 34 reached triage.

The other **28 are classified QUIET**. `is_followup` promotes on Ljung–Box,
max|ACF| or jitter only, and `score.sequence_bonus` "does not mint a
classification", so a star whose O–C decodes exactly to a recognised numeric
sequence cannot reach stage 2, cannot be vetted, and never appears as a
candidate. On this measurement that is the single largest loss of sensitivity in
the pipeline, and it is a routing fix rather than a physics one:

> a decode-quality promotion path — `seq_p` below the gate, or high bit
> confidence plus squareness, scored against the measured 0-in-1,631 null decode
> rate — would turn 28 discarded exact matches into 28 candidates.

The same gap shows on Cepheids without any catalog word: at dt/P = 0.05, decode
40% against triage 7%. It closes only where the gate works well — β Cep at 2×
rms, triage 60% against decode 56% — that is, where there are 118 timed cycles
and a near-white O–C.

### 3. Recovery is set by how many cycles were timed

Within one class at one amplitude, so that only the length of the O–C series
varies — TESS RR Lyrae, 8-bit word, at 3× and 5× the host's own O–C rms:

| O–C points | hosts | triage @3× | decoded @3× | triage @5× | decoded @5× |
|---|---:|---:|---:|---:|---:|
| < 25 | 9 | 11% | 33% | 11% | 89% |
| 25–40 | 32 | 19% | 41% | 44% | 81% |
| 40–60 | 29 | 59% | 83% | 72% | 97% |
| ≥ 60 | 10 | **100%** | 90% | **100%** | 100% |

The same monotonic climb appears at 2× (11 → 6 → 31 → 80%) and on β Cep. A star
observed for 60 pulsation cycles is found every time; the same star observed for
20 is missed nine times in ten at the same amplitude. Note also that at 5× the
decoder is already at 89% in the shortest bin while triage is at 11% — the gap
in finding 1 is widest exactly where the data is thinnest.

The median cached TESS Cepheid host has 18 timed cycles. That is a coverage
problem — more sectors, longer baselines, cross-survey stitching — not a
statistics problem, and no threshold change will fix it.

It also bounds what can be searched for: 307 of 1,980 injections could not run
at all because the word needs more cycles than the host was observed for (an
8-bit word needs 16; the catalog words need 40–56, which the median Cepheid
misses by a factor of three). A further 42 were refused for the opposite reason
— the amplitude needed to beat that host's O–C scatter exceeded a tenth of a
cycle, where the injection starts dragging `find_period` with it. Both are
recorded as out of range rather than as a miss.

## A note on cycles per bit

`decode_bit_code` folds on `cycle % L`, so one bit per pulsation cycle is baked
into it. A code at 4 cycles per bit returns as a word of 4× the length with each
bit repeated — still a recovery per cycle, but no longer a catalog word.
Measured on Cepheids at dt/P = 1e-2: 23% decoded at one cycle per bit against
14% at four. A search over bit rate needs the fold to search it too.

## After the decode route was added

`is_followup` gained the two routes finding 2 argues for: `p_decode_shuffle`
(the description-length statistic, calibrated against the same shuffle null as
the ACF) together with `seq_square`, and a catalog match at `seq_p < 0.01` on
its own. The grids were re-run unchanged on the same hosts.

**The failure this was built to close is closed.**

| | before | after |
|---|---:|---:|
| exact catalog matches (21-bit word, 168 injections) | 34 | 34 |
| ... of those promoted to stage 2 | **6** | **34** |
| ... of those discarded as QUIET | **28** | **0** |
| correctly decoded words promoted (all RR Lyrae cells) | 43% | **77%** |
| false positives on 80 uninjected controls | 0 | **0** |

Recovery at the triage gate, same hosts and amplitudes:

| class | word | A/rms | before | after |
|---|---|---:|---:|---:|
| rrlyrae | random:8 | 2.0 | 25% | **38%** |
| rrlyrae | random:8 | 3.0 | 42% | **60%** |
| rrlyrae | random:8 | 5.0 | 57% | **80%** |
| rrlyrae | counting (21 bits) | 2.0 | 2% | **12%** |
| rrlyrae | counting (21 bits) | 3.0 | 4% | **20%** |
| rrlyrae | counting (21 bits) | 5.0 | 5% | **45%** |
| cepheid | random:8 | 1.0 | 4% | **10%** |
| cepheid | random:8 | 1.5 | 4% | **11%** |
| cepheid | random:8 | 2.0 | 7% | **19%** |

The long-word case — the one a message actually looks like — gains a factor of
five to nine. The 8-bit case gains about a half, since autocorrelation was
already finding some of those. Cepheids gain a factor of three and remain poor
in absolute terms, which is the right outcome: the new route addresses the gate,
and nothing about the gate fixes an O–C whose own scatter is 25× its error bars
on a series 18 cycles long. Findings 2 and 3 are still open.

Two things the re-run also showed. The squareness side condition is load-bearing
rather than decorative: an injected *sinusoidal* O–C saturates
`p_decode_shuffle` — the shuffle null is white, so a smooth wave beats it — and
is kept out of the decode route by `seq_square` alone (it still reaches triage
through max|ACF|, as it should, and stage 2 vetoes it as astrophysics). And on
Cepheids only 31% of correctly decoded words are promoted against RR Lyrae's
77%, because red noise makes their two levels ragged enough to fail the square
test; that is the honest ceiling of this fix on that population.

## Searching for the amplitude instead of sampling it

*Added 2026-09-17.*

Everything above is a fixed grid: five `--snr` values, and a completeness limit
read off by interpolating between two adjacent points. That is a measurement
only when the grid happens to bracket the answer. For TESS Cepheids it never
did — grid G ran out to dt/P = 5e-2 and peaked at 19% triage — so the number
the whole stage exists to produce came back as `--`, and a null result with no
threshold attached is not a null result.

`inject --find-threshold` searches for the amplitude:

1. **Bracket.** A log-spaced scan over dt/P from 1e-5 to `--max-frac`, run as
   cells of a single `run_host` call per star, so the light curve is loaded,
   prepared and null-analysed once for the whole pass.
2. **Refine.** A second pass placed inside whichever interval the first one
   crossed the target rate in. Two passes over the hosts, not the eight a
   bisection would cost.
3. **Fit.** A two-parameter logistic in log amplitude,
   \(p(x) = [1 + e^{-(a + bx)}]^{-1}\), fitted by maximum likelihood to every
   injection from both passes, and the threshold read off it as
   \(x = (\mathrm{logit}(\mathrm{target}) - a)/b\).

The fit is the part that matters. Recovery at one amplitude over \(n\) hosts is
a binomial with \(\sigma \approx 0.5/\sqrt{n}\) — ±8% at \(n = 40\) — and
interpolating between two such points propagates all of that into the
threshold. Using every point instead is what the monotone shape of a
completeness curve is worth. A fit whose slope comes out non-positive is
refused rather than extrapolated: recovery that falls with amplitude is a
broken injection, not a threshold.

**Two gates, because they answer different questions.** `triage` is what
`is_followup` would have flagged; `decoded` is what the pipeline could actually
*read* (bit recall ≥ 0.9 against the injected word). The gap between them is
the interesting part, and it does not always point the way one would guess —
see below.

**The report prints two 50% numbers and they are not the same statistic.** The
"Amplitude required" table is the fitted threshold; the sentence under each
completeness curve is still the old linear interpolation between adjacent grid
points, which is worth keeping because it is a raw reading of the data with no
model in it. On the TESS RR Lyrae search they came out at 2.5e-3 and 1.7e-3:
the curve reads 50% at 1.7e-3 and then sits on a 50–52% plateau out to 2.5e-3
before jumping to 88%, so the interpolation takes the first crossing of a noisy
plateau and the fit takes the whole shape. The fit is the one to quote — 40
hosts per point is ±8% on each — but a large gap between the two is a sign the
curve is not logistic and is worth looking at rather than averaging away.

**Bounds are reported as bounds.** Above dt/P ≈ 0.1 the injection stops being a
phase code: the shift is a sizeable fraction of a cycle, `find_period` and the
template follow it, and the row measures the damage rather than the code. A
class that has not reached the target rate by then is recorded with
`is_bound`, and everything downstream — the report, the store, the viewer's
Kardashev class — carries the `≥`. The search establishes that the requirement
is *at least* that much and nothing more.

### Measured 2026-09-17

TESS, `random:8` at one bit per cycle, 40 hosts per class, scan 6 + refine 4.

| class | gate | dt/P for 50% | best rate seen | injections |
|---|---|---:|---:|---:|
| rrlyrae | triage | 2.5e-3 | 90% | 400 |
| rrlyrae | decoded | 1.7e-3 | 100% | 400 |
| cepheid | triage | **≥ 1.0e-1** | 53% | 320 |
| cepheid | decoded | 4.8e-2 | 62% | 320 |

(β Cep was not measured — the run was stopped part way, 161 rows under
`H-bcep-threshold` and no threshold written. That class still falls back to the
Learned scaling in the viewer.)

Two things in this table are worth more than the numbers.

**The decode gate needs *less* amplitude than triage, in both populations.**
RR Lyrae: 1.7e-3 to read a word against 2.5e-3 to be flagged. Cepheids:
4.8e-2 against a triage threshold the search cannot even establish. The
pipeline can decode a word at amplitudes where it would never have looked at
the star — which is the same gap commit 866788c opened the decode route to
close, still open at the 50% level. `is_followup`, not the decoder, is what
sets this search's sensitivity.

**The Cepheid triage bound is a plateau, not a cliff.** The curve goes 12% at
2.5e-3, 28% at 1.6e-2, then sits at 34% from 2.3e-2 to 6.9e-2 before reading
53% at the 0.1 ceiling. So the fitted crossing is at or past the ceiling and
the honest report is "≥ 1.0e-1", but the shape is a detector saturating
against a red-noise floor rather than one that simply needs more amplitude —
and the last point is one 32-host sample, ±9%. The `≥` is doing real work
here: the search establishes the requirement is at least a tenth of a cycle,
and the mechanism above says raising the amplitude further would not fix it.

Also: 80 of 400 Cepheid injections were dropped as `too-short` — an 8-bit word
at one bit per cycle needs 16 cycles to repeat twice, and a median TESS Cepheid
carries 18. The Cepheid numbers rest on 320 injections over the hosts long
enough to test, not 400.

## From a threshold to a transmitter power

The viewer's energy-per-bit box and Sagan–Kardashev class used to assume the
Learned et al. trigger scaling, dt/P = 1e-3, for every star:

\[
E_{\mathrm{bit}} = (\Delta t/P) 	imes c 	imes L 	imes P_{\mathrm{s}}
\]

with the whole right-hand side collapsed into one constant. That is a
statement about the *hypothesis*, not about this search: it says what a
transmitter would spend, not what it would have to spend for us to notice.
With `--publish` the measured dt/P for the star's own survey and class is used
instead, so the number on the page is the power this pipeline's sensitivity
demands. A class with no published threshold falls back to the Learned scaling
and the tooltip says so, so the two are never confused.

**One discrepancy surfaced and was left alone.** The coupling constant \(c\) in
the viewer is 0.1, while the tooltip beside it has always read "×10 for
acoustic coupling" — a factor of 100, or 0.2 in Sagan K. `TICKLE_COUPLING` now
names it instead of burying it in a magic constant, and keeps the value the
viewer has been showing, so this change restates no published figure. Which of
the two the model intended is a question for whoever wrote it down.

## The waveform was wrong, and the right one is harder

*Added 2026-09-17, branch `staircase-physics-and-recovery`.*

Everything above injected a **square wave**: a set bit shifted its own cycle
by dt and the next clear bit sat back on the original ephemeris. That is not
what a star does. A pulsator is a limit-cycle oscillator with no restoring
force along the phase direction, so a triggered pulse resets the clock and
every later cycle inherits the advance — the O–C of a phase code is the
running sum of its bits, a one-sided **staircase**, and the bits are its
*steps* ([modulation-physics.md](modulation-physics.md), §2). The injector,
the built-in controls and the pinned viewer fixtures now write the cumulative
shift, `decode_bit_code` thresholds the first difference of the O–C by
default, and stage 1's shuffle null is differenced the same way. The square
wave survives as `--waveform square` / `model="square"` for comparison.

The note predicted the staircase would be ~5× *easier* to detect, on the
grounds that the demeaned running sum of an L-bit word has amplitude ~L/4 dt.
The same hosts, re-measured, say the opposite:

| class | gate | square wave (before) | staircase | ratio |
|---|---|---:|---:|---:|
| rrlyrae | triage | 2.5e-3 | **5.5e-3** | 2.2× |
| rrlyrae | decoded | 1.7e-3 | **6.3e-3** | 3.7× |
| cepheid | triage | ≥ 1.0e-1 (best 53%) | **≥ 5.0e-2** (best 31%) | bound |
| cepheid | decoded | 4.8e-2 | **≥ 5.0e-2** (best 19%) | bound |
| bcep | triage | 3.7e-2 (grid D, interpolated) | **≥ 5.0e-2** (best 52%, at the ceiling) | bound |
| bcep | decoded | not measured | **≥ 5.0e-2** (best 20%) | bound |

TESS, `random:8` at one bit per cycle, 40 hosts per class, 360 (RR Lyrae), 160
(Cepheid) and 360 (β Cep) staircase injections; run ids `S-<class>-staircase`. The
Cepheid search never bracketed the target, so both gates are bounds at the
staircase ceiling, which is 0.05 rather than the square wave's 0.1 — the
detrended running sum of an 8-bit word swings ±2 dt where a square wave
swings ±dt/2, so the amplitude at which the injection starts smearing the
template is half as large. β Cep shows what that ceiling means: at 0.05 the
triage rate jumps to 52% while the decode rate falls to 0% and the O–C
statistics run for only 80% of the hosts — the template is following the
code, which is exactly the row the ceiling exists to exclude. Its threshold
is a bound, and the 52% is not a measurement of the code.

**Why the prediction failed.** The L/4 excursion belongs to a word with long
runs, `1111 0000`. A near-fair-coin word — which is what a message carrying
information looks like — keeps its running sum close to its own ramp: for the
injected `random:8` word the detrended staircase has an rms of 0.34 dt
against the square wave's 0.48 dt, so the *level* series the autocorrelation
routes read carries less periodic power, not more. And the bit itself now
lives in the difference of two noisy timings, √2 noisier than a level, which
is the decode gate's 2–4×. Only a low-entropy message gets the gain the note
described. The physical waveform is harder to find than the one this search
had been measuring, and every square-wave figure above was optimistic by
2–4×.

**Two things the diagnosis turned up that are not about the waveform.**

- *The per-cycle timing compresses the shift.* On a real TESS RR Lyrae host
  (TIC 152231997, 30 timed cycles, 17 s formal error) both waveforms come
  back at 0.6–0.7 of the injected amplitude from dt/P = 2.5e-3 to 2e-2,
  while the synthetic star in `tests/test_inject.py` returns 0.9–1.0. The
  thresholds are end-to-end and already include it, but a 35% loss of signal
  on the sky that is absent from the idealisation is unexplained and worth a
  note of its own.
- *Differencing whitens the noise, and the squareness test notices.*
  `seq_squareness ≥ 0.25` now holds for about twice as many real stars (424
  of 5,581 vetted ASAS-SN rows against 231), because the sinusoid that used
  to soak up red noise in the level series explains nothing in its
  increments. The predicate still separates a two-level word from Blazhko,
  and in stage 1 it only ever promotes in conjunction with `p_decode_shuffle`,
  but its base rate has moved and stage 2's `_ltte_dominates` reads it; the
  566 stored stage-2 rows that either decoder called square were re-vetted.

**What it does to the power figure.** With the coupling derived rather than
assumed (×1.67 for an annular beam at unit gain in place of the old 0.1,
+0.12 in Sagan K for every star) and the RR Lyrae threshold at 5.5e-3, a
transmitter at the reference RRab spends P_beam/L ≈ 4.6e-3 — 4.6% of a
Dyson sphere at η = 0.1, 46% at η = 0.01, and impossible at η = 1e-3. The
Cepheid bound at 0.05 is P_beam/L ≥ 0.042: 42% of the star's sky even at
η = 0.1. Both are within what a self-powered annular transmitter can supply
at the optimistic efficiency, and the report's "Power the threshold implies"
section carries the same table per run.

**The concrete next step** is a triage route that reads the increments, not
only the levels: for a staircase of a fair-coin word the differenced O–C is
exactly the square wave the autocorrelation routes were built for, at √2 the
noise, and a Ljung–Box or max|ACF| on it would give back most of what the
level statistics lost. It needs its own null and its own look-elsewhere
accounting, so it is not in this branch.

## Listening for a designed beacon

*Added 2026-09-17, later the same day.* The staircase result above says a
message-carrying word is harder to find than the search assumed. A beacon
*designed* to be found, under the two constraints a transmitter actually has —
distinguishability from anything a star does, and no place to store cycles
of energy for a big kick — is a long recognisable word, repeated, one bit
per cycle, at a dt/P a self-powered transmitter can kick every cycle:
`η f_coll / coupling` ≈ 3e-4 at η = 0.01 and 5% of the star's sky
(`physics.kick_ceiling`). That is an order of magnitude below every
threshold measured above, and integration over repeats is the only way down
to it. Three things were changed, and one bug found on the way.

**The timing was clipping the signal.** `prepare` fits the template with a
5σ clip against the mean ephemeris and `dense_timings` was fed only the kept
points, so a cycle that genuinely arrived early lost the steep-branch points
that carry its timing: on TIC 152231997 both waveforms came back at 0.55–0.72
of the injected shift. The timers now take every point and clip per cycle
against the cycle's own fine-grid shift (never more than half); recovery is
0.83–1.00 from dt/P = 2.5e-3 up. The same clip had been suppressing genuine
cycle-to-cycle wander, so the noise floor was understated too (that host's
null rms went from 20 s to 44 s).

**Fold first.** The decode statistic and `decode_bit_code` average the
steps per residue of `(cycle // k) % L` for k in 1, 2, 4, 8 and every L, then
split the *means* (Otsu, so a lopsided word is no longer undecodable), pricing
the word at ln 2 per bit and ln n per level. Two new triage routes read the
steps directly: **quantised steps** (Otsu's between-class share of the step
variance, 2/π for Gaussian steps, → 1 for two levels) and a **matched catalog
filter** (`matched.py`) that correlates the folded steps against every window
of the packed numeric catalogs at the fold's own length, divides by the noise
of the whole series, charges every window and rotation, and requires the
residuals to fit the noise and the window to fit as well as the fold's own
free word. Both routes, like the decode route, only promote when the two
levels beat a sinusoid: without that guard they fired on 7% and 9% of real
re-timed RR Lyrae, whose smooth modulation has an arcsine step distribution
and resembles square words at ρ ≈ 0.9. Stage 2 decodes each half of the
baseline and asks for the same word.

**Measured** (run ids `W-<class>-fold`; TESS, 40 quiet hosts per class drawn
from the *re-timed* store, staircase, one bit per cycle; 0 of 240 uninjected
runs reached triage):

| class | word | gate | before (staircase, clipped timer) | now | repeats seen |
|---|---|---|---:|---:|---:|
| rrlyrae | random:8 | triage | 5.5e-3 | **4.0e-3** | 8 (1.5–52) |
| rrlyrae | random:8 | decoded | 6.3e-3 | **4.2e-3** | |
| rrlyrae | random:8 | quantised | — | 2.6e-2 | |
| rrlyrae | random:8 | catalog | — | never (correct: it is not a catalog word) | |
| rrlyrae | primes (28 bits) | triage | — | **2.4e-3** | 2.3 (0.6–15) |
| rrlyrae | primes | decoded / catalog / quantised | — | ≥ 5e-2 (best 36% / 10% / 10%) | |
| cepheid | either | every gate | ≥ 5e-2 | ≥ 5e-2 (best 27% triage) | 3.0 / 1.0 |
| bcep | either | every gate | ≥ 5e-2 | ≥ 5e-2 (triage best 47–48% at the ceiling; decode 0%) | 40 hosts drawn from the 87 cached β Cep the re-timed store still calls quiet |

The comparison is not clean: the hosts are drawn from the store after the
cached TESS stars were re-timed, and the unclipped timer moved 274 of 465
cached RR Lyrae into triage, so the quiet pool the search draws from is a
different, quieter set. Within that caveat the random word is ~1.4× easier at
both gates and the long word reaches triage at 2.4e-3 — the staircase of a
word with long runs has the excursion the autocorrelation routes want — but
cannot be *read* or *named* on these hosts, which hold 2.3 repeats of it on
average and under one at the short end. Against the kick ceiling the gap is
now 8× (primes, triage) to 13× (random word); it was 18× before the fold.

**What the false positives taught.** The matched filter went through three
versions in an afternoon, each caught by an injection run rather than a unit
test. A t-test on the correlation called a noise-free differenced square wave
a certain match at ρ = 0.8. Folding at every length let a strong random word
alias into structure whose within-residue noise estimate was inflated by the
signal itself, and 65% of RR Lyrae injections "matched" the counting sequence.
Taking the noise from the fold that was *chosen* for compressing well
understated it, and pure noise fired six times in a hundred. The version that
stands: fold at the decoder's length only, noise from the whole series,
goodness of fit against that noise, at least 97% of the free word's fit, and
the squareness guard. Pure noise and weak random words: 0 of 200 each. It is
worth recording because the same reasoning applies to any receiver that names
a word: a p-value that does not include how the length was chosen, or a noise
estimate taken from the thing being tested, is not a p-value.

**Store.** Every stored O–C statistic was recomputed from the stored series
(`stage1 --recompute-routes --triage`, no light curves) so the new routes and
the corrected decode reach every star; the 1,375 TESS stars with cached light
curves were fully re-timed as well. Triage went 13,999 → 14,103; 66 stars
are in triage on the new routes alone (50 on the fold-first decode, 8 on
quantised steps, 3 on an exact catalog word, the rest on combinations), and
stage 2 vetoed every one of them (58 artifact, 7 red noise, 1 astrophysics).
That is the null the routes were built to make meaningful: at the amplitudes
the search now reaches, no star in the store carries a repeated two-level
word. The ~16,000 other TESS rows and every sparse row still carry O–C
compressed by the old clip and need a fetch.

## What this does not cover

- **Sparse surveys.** Every host here is TESS dense mode, and the block-averaged
  and jitter channels are untested even though two thirds of the sample lives
  there. The reason was mechanical rather than scientific: ASAS-SN and ZTF had
  no light-curve cache, `hosts_from_store` draws cached hosts by globbing FITS
  off disk, and re-downloading ~25 s of photometry per injection made a grid
  impossible. `lc_cache.py` now persists both surveys' raw rows and
  `cache_index` lists them, so `--survey asassn` / `--survey ztf` can draw
  hosts the same way TESS does. **The grid has not been run.** It needs a
  prefetch pass over the hosts first, and nothing below has been measured on a
  sparse light curve.
- **Stage 2.** `--stage2` vets every injection that reaches triage, but these
  grids ran without it. The numbers above stop at the stage-1 gate, and the
  vetting battery can only lose candidates from there.
- **Hosts that are not quiet.** Hosts are drawn from stars stage 1 already
  called quiet, so the background here is a reproducibility check and not a
  false-positive rate. `--any-host` measures that one.

## The sparse arm, measured (2026-09-24)

*Runs `R1-asassn-rrlyrae`, `R2-ztf-rrlyrae`, `R3-tess-rrlyrae-anyhost` and
`R4-tess-rrlyrae-beat`, made for the response to review 1 (M2, M3 and the beat
control promised in Table 3). 40 quiet hosts each; the sparse runs at 1, 8 and
32 cycles per bit with the kick and epoch-fold routes of
[kick-and-epoch-routes.md](kick-and-epoch-routes.md).*

**A one-bit-per-cycle code is invisible on a ground-based star** at any
amplitude that is still a phase code: the flag rate never exceeds 18% below
dt/P = 0.05 on either survey, for either word. Block averaging sums tens of
bits into every step, as predicted.

**A slow code is heard and not read.** Flagged-gate 50% thresholds (logistic
fits excluding the 0.05 ceiling row, where the template follows the code and
the rate collapses to 0-18%):

| survey | word | 8 cycles/bit | 32 cycles/bit |
|---|---|---:|---:|
| ASAS-SN | random:8 | >= 6e-3 (best 33%) | 2.6e-3 |
| ASAS-SN | primes | 4.4e-3 | 8.9e-4 |
| ZTF | random:8 | 5.8e-3 | 1.5e-3 |
| ZTF | primes | 2.2e-3 | 6.3e-4 |

Not one of 4,280 sparse injections at any bit period was decoded (bit recall
>= 0.9) or named by the catalog matcher at any gate. The fitted thresholds the
report itself prints for the k = 32 random word (3.5e-2 ASAS-SN, and 2.7e-2
for the ZTF primes word) are artefacts of the ceiling row bending the
logistic; the raw interpolations (2.5e-3, 5.7e-4) agree with the fits above.
`inject --find-threshold` should exclude the ceiling row from the fit by
default; until it does, read the raw crossing beside the fit.

**Completeness is not much worse on hosts that are not quiet.** `--any-host`
on 40 TESS RR Lyrae, of which 24 (60%) reach triage with nothing injected
(the cached pool is the re-timed one, which is why this exceeds the 21%
population flag rate): at dt/P = 2.4e-3 the 8-bit word decodes in 65% and the
primes word in 31% (41% at 4e-3), against 50% at 4.2e-3 on quiet hosts. The
flag gate is uninformative here because most hosts were flagged already.

**The beat control does what it should.** A second mode at nu = 0.618 or
0.707 written into the clock of 40 quiet TESS RR Lyrae at 3, 5 and 10x the
O-C rms reaches triage in 74-94% of 228 injections (through the whiteness
routes, as a non-white clock must), passes the squareness guard in 0, is
promoted by a decode route in 0, decodes as a word in 0, and is named by the
matcher at some p in 8% (counting and primes windows deep in the 8-bit
packings), none below the 0.01 gate. The decode routes' false-promotion rate
on a beating star is therefore measured at 0 of 228.

**The flag rate is the false-positive rate.** The uninjected controls above are
a reproducibility check on hosts pre-selected as quiet. The empirical
false-positive rate of stage 1 is the flag rate of the survey itself, 39%
overall (Table 2 of the paper): TESS 7 / 21 / 55% (Cepheid / RR Lyrae / beta
Cep), ASAS-SN 14 / 46 / 12%, ZTF 24 / 44%.
