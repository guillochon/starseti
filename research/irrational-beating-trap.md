# The irrational-beating trap, and how not to fall into it

*Written 2026-09-18.* The one previous search for a phase-coded pulsator
(Hippke et al. 2015, ApJ 798, 42) found that the cycle-to-cycle period sequence
of the RRc star KIC 5520878 was autocorrelated at prime-number lags at nominal
99.8% confidence, and then found the cause: the star pulsates in two modes whose
period ratio is close to 1/√2 (Hippke, Learned & Zee 2014, arXiv:1407.3442), and
sampling an incommensurate second mode once per cycle of the first produces a
quasi-periodic sequence whose autocorrelation is structured at the denominators
of the ratio's continued-fraction convergents. Primes among those lags were a
coincidence; the "code" was arithmetic. This note asks whether the present
pipeline would make the same mistake. **It would**, on one route, and the fix
is cheap.

## 1. What the beat looks like to us

An O–C series sampled once per pulsation cycle carries a second mode as

    x_c = A sin(2π ν c + φ),      ν = P₁/P₂ mod 1

with ν irrational. For KIC 5520878, P₂/P₁ ≈ 0.6177, so ν ≈ 0.618 ≡ 1/φ (the
golden ratio), whose convergents are the Fibonacci ratios 1/2, 2/3, 3/5, 5/8,
8/13, 13/21, 21/34. For a ratio near √2 the convergents are 1, 3/2, 7/5, 17/12,
41/29. Folding x_c on any convergent denominator q almost repeats: ν q is within
1/q of an integer, so every residue collects samples of nearly the same phase
and the fold is a clean sinusoid over q residues. **A two-level split of that
fold is a stable q-bit "word".**

Measured with `sequence.decode_bit_code` and `stage1_pulse.analyse_oc` on a
pure two-mode series, 120 cycles, no code injected (script at the end):

| series | word length L | squareness (current) | catalog match | Ljung–Box | max\|ACF\| | prime lags | decode MDL | quantised steps |
|---|---:|---:|---|---:|---:|---:|---:|---:|
| ν = 0.6177, A = 3σ | 21 | **0.76** | none | 1e-3 (sat.) | 1e-3 | 1e-3 | 1e-3 | 1e-3 |
| ν = 0.6177, A = 10σ | 21 | **0.79** | none | 1e-3 | 1e-3 | 1e-3 | 1e-3 | 1e-3 |
| ν = √2, A = 6σ | 29 | **0.78** | none | 1e-3 | 1e-3 | 1e-3 | 1e-3 | 1e-3 |

Every whiteness route saturates, which is correct (the series is coherent), and
stage 2 is supposed to sort that out. The problem is the **decode route**:
`p_decode_shuffle` saturates, `seq_square` is *True* at 0.76–0.79 against a
gate of 0.25, and `p_quantised_shuffle` saturates too. A star doing nothing but
beating in two modes therefore reaches triage as a *candidate phase code*, the
same thing an injected staircase does. The decoded 21-bit word for the golden
ratio is

    1 1 0 1 0 0 1 0 1 1 0 1 0 1 1 0 1 0 0 1 0

which is not in any catalog, so the matched filter stays quiet
(`p_catalog_matched = 1`) and no star is *named* this way. But the promotion
to vetting, the "coherent two-level phase code" verdict and the QUANTISED and
DECODE routes are all fooled, and stage 2's `classify_secondary` would label
the 2.6-cycle secondary period only as a soft "possible unresolved beat" for an
RR Lyrae star (a hard ASTRO veto needs the mode to have been prewhitened, and
RR Lyrae are not prewhitened by default). Concretely, KIC 5520878 would sit in
the WEAK or SURVIVOR class with a two-level "code" until a person looked.

## 2. Why the squareness test lets it through

`sequence._squareness` compares the two-level word's residual against a
sinusoid **at the word's own period L**. The beat is not a sinusoid at period
L = 21; it is a sinusoid at period 1/ν = 1.62 cycles, aliased to 1/(1−ν) = 2.62
cycles, which merely *folds* well at 21. A single sinusoid with period 21 fits
that fold poorly, the two levels fit it better, and the statistic reads "square".
The test answers the wrong question: not "is there a smooth wave that explains
this as well as the bits do" but "is there a smooth wave *at L*".

## 3. The fix: compare against the best sinusoid at any frequency

Replace the period-L sinusoid by the best-fitting sinusoid over all frequencies
0 < f ≤ 0.5 cycle⁻¹ (the Nyquist range for a once-per-cycle series), i.e. take
the periodogram peak of the *steps* and fit a + b cos + d sin there. The same
series, same decoded words:

| series | squareness vs sinusoid at L | **squareness vs best free sinusoid** | best f (cycle⁻¹) | \|fL − round(fL)\| |
|---|---:|---:|---:|---:|
| beat ν = 0.6177, 3σ | 0.76 | **0.00** | 0.3809 (= 1 − 0.6191) | 0.001 |
| beat ν = 0.6177, 10σ | 0.79 | **0.00** | 0.3812 | 0.005 |
| beat ν = √2, 6σ | 0.78 | **0.00** | 0.4143 (= √2 − 1) | 0.015 |
| **injected staircase**, synthetic RRab, counting word, 106 cycles | 0.94 | **0.93** | 0.0947 | 0.012 |

The beat's squareness collapses to zero: one sinusoid at the aliased beat
frequency explains everything the 21 levels explained. The real code keeps its
squareness, because no single sinusoid can reproduce a staircase's two-level
increments. This is a one-line change in what `_squareness` compares against,
and it costs a periodogram of a ≤480-point series.

Two remarks. The rightmost column, the *convergent test* (how close f·L is to an
integer), does flag the beats, but it flags the real code as well: any word that
folds cleanly at L has power near k/L by construction. It is a description of
why the fold worked, not a discriminator, and should not be used as a veto on
its own. And the recovered frequency is a physical measurement: 1 − 0.3809 =
0.6191 and 1 − 0.4143 = 0.5857 = 1/(1+√2−1)... reproduce the injected ratios to
three figures, so the same fit tells stage 2 *which* second mode to look for in
the light curve.

## 4. Two further discriminants worth carrying

**The word itself is a rotation coding.** Thresholding a sampled irrational
rotation gives the coding of `c → νc mod 1` by an interval: a *rotation word*.
(A first draft of this note called it Sturmian, i.e. balanced. That is only
true when the interval length equals ν; the decoder's two-level split sits
near the sinusoid's midpoint, so the density is 1/2, not ν, and the word above
contains both `00` and `11`, which a balanced word cannot. It is a Rote word,
not a Sturmian one.) What every rotation word does satisfy is a *factor
complexity* bound: at most 2m distinct substrings of each length m, against
the 2^m a message can use. Measured over cyclic substrings: the decoded words
of the golden-ratio beat (21 and 34 bits) and of the √2 beat (29 bits) all
satisfy the bound; the primes, Fibonacci and counting words all break it by
m = 3. A fair-coin word satisfies it with probability 80% at 8 bits, 8% at 16,
0.9% at 21 and 0.2% at 24 (20,000 draws each), so it is a usable label for a
long word and nothing more for a short one. A single noise-flipped bit can
break the bound (it did for the 34-bit word and not for the 21-bit one), so the
flag is a secondary label behind the amplitude test, not a replacement for it.

**Phase modulation puts sidebands at every harmonic, with amplitude ∝ k.** A
second pulsation mode is one independent peak in the light-curve spectrum, plus
combination frequencies at amplitudes that fall with order. A code that advances
the pulsation phase by ε(t) multiplies each harmonic k of the template by
exp(i k ε(t)), so the modulation sidebands around harmonic k have amplitude
k·A_k·ε: the *ratio* of sideband to carrier grows linearly with harmonic number.
That is a test on the light curve, not on the O–C, and no beat between two
radial modes reproduces it. It is the only diagnostic in this note that would
distinguish a genuine phase code from *any* additive second signal, and it is
not implemented anywhere in the pipeline yet.

**Prewhiten before timing, for every class.** β Cephei stars are prewhitened of
independent modes before the per-cycle fit precisely because their beats read as
O–C structure. RR Lyrae and Cepheids are not, on the argument that additional
modes there are rare and small. KIC 5520878 is the counterexample (its second
mode is at ~1% amplitude, enough for the trap), and the double-mode fraction
among RRc stars in Kepler and K2 is not small (Moskalik et al. 2015; Netzel et
al. 2023). Turning the prewhitening on for RR Lyrae, with the same "not a
harmonic of the dominant mode" rule, removes the beat before it reaches the
timer, and the removed mode is then a labelled input to `classify_secondary`
rather than a soft annotation.

## 5. What was changed (implemented 2026-09-18, same day)

1. **`sequence._squareness`** now compares the two-level word against the
   better of two sinusoids: the one at the word's period L (the old test) and
   the best one at any frequency up to the once-per-cycle Nyquist limit
   (`_best_sinusoid_rss`: a grid at a quarter of the resolution element, then a
   local refine). `decode_bit_code` returns `squareness` (the new, conservative
   value), `squareness_at_period` (the old one, kept as a diagnostic) and
   `beat_frequency` (the free sinusoid's frequency in cycle⁻¹). `seq_square`,
   which every decode route in `is_followup` and every gate in `inject.py`
   reads, follows the new value. `fields_from_code` carries the three
   diagnostics as `seq_squareness_at_period`, `seq_beat_frequency`,
   `seq_rotation_word` for the viewer; they are not persisted columns.
2. **`sequence.is_rotation_word`** implements the factor-complexity test of §4
   (at most 2m distinct cyclic substrings of each length m); `decode_bit_code`
   sets `rotation_word`, and `is_square_code` refuses a rotation word of
   `ROTATION_WORD_MIN_BITS = 21` bits or more, where a random word passes the
   bound 0.9% of the time. Below that the flag is reported but does not veto.
3. **RR Lyrae prewhitening.** `StarClass` gained `prewhiten_ground`; the
   `rrlyrae` preset is `prewhiten=3, prewhiten_ground=False`, so up to three
   independent modes are removed before timing on TESS data and never on
   ASAS-SN/ZTF, where nightly sampling puts daily aliases where a mode would
   be. `bcep` is unchanged (5 modes, every survey). Removed modes are recorded
   as before, so `classify_secondary` now reaches its hard `beat` veto for RR
   Lyrae with a detected second mode.
4. **`classify_secondary`** also vetoes (kind `beat`) a secondary period that
   matches the decoded word's `beat_frequency` or its alias `1/(1−ν)`, when the
   word is not square: that is the sinusoid the decoder was reading.
5. **Injection beat control.** `--words beat[:<nu>]` (default ν = 0.618)
   injects `frac·P·sin(2πνc)` into the time stamps instead of a code
   (`inject.inject_beat`); rows carry `n_bits = 0`, no recall, and the
   `triage`/`quantised`/`catalog` gates then measure the false-promotion rate
   of a two-mode star at each amplitude on real hosts.

Not implemented: the harmonic-sideband test on the light curve (§4), the
discriminant of last resort for a survivor.

### What the existing store needs

The squareness and rotation-word changes act on stored O–C series, so
`starseti stage1 --recompute-routes --triage` re-derives every sequence row,
triage flag and score without touching a light curve. The RR Lyrae
prewhitening acts before timing, so every TESS RR Lyrae row has to be re-timed
from its cached light curve (`starseti stage1 --survey tess --target-file
tess_rrlyrae.txt`); ground-based RR Lyrae rows are unaffected. Every stage-2
row was vetted with the old predicate and needs `starseti stage2 --force`. The
published thresholds were measured with the old gate and should be re-measured,
with the beat control alongside.

## 6. The real star (run 2026-09-19)

KIC 5520878 was in the target lists all along, under its other names: Gaia DR3
2100640814826289280 is `TIC 121329712` (TESS), `ASASSN 103080025403` and
`ZTF 727110300024203`. The TESS row had died on a transient MAST error, and the
only TESS products for it are two TASOC 30-min FFI sectors (no 2-min data), so
none of the three rows could see a per-cycle beat. To put the actual Kepler
data through the same cuts, stage 1 gained a `kepler` survey for named KIC
targets (`DENSE_SURVEYS = ("tess", "kepler")`: per-cycle timing, no alias test,
RR Lyrae prewhitening; stage 2 knows Kepler's month/quarter/year as its
sampling periods). It stitches the Kepler-pipeline short-cadence months where
they exist, else long cadence. A stage-0 row was copied from the TIC row so the
viewer links the four names, and `data/starlists/kepler_rrlyrae.txt` records
the target:

    ./starseti stage1 --survey kepler --target-file kepler_rrlyrae.txt
    ./starseti stage2 --targets "KIC 5520878" --force --catalogs

Stage 1 on 19 short-cadence months (723k points, 395 per cycle, 856 d
baseline, 1841 timed cycles):

| step | result |
|---|---|
| prewhitening | one mode removed: P = 0.17010 d, f1/f0 = 1.5825, A = 0.0071 (2.4% of the 0.30 peak-to-peak; the floor is 2%) |
| O–C | rms 111 s; Ljung–Box p = 5e-4, low permutation entropy, secondary periods at 1.289 d (4.79 cycles), 1.285 d and 205 d |
| decode | 24-bit word, `seq_square` False, squareness 0.00, unrecognized, `seq_p` 1.0 |
| triage | flagged (`TRIAGE`, score 0.32) for autocorrelation and the secondary period, not for a code |

So the trap is closed at the point it mattered: the residual beat (the mode is
removed only to the fit's precision, and its aliases 1−ν and the 4.79-cycle
subharmonic remain) still makes the O–C non-white, which is what stage 1 is
for, but the decoder no longer calls it a two-level word. Stage 2 returns **ASTRO** with two reasons:
the 205 d secondary period is in the Blazhko/multi-mode range, and the 4.79-
and 4.77-cycle secondary periods are, within a per cent, `1/beat_frequency`
of the decoded word (0.2088 per cycle): the §5.4 veto, "the sinusoid at 4.79
cycles that explains the decoded word's steps as well as two levels do: a
second mode sampled once per cycle". The 24-bit word itself
(`110101001010010101111011`) has squareness 0.00 against the free sinusoid and
0.13 against the sinusoid at its own period, and is not a rotation word,
because the residual after prewhitening is no longer a clean rotation. So on
the real star the trap is closed twice over, at the decoder and at the vet,
and the star sits in the dashboard as an ASTRO row next to its other three
names. The other three rows
after re-timing with the new code: TESS/TASOC sparse (13 points per cycle),
quiet; ASAS-SN quiet; ZTF flagged for a 1564 d secondary period with 74 cycles
per block, which is secular, not the beat.

## Script

```python
import sys, numpy as np

sys.path.insert(0, "src")
import sequence, stage1_pulse as s1


def best_sinusoid_rss(c, y):
    best = None
    for f in np.linspace(0.002, 0.5, 2000):
        A = np.column_stack(
            [np.ones_like(c, dtype=float), np.cos(2 * np.pi * f * c), np.sin(2 * np.pi * f * c)]
        )
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        r = y - A @ coef
        rss = float(r @ r)
        if best is None or rss < best[0]:
            best = (rss, f)
    return best


def two_level_rss(c, y, word, L):
    lv = {0: [], 1: []}
    for ci, yi in zip(c, y):
        lv[word[ci % L]].append(yi)
    m = {b: np.mean(v) for b, v in lv.items()}
    return float(sum((yi - m[word[ci % L]]) ** 2 for ci, yi in zip(c, y)))


def report(oc, cyc, tag):
    dec = sequence.decode_bit_code(list(oc), list(cyc))
    pts = sequence.code_points(list(oc), list(cyc))
    cs = np.array([p["c"] for p in pts])
    ys = np.array([p["y"] for p in pts])
    L = dec["n"]
    word = [int(b) for b in dec["bits"]]
    rss2 = two_level_rss(cs, ys, word, L)
    rss_sinL = sequence._sinusoid_rss(list(cs), list(ys), L)
    rss_best, fbest = best_sinusoid_rss(cs, ys)
    print(tag, L, dec.get("label"), "".join(map(str, word)))
    print("  squareness vs sinusoid at L :", round(1 - rss2 / rss_sinL, 2))
    print(
        "  squareness vs best sinusoid :",
        round(max(0.0, 1 - rss2 / rss_best), 2),
        "f =",
        round(fbest, 4),
    )


rng = np.random.default_rng(1)
c = np.arange(120)
report(30 * np.sin(2 * np.pi * c / 0.6177) + rng.normal(0, 10, 120), c, "beat 0.6177")
report(30 * np.sin(2 * np.pi * c * np.sqrt(2)) + rng.normal(0, 5, 120), c, "beat sqrt2")
lcd = s1.synthetic("tess", inject=True, star_class="rrlyrae")
res, df = s1.run_target("SYNTH rrlyrae", lcd, 200, s1.STAR_CLASSES["rrlyrae"])
report(list(df["oc_detrended_days"] * 86400.0), list(df["cycle_index"]), "staircase code")
```
