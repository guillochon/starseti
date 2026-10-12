"""matched.py -- the receiver a beacon designer expects, and the only recognizer.

A word is named by correlating the folded step means directly against every
catalog window -- the packed numeric catalogs (primes, Fibonacci, counting,
powers of two, squares, triangular, in minimal binary and 4/8-bit fixed
width) and the binary expansions of pi and e -- at the word length and bit
period the description-length search chose, and that length's tilings up to
40 bits, over every rotation. Nothing is thresholded first, so a word whose
bits are individually below the noise but collectively right still scores.
The exact-substring matcher this replaced needed every folded bit to land on
the right side of the split, and lost a 28-bit word to one bad bit.

It is a matched filter in the proper sense. The amplitude of each unit
pattern in the folded means is divided by the noise those means carry and
the best of N such amplitudes is charged N times its Gaussian tail, N
counting every window, rotation and word length searched. A hit also has to
leave residuals consistent with that noise: a periodic pattern that is not
the word correlates at 0.8 with something in any finite catalog, and a test
on the correlation alone would call a noise-free one certain. The
goodness-of-fit check is what makes this a match and not a resemblance.

Two things sit on top of the filter.

**Near misses.** A long word that differs from a catalog window in a bit or
two is the catalog word with bit errors, not a different word, and a beacon
is allowed a bit error. So every window is also tried with up to
`L // NEAR_MISS_DIVISOR` of its bits flipped to agree with the fold's own
two-level word. The pattern that results *is* the fold's own word, so it
needs no catalog charge; what it pays for is the bits (below). The hit
records how many bits differ, and the viewer shows them.

The amplitude of that word is not evidence the way a catalog window's is.
The window is a fixed pattern, so its projection on noise is one Gaussian
tail; the fold's own word is *fit* to the means, and the best two-level
split of L noise means captures most of their scatter -- at 33 residues it
stood six sigma above the noise with no signal in the data, and a Gaussian
tail charged L - 1 ways called that 1e-7. So the near-miss route's amplitude
p is the plain statement that the folded means are not all alike: their
centred sum of squares over the noise a mean carries, chi-square on L - 1
under the null. It says there is structure at this fold; the bits null says
whether that structure is a catalog word.

**The bits null.** The amplitude p says the folded means are not noise. It
does not say the word is unlikely to resemble the catalog by chance: at
twelve bits almost every word is a catalog window, and a strong random word
correlates perfectly with itself. So each hit also carries the p-value of a
chance word of the decoded length landing within its mismatch count of some
window, and reports the larger of the two. A claim has to survive both
nulls. The transcendental expansions are charged by where the match begins
(`sequence.stream_trials`), the numeric packings flat, as before.

The chance word is not a fair coin. The bits come from a two-level split of
the folded means, and a split that isolates one residue makes a word of one
set bit and L - 1 clear ones out of any fold with an outlier step. Under a
fair coin that word is as unlikely as any other, 2^-L, and the powers of two
in minimal binary -- lone ones between growing runs of zeros -- have a window
for every position the one could sit in, so a single outlier was named as a
power of two at the p of a full 28-bit primes match. The null therefore
conditions on the word's composition: given `ones` set bits in `word_len`,
the word is uniform over the C(word_len, ones) arrangements, and the chance
of landing within m flips of an alignment is the number of arrangements
within reach of that window (`reachable`) over C(word_len, ones), summed
over alignments. A spike word has 28 arrangements and 28 rotations of a
lone-one window to land on, so it claims nothing; the primes and counting
controls move by less than a factor of two.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy import stats

import sequence

MIN_LEN = sequence.SEQ_MIN_BITS_STREAM  # 12
MAX_LEN = 40

# What every trial count is multiplied by: half the prior mass to each
# corpus, and each window tried in both polarities (sequence.SEQ_POLARITIES).
_SPLIT = sequence.SEQ_CORPUS_SPLIT * sequence.SEQ_POLARITIES

# A hit has to leave residuals consistent with the measured noise: below this
# goodness-of-fit p the folded means are a pattern, but not this pattern.
GOF_MIN_P = 1e-3

# And it has to fit the means at least this fraction as well as the fold's
# own free two-level word does. Smooth modulation folded at some length
# correlates at ~0.9 with a square catalog window while its own word fits it
# far better; this is what separates naming a word from resembling one.
FREE_FRACTION = 0.97

# One bit error per this many bits of window: a 12-bit window may miss one
# catalog bit, a 28-bit window two, a 40-bit window four. A short word cannot
# afford an error because the bits null (see the module docstring) would not
# let it claim anything anyway.
NEAR_MISS_DIVISOR = 10

# An analytic tail underflows long before a strong hit stops being real;
# stored p-values are floored here so a saturated one stays distinguishable
# from "no evidence" and has a logarithm.
P_FLOOR = 1e-300

# A set bit is an early cycle: a negative step. The pattern a word predicts
# in the folded means is -1 where the bit is set and +1 where it is clear.
_SIGN = {"1": -1.0, "0": 1.0}

STREAM_ENCODING = "bitstream"


@dataclass(frozen=True)
class Windows:
    """Every distinct catalog window of one length, as +-1 rows.

    `sign` is (rows, L) in {-1, +1}; `norm` the length of each centred row,
    so `sign @ x / norm` is the amplitude of the unit pattern in a centred x.
    Numeric windows come first (`n_numeric` of them) and carry `start = -1`;
    expansion windows follow with the offset their match begins at. `labels`
    are (catalog index into sequence.KNOWN_BIT_SEQ, encoding, start).
    """

    sign: np.ndarray
    norm: np.ndarray
    start: np.ndarray
    labels: tuple[tuple[int, str, int], ...]
    n_numeric: int

    def word(self, row: int) -> str:
        return "".join("1" if s < 0 else "0" for s in self.sign[row])


def _numeric_rows(L: int) -> list[tuple[str, tuple[int, str, int]]]:
    out: list[tuple[str, tuple[int, str, int]]] = []
    seen: set[str] = set()
    for enc_name, enc in sequence.catalog_encoders():
        for idx, seq in enumerate(sequence.KNOWN_BIT_SEQ):
            if not seq["values"]:
                continue
            packed = "".join(b for b in (enc(v) for v in seq["values"]) if b)
            for start in range(0, len(packed) - L + 1):
                w = packed[start : start + L]
                if w in seen or "0" not in w or "1" not in w:
                    continue
                seen.add(w)
                out.append((w, (idx, enc_name, start)))
    return out


def _stream_rows(L: int) -> list[tuple[str, tuple[int, str, int]]]:
    """Distinct windows of the expansions, each kept at the shallowest offset
    it occurs at -- the cheapest under the offset prior, so there is no better
    alignment to look for once it is found."""
    best: dict[str, tuple[int, str, int]] = {}
    for idx, seq in enumerate(sequence.KNOWN_BIT_SEQ):
        if seq["values"]:
            continue
        bits = seq["bits"]
        for start in range(0, len(bits) - L + 1):
            w = bits[start : start + L]
            if "0" not in w or "1" not in w:
                continue
            cur = best.get(w)
            if cur is None or start < cur[2]:
                best[w] = (idx, STREAM_ENCODING, start)
    return list(best.items())


def _build_windows() -> dict[int, Windows]:
    out = {}
    for L in range(MIN_LEN, MAX_LEN + 1):
        numeric = _numeric_rows(L)
        streams = _stream_rows(L)
        rows = numeric + streams
        if not rows:
            continue
        sign = np.array([[_SIGN[ch] for ch in w] for w, _ in rows], dtype=np.int8)
        centred = sign.astype(float) - sign.mean(axis=1, keepdims=True)
        norm = np.linalg.norm(centred, axis=1)
        start = np.array([lab[2] if lab[1] == STREAM_ENCODING else -1 for _, lab in rows])
        out[L] = Windows(
            sign=sign,
            norm=norm,
            start=start,
            labels=tuple(lab for _, lab in rows),
            n_numeric=len(numeric),
        )
    return out


_WINDOWS = _build_windows()

# Every (numeric window, rotation, L) the search below looks at: the flat
# trial count a numeric hit is charged with. Expansion windows are charged
# per offset instead (sequence.stream_trials), so they are not counted here.
_ROTATED_TRIALS = sum(w.n_numeric * L for L, w in _WINDOWS.items())


def n_trials(bit_cycles: Sequence[int] = sequence.CODE_BIT_CYCLES) -> int:
    """Every (numeric window, rotation, L, k) a blind search would look at."""
    return int(_ROTATED_TRIALS) * len(bit_cycles)


def lengths_for(word_len: int) -> list[int]:
    """The window lengths tried for a word of length `word_len`: the word and
    its tilings up to MAX_LEN, since a catalog window may span several
    repeats of a short word."""
    if word_len < 1:
        return []
    return [L for L in range(word_len, MAX_LEN + 1, word_len) if L in _WINDOWS]


def trials_for(word_len: int) -> int:
    """Numeric windows times rotations at the lengths tried for `word_len`."""
    return sum(_WINDOWS[L].n_numeric * L for L in lengths_for(word_len))


def near_miss_bits(L: int) -> int:
    """How many catalog bits a window of length L may disagree with."""
    return max(0, int(L) // NEAR_MISS_DIVISOR)


def hamming_ball(n: int, radius: int) -> int:
    """Words of n bits within `radius` flips of a given one."""
    r = max(0, min(int(radius), int(n)))
    return sum(math.comb(int(n), i) for i in range(r + 1))


def reachable(n: int, window_ones: int, word_ones: int, radius: int) -> int:
    """Words of n bits with `word_ones` set bits that lie within `radius`
    flips of a fixed word with `window_ones` set bits: a set bits cleared and
    b clear bits set, a + b <= radius, window_ones - a + b = word_ones."""
    n, wp, wt, r = int(n), int(window_ones), int(word_ones), int(radius)
    total = 0
    for a in range(0, min(wp, r) + 1):
        b = wt - wp + a
        if 0 <= b <= n - wp and a + b <= r:
            total += math.comb(wp, a) * math.comb(n - wp, b)
    return total


def _prefix_ones(win: Windows, word_len: int) -> np.ndarray:
    """Set bits in the first `word_len` bits of each window: the part of a
    tiled window a word of that length has to match."""
    return (win.sign[:, :word_len] < 0).sum(axis=1)


def _bits_length(word_len: int) -> int | None:
    """The one length the numeric bits null counts windows at: the shortest
    searched for the word. A tiled window can only match a tiled word if its
    first `word_len` bits do, and those are a window at the shorter length,
    so the longer lengths add no alignments the word could land on."""
    lengths = lengths_for(int(word_len))
    return lengths[0] if lengths else None


def bits_alignments(word_len: int, ones: int, mismatch: int) -> float:
    """Numeric alignments a word of `word_len` bits with `ones` set bits can
    land within `mismatch` flips of, weighted by how many of its arrangements
    reach each: the distinct windows at the shortest length searched, each
    counted `reachable` times, times the word's rotations and the corpus
    split. Over C(word_len, ones) this is the bits p (`bits_pvalue`).
    """
    L0 = _bits_length(word_len)
    if L0 is None:
        return 0.0
    win = _WINDOWS[L0]
    hist = np.bincount(_prefix_ones(win, word_len)[: win.n_numeric], minlength=word_len + 1)
    hits = sum(int(h) * reachable(word_len, wp, ones, mismatch) for wp, h in enumerate(hist) if h)
    return _SPLIT * int(word_len) * hits


def bits_trials(word_len: int) -> float:
    """Every numeric alignment a word of `word_len` bits is tried against,
    composition aside: distinct windows at the shortest length searched times
    rotations times the corpus split. The flat count the null used before it
    conditioned on the word's set bits; kept for the trial bookkeeping."""
    L0 = _bits_length(word_len)
    if L0 is None:
        return 0.0
    return _SPLIT * _WINDOWS[L0].n_numeric * int(word_len)


def bits_pvalue(word_len: int, ones: int, alignments: float) -> float:
    """P(a chance word of `word_len` bits with `ones` set bits lands on one of
    `alignments` weighted alignments), the word uniform over its
    C(word_len, ones) arrangements: weighted Bonferroni, `alignments`
    already carrying the reach of each window, the corpus split and the
    offset prior where they apply. A word with every bit alike is no word.
    """
    word_len, ones = int(word_len), int(ones)
    if word_len < 1 or alignments <= 0 or ones <= 0 or ones >= word_len:
        return 1.0
    return min(1.0, float(alignments) / math.comb(word_len, ones))


def _fold(y: np.ndarray, cyc: np.ndarray, L: int, k: int) -> np.ndarray | None:
    j = (cyc // k) % L
    cnt = np.bincount(j, minlength=L)
    if np.any(cnt < 1):
        return None
    return np.bincount(j, weights=y, minlength=L) / cnt


def _empty() -> dict:
    return {"p": 1.0, "rho": None, "label": None, "len": 0, "k": 0, "z": None}


def format_label(
    name: str, encoding: str, start: int, mismatch: int = 0, inverted: bool = False
) -> str:
    """The stage-1 label: catalog, packing, where the window begins, how many
    bits the word disagrees with it on, and whether it matched inverted."""
    where = f"{encoding}, bit {start + 1}"
    if mismatch:
        where += f", {mismatch} bit{'s' if mismatch > 1 else ''} off"
    if inverted:
        where += ", inverted"
    return f"{name} ({where})"


def matched_catalog(y: np.ndarray, cyc: np.ndarray, word_len: int, k: int = 1) -> dict:
    """Correlate the folded steps against the catalog windows of the word's length.

    `y` are the O-C steps (increments) and `cyc` their cycles, in any order;
    `word_len` and `k` are the word length and bit period the description-
    length search chose (`_fold_cost_nb` / `sequence._fold_word`). The filter
    runs at that fold and its tilings up to MAX_LEN only. Folding at a length
    that is not the word's mixes bits within a residue, which inflates the
    within-residue noise estimate with signal and lets an aliased pattern pass
    the goodness-of-fit check while still standing eight sigma above the
    inflated noise -- a tiled random word "matched" the counting sequence that
    way. At the word's own length the fold is clean and the noise honest.

    Returns the best hit: its p-value (the larger of the amplitude p over
    every window, rotation and length tried, and the bits p of the word
    resembling the catalog by chance; see the module docstring), the
    correlation, which catalog window it is and where, and how many bits the
    fold's own word disagrees with that window on. `p` is 1.0 when nothing
    could be folded or the word is too short to name.
    """
    y = np.asarray(y, dtype=float)
    cyc = np.asarray(cyc, dtype=np.int64)
    best = _empty()
    k = max(1, int(k))
    word_len = int(word_len)
    lengths = lengths_for(word_len)
    if not lengths or y.size < 2 * MIN_LEN:
        return best
    n_len = len(lengths)
    numeric_trials = _SPLIT * trials_for(word_len)
    radii = range(word_len + 1)
    exact_hits = [bits_alignments(word_len, wp, 0) for wp in radii]
    for L in lengths:
        win = _WINDOWS[L]
        if y.size < 2 * L * k:
            continue
        m = _fold(y, cyc, L, k)
        if m is None:
            continue
        # The noise a folded mean carries. Not the within-residue scatter at
        # this fold: the fold was *chosen* for compressing well, so its
        # residues are the ones that happened to be tight, and the scatter
        # they leave understates the noise (pure noise fired at p < 1e-3 six
        # times in a hundred that way). The total variance of the steps is an
        # honest upper bound -- signal only adds to it -- and costs a strong
        # code some z it can spare.
        n_bar = y.size / L
        s2 = float(np.var(y))
        sigma_a = math.sqrt(s2 / n_bar) if s2 > 0 else 0.0
        mc0 = m - m.mean()
        norm = float(np.linalg.norm(mc0))
        if not norm > 0 or not sigma_a > 0:
            continue
        # The fold's own two-level word, with no catalog in sight: a catalog
        # window has to explain the means about as well as that, or it is a
        # resemblance to structure the word does not have.
        cut = sequence.split_means(list(m))
        if cut is None:
            continue
        free0 = np.where(m < cut, -1.0, 1.0)
        fc0 = free0 - free0.mean()
        fnorm = float(np.linalg.norm(fc0))
        if not fnorm > 0:
            continue
        # The free word's fit does not depend on which symbol is called set.
        amp_free = float(np.dot(mc0, fc0) / fnorm)
        rho_free = amp_free / norm
        dof = max(L - 2, 1)
        idx = (np.arange(L)[:, None] - np.arange(L)[None, :]) % L
        S = win.sign.astype(float)
        tiles = L // word_len
        # What each alignment's amplitude is charged: flat over the numeric
        # packings, by offset over the expansions, both times the corpus
        # split and the two polarities, over every rotation and length
        # tried. The bits null charges the word's own length only
        # (`bits_trials`), and the expansions by offset times the rotations.
        is_stream = win.start >= 0
        stream_offsets = sequence.stream_trials(np.maximum(win.start, 0))
        charge = np.where(is_stream, _SPLIT * stream_offsets * L * n_len, numeric_trials)[:, None]
        # The bits null conditions on composition (module docstring). Each
        # window's set bits over the word's length say what a word of a given
        # composition can reach; the numeric alignments are summed once over
        # the windows at the word's own length, the expansions charged per
        # window by offset. The word's set bits on route B are the fold's own
        # word's, per repeat when the window is a tiling of it.
        win_ones = _prefix_ones(win, word_len)
        stream_charge = _SPLIT * stream_offsets * word_len
        p_bits_a = np.array(
            [
                bits_pvalue(word_len, wp, stream_charge[i] if is_stream[i] else exact_hits[wp])
                for i, wp in enumerate(win_ones)
            ]
        )[:, None]
        allowed = near_miss_bits(L)
        chi2_free = max(norm * norm - amp_free * amp_free, 0.0) / sigma_a**2
        gof_free = float(stats.chi2.sf(chi2_free, dof))
        ok_free = amp_free > 0 and gof_free >= GOF_MIN_P
        p_amp_b = float(stats.chi2.sf(norm * norm / sigma_a**2, L - 1))
        # Both polarities: the same means with every bit flipped are the same
        # catalog word sent by a transmitter that kicks on the other symbol.
        # Everything above is polarity-free; below, the means and the fold's
        # own word turn over together and the window's mismatch to the word
        # becomes L minus what it was.
        for polarity in (1.0, -1.0):
            mc = polarity * mc0
            free = polarity * free0
            # Every rotation of the folded word at once: column r is mc rolled by r.
            R = mc[idx]  # (L, L): R[:, r] = roll(mc, r)
            F = free[idx]  # the same rotations of the fold's own word
            amp = (S @ R) / win.norm[:, None]  # amplitude of each unit pattern
            mismatch = np.rint((L - S @ F) / 2).astype(int)  # bits vs the fold's word
            word_mismatch = -(-mismatch // tiles)  # per word, not per tiled window
            free_ones = int(np.clip(round(float((free < 0).sum()) / tiles), 0, word_len))
            # Route A: the window itself. Amplitude against the noise,
            # residuals consistent with it, and at least a fraction of the fit
            # the fold's own word gets. Its bits null is the exact alignment
            # of a word with the window's own composition.
            z = amp / sigma_a
            chi2 = np.maximum(norm * norm - amp * amp, 0.0) / sigma_a**2
            gof = stats.chi2.sf(chi2, dof)
            ok_a = (amp > 0) & (gof >= GOF_MIN_P) & (amp / norm >= FREE_FRACTION * rho_free)
            p_amp_a = charge * stats.norm.sf(z)
            p_a = np.where(ok_a, np.maximum(p_amp_a, p_bits_a), np.inf)
            # Route B: the window with its disagreeing bits flipped, which is
            # the fold's own word. That word is fit to the means, so its
            # projection is no Gaussian tail (module docstring): the
            # amplitude p is the chi-square of the centred means over the
            # noise, L - 1 degrees of freedom -- the means are not all alike
            # -- and the word has to be an unlikely thing to lie this close
            # to the catalog by chance.
            ok_b = ok_free & (mismatch >= 1) & (mismatch <= allowed)
            wm = np.clip(word_mismatch, 0, word_len)
            numeric_hits = np.array([bits_alignments(word_len, free_ones, mm) for mm in radii])
            reach_table = np.array(
                [[reachable(word_len, wp, free_ones, mm) for mm in radii] for wp in radii],
                dtype=float,
            )
            aligns_b = np.where(
                is_stream[:, None],
                stream_charge[:, None] * reach_table[win_ones[:, None], wm],
                numeric_hits[wm],
            )
            # bits_pvalue, vectorised: a word with every bit alike is no word.
            arrangements = math.comb(word_len, free_ones) if 0 < free_ones < word_len else 0
            p_bits_b = (
                np.where(aligns_b > 0, np.minimum(1.0, aligns_b / arrangements), 1.0)
                if arrangements
                else np.ones_like(aligns_b)
            )
            p_b = np.where(ok_b, np.maximum(p_amp_b, p_bits_b), np.inf)
            p_all = np.minimum(p_a, p_b)
            flat = int(np.argmin(p_all))
            i, r = np.unravel_index(flat, p_all.shape)
            p_here = float(p_all[i, r])
            if not math.isfinite(p_here) or p_here >= best["p"]:
                continue
            via_a = p_a[i, r] <= p_b[i, r]
            cat_idx, enc, start = win.labels[i]
            name = sequence.KNOWN_BIT_SEQ[cat_idx]["name"]
            mm = int(mismatch[i, r])
            inverted = polarity < 0
            best = {
                "p": max(P_FLOOR, min(1.0, p_here)),
                "p_amplitude": float(p_amp_a[i, r] if via_a else p_amp_b),
                "p_bits": float(p_bits_a[i, 0] if via_a else p_bits_b[i, r]),
                "rho": float(amp[i, r] / norm) if via_a else float(rho_free),
                "label": format_label(name, enc, start, 0 if via_a else mm, inverted),
                "name": name,
                "encoding": enc,
                "catalog": int(cat_idx),
                "start": int(start),
                "len": int(L),
                "k": int(k),
                "rotation": int(r),
                "z": float(z[i, r] if via_a else amp_free / sigma_a),
                "gof": float(gof[i, r] if via_a else gof_free),
                # The catalog's own bits; `sequence.recognize_hit` turns them
                # into the word's polarity for the viewer.
                "window": win.word(i),
                # -1: the window matched the word with every bit flipped.
                "polarity": -1 if inverted else 1,
                # How many of the window's bits the fold's own word disagrees
                # with, in the polarity that matched. Zero on an exact word;
                # on route A it can be positive while the correlation still
                # carries the match, which is the below-the-noise case the
                # filter exists for.
                "mismatch": mm,
                "corrected": not via_a,
            }
    return best
