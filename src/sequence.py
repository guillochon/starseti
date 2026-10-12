"""O–C bit-sequence recognition. Matches candidates.html.

The waveform
------------
A pulsator is a limit-cycle oscillator, and phase is the marginal direction
along the cycle: an amplitude kick damps back, a phase kick stays. So a cycle
triggered early by dt does not snap back when the next bit is clear -- the
clock has been reset, and every later pulse arrives dt earlier than it would
have. The O-C of a one-sided phase code is therefore the running sum of the
bits, a **staircase**, not a square wave; see research/modulation-physics.md
section 2. The bits live in the *first difference* of the O-C, one step per
cycle, and that is what `decode_bit_code` thresholds by default
(`model="staircase"`). The level-thresholding decoder is kept as
`model="square"` for the waveform the pipeline assumed before, and for tests
that want to see the two side by side.

Naming the word
---------------
The decoded word is named by `matched.matched_catalog`, a matched filter over
every catalog window that never asks the bits to be exactly right: a word a
bit or two off a long catalog window is that window with bit errors, and a
word whose bits are individually below the noise but collectively right still
scores. This module turns its hit into the label, digits and spans the viewer
draws, and rotates the word so it reads in catalog order.
"""

from __future__ import annotations

import math
from bisect import bisect_left
from collections.abc import Callable, Iterable, Sequence
from itertools import pairwise
from typing import Any

import numpy as np

MAX_OC = 480

# How the O-C carries the bits. `staircase` -- the physical model -- reads the
# bit of cycle c off the step O-C(c) - O-C(c-1); `square` reads it off the
# level of O-C(c) itself. `fold` is the staircase read from a series with no
# consecutive cycles (one timing per observed cycle of a ground-based survey):
# the *levels* are folded on the word first and the bits are the steps
# between neighbouring residue means, so a word repeated r times is still read
# at sqrt(r) the per-cycle signal-to-noise and no two adjacent cycles ever
# have to have been observed.
CODE_MODELS = ("staircase", "square", "fold")
DEFAULT_CODE_MODEL = "staircase"
SEQ_MIN_BITS_STREAM = 12
SEQ_MAX_BITS = 64
SPARK_W = 72
SPARK_H = 18

# Which decoder wrote a stored word. Reading the O-C as a bit word is the
# whole cost of a viewer request -- seconds on a long series -- so the word is
# decoded once, when stage 1 writes the star, and the server serves what was
# stored. That is only honest while the stored word is the word this code
# would produce, so every row carries this number and a row stamped with an
# older one is decoded live instead.
#
# Bump it in the same commit as any change to `decode_bit_code`, to what it
# reads (`code_points`, `clip_steps`), to the catalogs it matches against, to
# `matched.matched_catalog` or to `recognize_hit`.
#
# Forgetting is not a slow path, it is a wrong one, and the asymmetry is worth
# stating plainly. A row carrying an *older* stamp is decoded live: slow, and
# right. A row carrying *this* stamp is served as authoritative -- so the rows
# already in the store, which keep whatever number they were written with, go
# on being served as a word the current decoder would not produce, and nothing
# anywhere says so. The bump is what turns that into the slow, right case.
#
# 3 is the winsorized steps: `clip_steps` changed what `code_points` hands the
# decoder without a bump, which is exactly the trap above, and 14,686 stored
# rows claimed to be current while holding a word from before it.
# 4: `clip_steps` leaves a one-sided cluster of outlying steps alone (a
# low-duty-cycle code's set bits), the bit-period grid grows with the series
# (`bit_cycles_for`), and the `fold` model reads a level series that has no
# consecutive cycles.
# 5: the matched filter names a word in either polarity (SEQ_POLARITIES),
# and every catalog p-value carries the factor of two for it.
SEQ_CODE_VERSION = 5

EMPTY = {
    "seq_matched": False,
    "seq_name": "",
    "seq_label": "unrecognized",
    "seq_digit_count": 0,
    "seq_sort_n": 0.0,
    "oc_spark": "",
    "seq_p": 1.0,
    "seq_bit_len": 0,
    "seq_squareness": None,
    "seq_square": False,
    # The same test asked of the catalog window the word is named for, which
    # a fixed pattern earns at the two repeats the matched filter needs rather
    # than the three the fold's own word does (is_square_match).
    "seq_window_squareness": None,
    "seq_window_square": False,
    # Diagnostics behind seq_square: the squareness against a sinusoid at the
    # word's own period (what the test used to compare with), the frequency of
    # the free sinusoid that explains the steps, and the rotation-word flag.
    "seq_squareness_at_period": None,
    "seq_beat_frequency": None,
    "seq_rotation_word": False,
    # One character per pip the viewer draws (one per catalog value, or per
    # four bits of an expansion): "1" where a near-miss bit sits inside it.
    "seq_off_mask": "",
    # The decoded word itself, kept beside the fields derived from it so the
    # caller that persists a star can store the word and the viewer need never
    # decode it again. None when the series carries no word.
    "code": None,
}

# A two-level model needs this much of the sinusoid's residual variance to
# count as squarer than a sinusoid, and the folded word this much
# self-consistency, before stage 2 will treat a coherent secondary period as a
# candidate phase code rather than as astrophysics. At these values the
# two-level model has to beat the sinusoid by roughly a 1.3-sigma amplitude.
CODE_MIN_SQUARENESS = 0.25
CODE_MIN_AGREEMENT = 0.75

# Pulsation cycles per bit the fold searches. A beacon designer facing unknown
# survey cadences may hold each bit for several cycles; the fold indexes the
# word on cycle // k, so k = 1 is the one-bit-per-cycle picture and the others
# come back as words of their own length rather than stretched copies. The
# search over k, like the one over L, is charged inside stage 1's shuffle null
# of the decode statistic, which redraws it.
#
# This is the floor of the grid: `bit_cycles_for` extends it in powers of two
# for as long as the shortest word (4 bits) can still repeat twice, so a
# 3,000-cycle series is searched out to 256 cycles per bit and a 40-cycle one
# stops where it always did. A code held for a night, a week or a month per
# bit -- what a transmitter facing nightly ground-based sampling would send --
# was invisible at 8.
CODE_BIT_CYCLES = (1, 2, 4, 8)
CODE_MIN_WORD_BITS = 4


def bit_cycles_for(n_points: int) -> tuple[int, ...]:
    """The bit periods the fold searches on a series of `n_points` steps.

    Powers of two from 1 while a `CODE_MIN_WORD_BITS`-bit word at that period
    still fits twice into the series; never shorter than `CODE_BIT_CYCLES`,
    so a short series searches exactly what it did before (the extra periods
    fold nothing there and were dropped by the fold anyway).
    """
    ks = list(CODE_BIT_CYCLES)
    k = ks[-1] * 2
    while 2 * CODE_MIN_WORD_BITS * k <= int(n_points):
        ks.append(k)
        k *= 2
    return tuple(ks)


def _u32_bits(v: float | int) -> str:
    return format(int(v) & 0xFFFFFFFF, "b")


def _primes_to(n: int) -> list[int]:
    mark = [False] * (n + 1)
    out = []
    for i in range(2, n + 1):
        if mark[i]:
            continue
        out.append(i)
        for j in range(i * i, n + 1, i):
            mark[j] = True
    return out


def _fib_list(n: int) -> list[int]:
    out = [1, 1]
    while len(out) < n:
        out.append(out[-1] + out[-2])
    return out


def _concat_binary(values: Iterable[int], max_bits: int) -> str:
    s = ""
    for v in values:
        s += _u32_bits(v)
        if len(s) >= max_bits:
            break
    return s


# ---- binary expansions of the transcendental constants --------------------
# These were 96-bit literals, which is an arbitrary depth and it was the
# *search* depth: a transmission of pi caught past its 96th bit could not match
# at all. They are computed instead, as deep as SEQ_STREAM_BITS, which the
# offset prior in matched.matched_catalog makes free -- see research/sequence-pvalue.md.
# Integer arithmetic only, so there is no float precision ceiling, and the
# tests pin the first 96 bits of each against the literals they replace.

SEQ_STREAM_BITS = 4096
_GUARD_BITS = 64  # absorbs the floor() of every series term


def _atan_inv_scaled(x: int, bits: int) -> int:
    """floor(atan(1/x) * 2**bits) for integer x >= 2, by the Gregory series."""
    total = 0
    sign = 1
    k = 0
    x2 = x * x
    pow_x = x  # x**(2k+1)
    while True:
        term = (1 << bits) // ((2 * k + 1) * pow_x)
        if term == 0:
            return total
        total += sign * term
        sign = -sign
        k += 1
        pow_x *= x2


def _pi_scaled(bits: int) -> int:
    """floor(pi * 2**bits), by Machin's formula."""
    b = bits + _GUARD_BITS
    return (16 * _atan_inv_scaled(5, b) - 4 * _atan_inv_scaled(239, b)) >> _GUARD_BITS


def _e_scaled(bits: int) -> int:
    """floor(e * 2**bits), by the factorial series."""
    b = bits + _GUARD_BITS
    total = 0
    term = 1 << b
    k = 0
    while term:
        total += term
        k += 1
        term //= k
    return total >> _GUARD_BITS


def _expansion_bits(scaled: int, nbits: int) -> str:
    """The leading `nbits` binary digits of a constant, integer part included."""
    return bin(scaled)[2:][:nbits]


def _transcendental_streams(nbits: int = SEQ_STREAM_BITS) -> dict[str, str]:
    """pi with and without its integer part, and e. The fractional-only variant
    exists because a transmitter has no "binary point" to send, so a beacon
    sending 3.14159... and one sending .14159... are both plausible alignments.
    """
    pi = _expansion_bits(_pi_scaled(nbits), nbits)
    return {"pi": pi, "pi_frac": pi[2:], "e": _expansion_bits(_e_scaled(nbits), nbits)}


_STREAMS = _transcendental_streams()


def _known_bit_seq(nbits: int | None = None) -> list[dict]:
    """The catalogs the matcher searches. `nbits` is the depth of the
    transcendental expansions; it is a parameter so a test can show that
    changing it moves no p-value, which is the point of the offset prior."""
    streams = _STREAMS if nbits is None else _transcendental_streams(nbits)
    primes = _primes_to(300)
    fib = _fib_list(40)
    fib_from2 = _fib_list(40)[1:]
    count = list(range(1, 81))
    pow2 = [1 << i for i in range(16)]
    squares = [(i + 1) * (i + 1) for i in range(40)]
    tri = [(i + 1) * (i + 2) // 2 for i in range(40)]
    return [
        {"name": "primes", "bits": _concat_binary(primes, 96), "values": primes},
        {"name": "Fibonacci", "bits": _concat_binary(fib, 96), "values": fib},
        {"name": "Fibonacci", "bits": _concat_binary(fib_from2, 96), "values": fib_from2},
        {"name": "counting", "bits": _concat_binary(count, 96), "values": count},
        {"name": "powers of two", "bits": _concat_binary(pow2, 96), "values": pow2},
        {"name": "squares", "bits": _concat_binary(squares, 96), "values": squares},
        {"name": "triangular", "bits": _concat_binary(tri, 96), "values": tri},
        {"name": "π", "bits": streams["pi"], "values": None},
        {"name": "π", "bits": streams["pi_frac"], "values": None},
        {"name": "e", "bits": streams["e"], "values": None},
    ]


KNOWN_BIT_SEQ = _known_bit_seq()


def spark_path(values: Iterable[float | None], w: int = SPARK_W, h: int = SPARK_H) -> str:
    ys = [v for v in values if v is not None and math.isfinite(v)]
    if len(ys) > 48:
        step = math.ceil(len(ys) / 48)
        ys = [y for i, y in enumerate(ys) if i % step == 0]
    if len(ys) < 2:
        return ""
    lo, hi = min(ys), max(ys)
    span = 1.0 if hi == lo else hi - lo
    parts = []
    n = len(ys)
    for i, y in enumerate(ys):
        x = 0.5 + (w - 1) * i / (n - 1)
        yy = h - 0.5 - (h - 1) * (y - lo) / span
        parts.append("{}{:.1f},{:.1f}".format("L" if i else "M", x, yy))
    return " ".join(parts)


def _solve_sym3(m: Sequence[Sequence[float]], v: Sequence[float]) -> list[float] | None:
    """Solve a 3x3 system by Gaussian elimination with partial pivoting."""
    a = [[*list(row), v[i]] for i, row in enumerate(m)]
    for col in range(3):
        piv = max(range(col, 3), key=lambda r: abs(a[r][col]))
        if abs(a[piv][col]) < 1e-300:
            return None
        a[col], a[piv] = a[piv], a[col]
        for r in range(col + 1, 3):
            f = a[r][col] / a[col][col]
            for k in range(col, 4):
                a[r][k] -= f * a[col][k]
    out = [0.0, 0.0, 0.0]
    for r in (2, 1, 0):
        s = a[r][3] - sum(a[r][k] * out[k] for k in range(r + 1, 3))
        out[r] = s / a[r][r]
    return out


def _sinusoid_rss(cs: Sequence[float], ys: Sequence[float], period: float) -> float | None:
    """Residual sum of squares of a + b cos(2 pi c / L) + d sin(2 pi c / L)."""
    n = len(ys)
    if n < 4 or period < 2:
        return None
    basis = []
    for c in cs:
        ang = 2.0 * math.pi * (c % period) / period
        basis.append((1.0, math.cos(ang), math.sin(ang)))
    m = [[0.0] * 3 for _ in range(3)]
    v = [0.0] * 3
    for b, y in zip(basis, ys):
        for i in range(3):
            v[i] += b[i] * y
            for j in range(3):
                m[i][j] += b[i] * b[j]
    coef = _solve_sym3(m, v)
    if coef is None:
        return None
    rss = 0.0
    for b, y in zip(basis, ys):
        r = y - (coef[0] * b[0] + coef[1] * b[1] + coef[2] * b[2])
        rss += r * r
    return rss


def _residue(c: int, L: int, k: int = 1) -> int:
    """Which bit of an L-bit word at k cycles per bit cycle c carries."""
    return int((((c // k) % L) + L) % L)


def _best_sinusoid_rss(
    cs: Sequence[float], ys: Sequence[float], fmax: float = 0.5
) -> tuple[float, float] | None:
    """(RSS, frequency) of the best a + b cos + d sin at *any* frequency.

    The series has one point per pulsation cycle, so frequencies run from one
    cycle over the span up to the Nyquist limit of half a cycle^-1. A coarse
    grid at a quarter of the resolution element is scanned, then the best peak
    is refined. This is what a second pulsation mode sampled once per cycle
    looks like -- a sinusoid at the aliased beat frequency -- and it is the
    comparison a two-level word has to beat (research/irrational-beating-trap.md).
    """
    c = np.asarray(cs, dtype=float)
    y = np.asarray(ys, dtype=float)
    n = c.size
    span = float(c.max() - c.min()) if n else 0.0
    if n < 6 or span <= 0:
        return None
    fmin = 1.0 / span
    if not (fmax > fmin):
        return None
    step = 0.25 / span
    freqs = np.arange(fmin, fmax + 0.5 * step, step)

    def rss_at(f: float) -> float:
        ang = 2.0 * math.pi * f * c
        basis = np.column_stack([np.ones(n), np.cos(ang), np.sin(ang)])
        coef, *_ = np.linalg.lstsq(basis, y, rcond=None)
        r = y - basis @ coef
        return float(r @ r)

    coarse = np.array([rss_at(f) for f in freqs])
    i = int(np.argmin(coarse))
    best_f, best_rss = float(freqs[i]), float(coarse[i])
    lo = freqs[max(0, i - 1)]
    hi = freqs[min(freqs.size - 1, i + 1)]
    for f in np.linspace(lo, hi, 25):
        r = rss_at(float(f))
        if r < best_rss:
            best_rss, best_f = r, float(f)
    return best_rss, best_f


def _squareness(
    pts: Sequence[dict], word: Sequence[int], period: int, k: int = 1, min_repeats: int = 3
) -> dict | None:
    """How much better a two-level word explains the series than a sinusoid.

    `min_repeats` is how many times the word has to be sampled before the
    comparison is fair: three for the fold's own word, whose bit assignment
    was fit to the same means it is then judged on and over-explains them at
    fewer; two for a catalog window (`window_squareness`), whose bits were
    fixed before the data were seen.

    A phase code is a square wave: the O-C sits at two levels with fast
    transitions between them. Blazhko modulation, light-travel time in a binary
    and beating between pulsation modes are all smooth, so a sinusoid describes
    them at least as well. Returns the fraction of the sinusoid's residual
    variance that the two-level model removes -- 0 means the sinusoid is no
    worse, and a pure square wave approaches 1.

    The sinusoid is the better of two: one at the word's own period, and the
    best one at *any* frequency up to the once-per-cycle Nyquist limit. The
    second is what closes the KIC 5520878 trap: a second mode with an
    irrational period ratio nu folds cleanly at every convergent denominator q
    of nu (nu q is within 1/q of an integer), so its fold at q is a clean
    sinusoid over q residues that two levels fit better than a *period-q*
    sinusoid does -- but a sinusoid at the aliased beat period 1/(1 - nu)
    explains all of it. Measured: squareness 0.76-0.79 against the period-L
    sinusoid alone, 0.00 against the free one; an injected staircase keeps 0.93
    either way (research/irrational-beating-trap.md).

    Returns {"squareness", "squareness_at_period", "beat_frequency"}; the last
    is the frequency (cycle^-1) of the free sinusoid, which stage 2 compares
    with the star's secondary periods.
    """
    if not word or period < 2 or len(pts) < min_repeats * period:
        return None
    groups = {0: [], 1: []}
    for pt in pts:
        groups[int(word[_residue(pt["c"], period, k)])].append(pt["y"])
    if not groups[0] or not groups[1]:
        return None
    level = {b: sum(vals) / len(vals) for b, vals in groups.items()}
    rss_sq = 0.0
    for pt in pts:
        r = pt["y"] - level[int(word[_residue(pt["c"], period, k)])]
        rss_sq += r * r
    cs = [pt["c"] for pt in pts]
    ys = [pt["y"] for pt in pts]
    rss_sin = _sinusoid_rss(cs, ys, period * k)
    if rss_sin is None or rss_sin <= 0:
        return None
    at_period = _clamp01(1.0 - rss_sq / rss_sin)
    free = _best_sinusoid_rss(cs, ys)
    if free is None or free[0] <= 0:
        return {"squareness": at_period, "squareness_at_period": at_period, "beat_frequency": None}
    rss_free, f_free = free
    return {
        "squareness": _clamp01(1.0 - rss_sq / min(rss_sin, rss_free)),
        "squareness_at_period": at_period,
        "beat_frequency": f_free,
    }


# A rotation word this long or longer is called a beat rather than a message.
# Random words satisfy the complexity bound below with probability 80% at 8
# bits, 8% at 16, 0.9% at 21 and 0.2% at 24 (20,000 draws each), so under the
# floor the flag is reported but does not veto.
ROTATION_WORD_MIN_BITS = 21


def is_rotation_word(bits: Sequence[int] | str) -> bool:
    """Is the word the coding of an irrational rotation, i.e. a sampled beat?

    Thresholding a sinusoid sampled at an irrational frequency nu gives the
    coding of the rotation c -> nu c mod 1 by an interval: a *rotation word*
    (a Sturmian word when the interval length equals nu, a Rote word
    otherwise). Every such word has at most 2m distinct substrings of each
    length m, against the 2^m a message can use, and the fold is cyclic, so
    the count is over cyclic substrings. The decoded words of the KIC 5520878
    analogue (21 and 34 bits, golden ratio) and of a sqrt(2) beat (29 bits)
    all satisfy the bound; the primes, Fibonacci and counting words all break
    it by m = 3 (research/irrational-beating-trap.md). A constant word is
    excluded: that is one level, not a rotation.
    """
    b = "".join("1" if str(x) in ("1", "True") else "0" for x in bits)
    L = len(b)
    if L < 3 or b.count("1") in (0, L):
        return False
    doubled = b + b
    for m in range(1, L):
        if len({doubled[i : i + m] for i in range(L)}) > 2 * m:
            return False
    return True


def is_square_code(code: dict | None) -> bool:
    """Does a decoded word look like a phase code rather than a smooth wave?

    Two levels have to beat the best sinusoid at any frequency, the per-cycle
    votes have to agree with the fold, and a long rotation word -- the coding
    of an irrational period ratio, with at most 2m distinct substrings of each
    length m -- is a beat, whatever its squareness came out as.
    """
    if not code:
        return False
    sq = code.get("squareness")
    if code.get("rotation_word") and int(code.get("n") or 0) >= ROTATION_WORD_MIN_BITS:
        return False
    return bool(
        sq is not None
        and sq >= CODE_MIN_SQUARENESS
        and (code.get("score") or 0) >= CODE_MIN_AGREEMENT
    )


# The matched filter names a word once it has been sampled twice
# (`matched.matched_catalog` folds nothing shorter than 2 L k steps), so this
# is the floor at which every named word has a window squareness.
MATCH_MIN_REPEATS = 2


def window_squareness(pts: Sequence[dict], hit: dict | None, k: int = 1) -> float | None:
    """How much better the catalog window a word is named for explains the
    series than the best sinusoid.

    `_squareness` on the fold's own word needs the word repeated three times,
    because the two-level split was fit to the same folded means it is then
    judged on and over-explains noise at fewer. A catalog window is a fixed
    pattern -- the search over windows and rotations is charged in `seq_p`,
    not fit here -- so two levels against a sinusoid is a fair contest as soon
    as the matched filter has enough of the word to name it, which is two
    repeats. Window bit i is fold residue (i - rotation) mod L (`_rotate_word`),
    so the window is turned back into fold order first. None when the hit
    names nothing or the series is too short.
    """
    if not hit or not hit.get("window"):
        return None
    window = str(hit["window"])
    L = int(hit.get("len") or len(window))
    if L < 2 or len(window) < L:
        return None
    r = int(hit.get("rotation") or 0) % L
    word = [1 if window[(i + r) % L] == "1" else 0 for i in range(L)]
    sq = _squareness(pts, word, L, max(1, int(k)), min_repeats=MATCH_MIN_REPEATS)
    return None if sq is None else sq["squareness"]


def is_square_match(code: dict | None) -> bool:
    """Does the catalog window a word is named for beat a sinusoid?

    The guard the catalog triage route needs. `is_square_code` asks the same
    of the fold's own word and cannot answer below three repeats, which is
    where most named words on block-averaged series sit (2.0-2.6 repeats in
    the store); a smooth wave that happens to correlate with a catalog window
    would otherwise reach triage with no two-level test at all. A long
    rotation word is a beat whatever the window's squareness came out as.
    """
    if not code or not code.get("label") or code.get("label") == "unrecognized":
        return False
    if code.get("rotation_word") and int(code.get("n") or 0) >= ROTATION_WORD_MIN_BITS:
        return False
    ws = code.get("window_squareness")
    return ws is not None and ws >= CODE_MIN_SQUARENESS


def _percentile(sorted_vals: Sequence[float], p: float) -> float | None:
    n = len(sorted_vals)
    if not n:
        return None
    i = (n - 1) * p
    lo = math.floor(i)
    hi = math.ceil(i)
    if lo == hi:
        return sorted_vals[lo]
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (i - lo)


def two_level_mid(values: Iterable[float | None]) -> float | None:
    """The level a two-level word is split at: midway between the 20th and 80th
    percentiles, so a lopsided duty cycle does not drag the threshold onto one
    of the levels. `None` when the series has no spread to split.

    Public because stage 1 calibrates the decode against a shuffle null and has
    to threshold the draws exactly the way `decode_bit_code` thresholds the
    observation, or the null is of a different statistic.
    """
    ys = sorted(float(v) for v in values if v is not None and math.isfinite(float(v)))
    lo, hi = _percentile(ys, 0.2), _percentile(ys, 0.8)
    if lo is None or hi is None or not (hi > lo):
        return None
    return (lo + hi) / 2


def _clamp01(x: float) -> float:
    if x < 0:
        return 0.0
    if x > 1:
        return 1.0
    return x


def _format_seq_label(
    name: str, nums: Sequence[int] | None, partial_start: bool, partial_end: bool
) -> str:
    items = []
    if partial_start:
        items.append("u")
    for n in nums or []:
        items.append(str(n))
    if partial_end:
        items.append("u")
    return name + " [" + ", ".join(items) + "]"


def _encode_minimal(v: int) -> str:
    return _u32_bits(v)


def _encode_fixed(w: int) -> Callable[[int], str | None]:
    max_v = 1 << w

    def enc(v: int) -> str | None:
        if v < 0 or v >= max_v:
            return None
        return _u32_bits(v).zfill(w)

    return enc


def catalog_encoders() -> tuple[tuple[str, Callable[[int], str | None]], ...]:
    """The packings a numeric catalog is searched in, by name."""
    return (
        ("minimal binary", _encode_minimal),
        ("4-bit", _encode_fixed(4)),
        ("8-bit", _encode_fixed(8)),
    )


# How the prior mass is spread over where a transmission of a transcendental
# expansion begins. A beacon sending pi sends it from the start, so offset 0
# carries most of the mass and offset k gets 1/(k+1)**2 of it; the sum
# converges (to zeta(2)), which is the property that matters. Under a flat
# count the correction is proportional to the stored depth, so every p-value
# depends on an arbitrary constant and searching deeper makes every existing
# hit look worse. Under a convergent prior the expansion can be searched as
# far as anyone likes, and a hit found 74 bits in is charged for having been
# found 74 bits in.
SEQ_STREAM_OFFSET_EXPONENT = 2
_ZETA_2 = math.pi**2 / 6

# Half the prior mass on the numeric catalogs and half on the expansions. It
# doubles every numeric p-value against the old flat bound, which moves no
# candidate anywhere and is the price of saying where the mass went.
SEQ_CORPUS_SPLIT = 2.0
# A set bit is an advance, so a word's steps are negative where the bit is
# set -- but which of the two symbols the transmitter chose to kick on is a
# convention, and a code that kicks on the other one is the same catalog word
# with every bit flipped. The matched filter therefore tries each window
# against both polarities of the folded means and charges every trial count
# this factor. Measured before it was added (2026-09-22, 12,054 stars at
# score >= 0.5): the flipped catalogs named stars at the same rate as the
# upright ones, so the factor is a real cost paid for a real freedom, not a
# bonus. Reading the word backwards was tested too and has no physical
# reading (the fold already searches every rotation), so it is not tried.
SEQ_POLARITIES = 2.0


N_STREAMS = sum(1 for seq in KNOWN_BIT_SEQ if not seq["values"])


def stream_trials(start):
    """Alignments charged to a transcendental hit whose match begins at `start`.

    `N_STREAMS` because the expansions are searched in parallel and any of
    them could have produced the hit; `zeta(2)` normalises the offset prior
    over an unbounded expansion. Takes a scalar or an array of offsets.
    """
    k = np.maximum(0, np.asarray(start, dtype=float))
    out = N_STREAMS * _ZETA_2 * (k + 1) ** SEQ_STREAM_OFFSET_EXPONENT
    return float(out) if out.ndim == 0 else out


def _bit_confidence(
    ys: Sequence[float], ones: int, n: int, mid: float, global_score: float, n_rep: float
) -> float:
    if not n:
        return 0.0
    vote = abs(2 * (ones / n) - 1)
    mean = sum(ys) / n
    var_sum = sum((y - mean) ** 2 for y in ys)
    sd = math.sqrt(var_sum / max(1, n - 1))
    sep = abs(mean - mid) / (sd + abs(mid) * 0.05 + 1e-9)
    snr = 1 / (1 + math.exp(-(sep - 0.9) * 1.6))
    structure = _clamp01((global_score - 0.5) / 0.45)
    rep = 1 - math.exp(-n_rep / 2.2)
    return _clamp01(
        0.06 + 0.94 * vote * (0.2 + 0.8 * structure) * (0.4 + 0.6 * snr) * (0.45 + 0.55 * rep)
    )


def _binary_entropy(p: float) -> float:
    if not (p > 0) or not (p < 1):
        return 0.0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


def _score_bit_period(pts: Sequence[dict], bits: Sequence[int], mid: float, L: int) -> tuple:
    acc = [0] * L
    cnt = [0] * L
    ys = [[] for _ in range(L)]
    for i, bit in enumerate(bits):
        k = int(((pts[i]["c"] % L) + L) % L)
        acc[k] += bit
        cnt[k] += 1
        ys[k].append(pts[i]["y"])
    if any(n < 1 for n in cnt):
        return None
    word = [1 if acc[i] / cnt[i] > 0.5 else 0 for i in range(L)]
    agree = 0
    for i, bit in enumerate(bits):
        k = int(((pts[i]["c"] % L) + L) % L)
        if bit == word[k]:
            agree += 1
    score = agree / len(bits)
    n_rep = len(bits) / L
    conf = [_bit_confidence(ys[k], acc[k], cnt[k], mid, score, n_rep) for k in range(L)]
    info = [{"p": conf[k], "n": cnt[k], "ones": acc[k]} for k in range(L)]
    residual = 0.0
    for k in range(L):
        e = min(acc[k], cnt[k] - acc[k]) / cnt[k]
        residual += cnt[k] * _binary_entropy(e)
    return {
        "bits": "".join(str(b) for b in word),
        "score": score,
        "n": L,
        "conf": conf,
        "info": info,
        "mdl": L + residual,
    }


def split_means(means: Sequence[float]) -> float | None:
    """Where to cut L folded means into two levels: the split that explains
    the most of their variance (Otsu), with at least one mean on each side.

    `two_level_mid`'s 20/80 rule is right for a long series of noisy points
    but wrong for a handful of folded means: a word with one or two set bits
    in eight has both percentiles on the same level and no split at all,
    which is why lopsided words used to be undecodable. None when the means
    have no spread.
    """
    xs = sorted(float(m) for m in means)
    n = len(xs)
    if n < 2 or not xs[-1] > xs[0]:
        return None
    total = sum(xs)
    best = -1.0
    cut = None
    run = 0.0
    for i in range(1, n):
        run += xs[i - 1]
        m0 = run / i
        m1 = (total - run) / (n - i)
        between = (i / n) * (1 - i / n) * (m0 - m1) ** 2
        if between > best:
            best = between
            cut = (xs[i - 1] + xs[i]) / 2.0
    return cut


def fold_means(pts: Sequence[dict], L: int, k: int = 1) -> tuple[list[float], list[int]] | None:
    """Mean of the series in each of the L residues of `(cycle // k) % L`.

    This is where repetition becomes signal-to-noise: a word repeated r times
    has each bit measured r times, and averaging those before thresholding is
    worth sqrt(r), which thresholding each cycle first and voting afterwards
    throws away when the per-cycle signal is below the noise. None when some
    residue was never observed, because a word with an unmeasured bit is not
    a word.
    """
    acc = [0.0] * L
    cnt = [0] * L
    for pt in pts:
        j = _residue(pt["c"], L, k)
        acc[j] += pt["y"]
        cnt[j] += 1
    if any(n < 1 for n in cnt):
        return None
    return [a / n for a, n in zip(acc, cnt)], cnt


def _fold_word(pts: Sequence[dict], L: int, k: int = 1) -> dict | None:
    """The two-level word of length L at k cycles per bit, and what it costs.

    The word is the folded means split at `two_level_mid`; the cost is a
    description length in amplitude space -- n ln(RSS/n) for the residuals
    about the two levels, ln 2 per bit of the word and ln n per level -- so a
    longer word has to buy its extra bits with a better fit. `sep_z` is how
    many standard errors apart the two levels are.
    """
    fm = fold_means(pts, L, k)
    if fm is None:
        return None
    means, cnt = fm
    mid = split_means(means)
    if mid is None:
        return None
    word = [1 if m < mid else 0 for m in means]
    if not 0 < sum(word) < L:
        return None
    groups: dict[int, list[float]] = {0: [], 1: []}
    for pt in pts:
        groups[word[_residue(pt["c"], L, k)]].append(pt["y"])
    level = {b: sum(v) / len(v) for b, v in groups.items()}
    n = len(pts)
    rss = 0.0
    for pt in pts:
        r = pt["y"] - level[word[_residue(pt["c"], L, k)]]
        rss += r * r
    # A bit is a binary choice and costs ln 2; the two levels are real
    # parameters and cost ln n each. Charging every bit ln n, as a BIC would,
    # makes a long word unable to pay for itself against a short one at
    # moderate noise, which is the case this decoder exists for.
    cost = n * math.log(max(rss, 1e-300) / n) + L * math.log(2.0) + 2 * math.log(n)
    sd = math.sqrt(rss / max(n - 2, 1))
    n0, n1 = len(groups[0]), len(groups[1])
    sep = (level[0] - level[1]) / (sd * math.sqrt(1.0 / n0 + 1.0 / n1)) if sd > 0 else float("inf")
    return {
        "word": word,
        "means": means,
        "counts": cnt,
        "cost": cost,
        "sep_z": sep,
        "levels": (level[0], level[1]),
        "n": L,
        "k": k,
    }


def _fold_level_word(pts: Sequence[dict], L: int, k: int = 1) -> dict | None:
    """The `fold` model's word of length L at k cycles per bit, and its cost.

    The levels are averaged per residue, the bits are read off the cyclic
    first difference of those means (bit j is set where mean j sits below
    mean j - 1, Otsu's cut deciding what "below" is), and the cost is the
    description length of the *levels* about the staircase the word
    predicts: a repeated one-sided word walks down by dt per set bit and the
    ramp of its duty cycle is what stage 1's detrend removed, so residue j
    sits at a - s (S_j - W j / L) with S_j the set bits through j and W the
    set bits in the word. n ln(RSS/n) for the residuals about that model, ln 2
    per bit and ln n for each of the two real parameters -- the same currency
    `_fold_word` prices the step model in, so the two are comparable.
    """
    fm = fold_means(pts, L, k)
    if fm is None:
        return None
    means, cnt = fm
    diffs = [means[j] - means[j - 1] for j in range(L)]
    cut = split_means(diffs)
    if cut is None:
        return None
    word = [1 if d < cut else 0 for d in diffs]
    w_total = sum(word)
    if not 0 < w_total < L:
        return None
    run = 0
    u = []
    for j in range(L):
        run += word[j]
        u.append(run - w_total * j / L)
    # Weighted least squares of the residue means on the model shape.
    n = len(pts)
    sw = float(n)
    su = sum(c * x for c, x in zip(cnt, u))
    sm = sum(c * m for c, m in zip(cnt, means))
    suu = sum(c * x * x for c, x in zip(cnt, u))
    sum_ = sum(c * x * m for c, x, m in zip(cnt, u, means))
    det = sw * suu - su * su
    if not det > 0:
        return None
    # level_j = a - s u_j
    s = -(sw * sum_ - su * sm) / det
    a = (sm + s * su) / sw
    level = [a - s * x for x in u]
    rss = 0.0
    for pt in pts:
        r = pt["y"] - level[_residue(pt["c"], L, k)]
        rss += r * r
    cost = n * math.log(max(rss, 1e-300) / n) + L * math.log(2.0) + 2 * math.log(n)
    sd = math.sqrt(rss / max(n - 2, 1))
    n1 = sum(c for c, b in zip(cnt, word) if b)
    n0 = n - n1
    sep = s / (sd * math.sqrt(1.0 / max(n0, 1) + 1.0 / max(n1, 1))) if sd > 0 else float("inf")
    return {
        "word": word,
        "means": diffs,
        "level_means": means,
        "counts": cnt,
        "cost": cost,
        "sep_z": sep,
        "levels": (0.0, -s),
        "n": L,
        "k": k,
    }


def pseudo_steps(pts: Sequence[dict], L: int, k: int = 1) -> list[dict]:
    """The step into each observed cycle, read against the nearest observation
    of the previous bit.

    A level series with no consecutive cycles has no steps to difference, but
    once a word length and bit period are chosen every point knows which
    residue holds the bit before it. Each point is paired with the point of
    residue j - 1 nearest to it in cycle, and the difference is the step; its
    fold means are the cyclic differences of the level means, its variance is
    twice the level noise, and it carries the original cycle, so everything
    written for a step series -- the squareness test, the matched filter, the
    per-cycle votes -- reads it as one.
    """
    by_res: dict[int, list[tuple[int, float]]] = {}
    for pt in pts:
        by_res.setdefault(_residue(pt["c"], L, k), []).append((pt["c"], pt["y"]))
    for lst in by_res.values():
        lst.sort()
    out = []
    for pt in pts:
        prev = by_res.get((_residue(pt["c"], L, k) - 1) % L)
        if not prev:
            continue
        cs = [c for c, _ in prev]
        i = bisect_left(cs, pt["c"])
        cands = [x for x in (i - 1, i) if 0 <= x < len(prev)]
        _c_nn, y_nn = min((prev[x] for x in cands), key=lambda cy: abs(cy[0] - pt["c"]))
        out.append({"y": pt["y"] - y_nn, "c": pt["c"]})
    return out


def _describe_word(pts: Sequence[dict], fold: dict, mid: float) -> dict:
    """Per-cycle agreement, confidence and vote counts for a folded word."""
    L, k, word = fold["n"], fold["k"], fold["word"]
    ys: list[list[float]] = [[] for _ in range(L)]
    ones = [0] * L
    agree = 0
    for pt in pts:
        j = _residue(pt["c"], L, k)
        bit = 1 if pt["y"] < mid else 0
        ones[j] += bit
        ys[j].append(pt["y"])
        if bit == word[j]:
            agree += 1
    score = agree / len(pts)
    n_rep = len(pts) / L
    conf = [_bit_confidence(ys[j], ones[j], len(ys[j]), mid, score, n_rep) for j in range(L)]
    info = [{"p": conf[j], "n": len(ys[j]), "ones": ones[j]} for j in range(L)]
    return {"bits": "".join(str(b) for b in word), "score": score, "conf": conf, "info": info}


def _trim_repeating_bits(
    bit_str: str, score: float | None, conf: Sequence[float] | None = None, info: dict | None = None
) -> dict | None:
    s = "".join(ch for ch in str(bit_str) if ch in "01")
    if not s:
        return None
    if conf and len(conf) == len(s):
        c = list(conf)
    else:
        fill = 1 if score is None else score
        c = [fill] * len(s)
    inf = list(info) if info and len(info) == len(s) else None
    for p in range(2, len(s) // 2 + 1):
        if len(s) % p:
            continue
        word = s[:p]
        if not ("0" in word and "1" in word):
            continue
        if s == word * (len(s) // p):
            folded = []
            folded_info = [] if inf else None
            for i in range(p):
                sm = 0.0
                n = 0
                for j in range(i, len(c), p):
                    sm += c[j]
                    n += 1
                folded.append(sm / n)
                if folded_info is not None:
                    votes = ones = psum = k = 0
                    for j in range(i, len(inf), p):
                        votes += inf[j].get("n") or 0
                        ones += inf[j].get("ones") or 0
                        psum += inf[j].get("p") or 0
                        k += 1
                    folded_info.append(
                        {"p": psum / k if k else folded[i], "n": votes, "ones": ones}
                    )
            return {
                "bits": word,
                "score": 1 if score is None else score,
                "n": p,
                "conf": folded,
                "info": folded_info,
            }
    return {"bits": s, "score": 1 if score is None else score, "n": len(s), "conf": c, "info": inf}


def oc_points(oc_s: Sequence[float | None], cycle: Sequence[int] | None) -> list[dict]:
    """Finite (cycle, value) points in cycle order; index order without a cycle column."""
    pts = []
    for i, y in enumerate(oc_s):
        if y is None or not math.isfinite(y):
            continue
        c = cycle[i] if cycle and i < len(cycle) and cycle[i] is not None else i
        pts.append({"y": y, "c": c})
    pts.sort(key=lambda p: p["c"])
    return pts


def increments(pts: Sequence[dict]) -> list[dict]:
    """The first difference of the O-C, one point per observed step.

    Point `c` carries `y(c) - y(c - 1)`: on a staircase the bit of cycle c is
    the step *into* c, so the increment sits on the later cycle. A step across
    a gap is the sum of the bits nobody observed plus whatever the gap did to
    the ephemeris -- a sector boundary's level jump, a season's worth of
    unseen bits -- and it is dropped rather than averaged: an average lands
    a mid-level value in one residue of the fold, and it is what stitching
    sectors used to cost. Duplicate cycle indices contribute nothing. Public
    because stage 1 calibrates the decode against a shuffle null and has to
    difference the draws exactly the way the observation is differenced.
    """
    out = []
    for prev, cur in pairwise(pts):
        if cur["c"] - prev["c"] != 1:
            continue
        out.append({"y": cur["y"] - prev["y"], "c": cur["c"]})
    return out


# Steps further than this many robust sigmas from the median are pulled in to
# that distance before the fold. One glitch, sector boundary or bad cycle
# otherwise dominates the step variance, and the cheapest two-level
# description of such a series is "that one residue is different": a long
# word of zeros with a single one, at whatever length isolates the outlier.
# Three sigma, not five: a glitch pulled in to five sigma is still the most
# extreme step in the series, and three of them on one residue of some fold
# still bought a ten-bit word; at three sigma a clipped step is one ordinary
# tail point among the 0.3% Gaussian noise puts there anyway.
#
# An outlier is rare by definition. If more than STEP_CLIP_MAX_FRACTION of
# the steps lie beyond the limit, they are not outliers but the series' own
# structure -- a code's set bits, or the transitions of a square wave -- and
# nothing is clipped: pulling a fifth of the steps to one common value would
# manufacture exactly the two-level structure the fold is looking for, in the
# observation and in every shuffle draw whose robust scale happened to
# collapse. A tenth is generous to glitches (one bad cycle is two outlying
# steps, so this is six glitches in 120 cycles) and still below the sparsest
# catalog packing (powers of two in 8-bit is one set bit in eight).
STEP_CLIP_SIGMA = 3.0
STEP_CLIP_MAX_FRACTION = 0.1

# The other way outliers can be structure rather than glitches: a code with a
# duty cycle below STEP_CLIP_MAX_FRACTION. A transmitter with no energy store
# that saves for N cycles and kicks once puts N times the per-cycle amplitude
# into one step, so its set bits are exactly the rare, large, *one-sided*
# steps the clip was built to flatten. A glitch has no sign preference and no
# characteristic size; a bad cycle is two outlying steps of opposite sign. So
# when at least KICK_MIN_STEPS outlying steps all lie on the same side of the
# median and agree in size to within half their median distance
# (KICK_SPREAD), they are kept as they are. Half, not all of it: two 9-sigma
# kicks and one ordinary 3.2-sigma tail point are not a cluster, and at unity
# they were. The rule is applied to every shuffle draw too
# (stage1_pulse._clip_steps_nb), so the null it changes is the null of the
# same statistic.
KICK_MIN_STEPS = 3
KICK_SPREAD = 0.5


def one_sided_cluster(dist: np.ndarray, sign: np.ndarray) -> bool:
    """Are these outlying steps a one-sided cluster of kicks, not glitches?

    `dist` is each outlier's distance from the median and `sign` which side it
    lies on. True when there are KICK_MIN_STEPS or more, every one on the same
    side, and the spread of their distances (largest minus smallest) is no
    more than KICK_SPREAD times their median distance.
    """
    n = dist.size
    if n < KICK_MIN_STEPS:
        return False
    if np.any(sign != sign[0]):
        return False
    med = float(np.median(dist))
    return bool(med > 0 and (float(dist.max()) - float(dist.min())) <= KICK_SPREAD * med)


def clip_steps(values) -> np.ndarray:
    """Winsorize a step series at STEP_CLIP_SIGMA robust sigmas of the median.

    The scale is the normal-equivalent of the interquartile range. The clip
    is applied only when at most STEP_CLIP_MAX_FRACTION of the steps lie
    beyond the limit and those steps are not a one-sided cluster of the same
    size (`one_sided_cluster`, a sparse code's kicks); with fewer than eight
    steps, or a zero range, or more than that many outlying steps, the series
    comes back unchanged. Public because stage 1 differences its shuffle draws
    the way the observation is differenced, and has to clip them the same way
    too.
    """
    y = np.asarray(values, dtype=float)
    if y.size < 8:
        return y
    q1, med, q3 = np.percentile(y, [25, 50, 75])
    scale = (q3 - q1) / 1.349
    if not scale > 0:
        return y
    lim = STEP_CLIP_SIGMA * scale
    dev = y - med
    beyond = np.abs(dev) > lim
    if not beyond.any() or beyond.mean() > STEP_CLIP_MAX_FRACTION:
        return y
    if one_sided_cluster(np.abs(dev[beyond]), np.sign(dev[beyond])):
        return y
    return np.clip(y, med - lim, med + lim)


# ---- kicks: a message written in the spacing of rare pulses -----------------
# The catalogs as integer sequences, for a code whose set bits are rare kicks
# and whose message is where they fall. Two encodings are read. `positions`:
# the kicks land on the catalog's values as cycle numbers from some origin,
# so the gaps between them are the catalog's first differences (kicking on
# every prime-numbered cycle gives 1, 2, 2, 4, 2, 4, ...). `gaps`: the gaps
# themselves are the values (a kick, then 2 cycles, a kick, then 3, then 5,
# ...). Deeper than the packed-binary catalogs because a gap is a cycle count
# and can run to the length of the series.
KICK_GAP_MAX = 4096
KICK_ENCODINGS = ("positions", "gaps")


def _kick_value_catalogs(max_value: int = KICK_GAP_MAX) -> list[tuple[str, list[int]]]:
    primes = _primes_to(max_value)
    fib = [v for v in _fib_list(60) if v <= max_value]
    count = list(range(1, max_value + 1))
    pow2 = [1 << i for i in range(64) if (1 << i) <= max_value]
    squares = [v * v for v in range(1, max_value) if v * v <= max_value]
    tri = [v * (v + 1) // 2 for v in range(1, max_value) if v * (v + 1) // 2 <= max_value]
    return [
        ("primes", primes),
        ("Fibonacci", fib),
        ("counting", count),
        ("powers of two", pow2),
        ("squares", squares),
        ("triangular", tri),
    ]


KICK_VALUE_CATALOGS = _kick_value_catalogs()


def kick_gap_sequences() -> list[tuple[str, str, list[int]]]:
    """Every (catalog, encoding, gap sequence) a run of kick gaps is read
    against: the catalog's first differences for `positions`, its values for
    `gaps`. The order is the search order, and stage 1's compiled null walks
    the same list, so a tie in run length lands on the same entry."""
    out = []
    for name, vals in KICK_VALUE_CATALOGS:
        out.append((name, "positions", [b - a for a, b in pairwise(vals)]))
        out.append((name, "gaps", list(vals)))
    return out


KICK_GAP_SEQUENCES = kick_gap_sequences()

# Fewer matching gaps than this names nothing: two gaps agree with some window
# of the counting numbers or the primes far too often to mean anything.
KICK_MIN_RUN = 3


def kick_gap_match(gaps: Sequence[int]) -> dict:
    """The longest run of consecutive kick gaps that reads as a window of one
    of the integer catalogs, in either encoding.

    `gaps` are the cycle counts between consecutive kicks, in order; a gap that
    could hide an unobserved kick should be passed as a non-positive value,
    which matches nothing and so ends any run. Returns the run length, the
    catalog and encoding, where in the catalog the run begins, where in the
    gap sequence it begins, and a label -- or run 0 and no label when nothing
    of KICK_MIN_RUN gaps matches. The catalog's own trials are not charged
    here: stage 1 calibrates the run against random kick sets of the same size
    on the same observed cycles, which searches exactly this.
    """
    empty = {"run": 0, "catalog": None, "encoding": None, "start": None, "at": None, "label": None}
    best = dict(empty)
    g = [int(v) for v in gaps]
    if len(g) < KICK_MIN_RUN:
        return best
    for name, enc, seq in KICK_GAP_SEQUENCES:
        pos: dict[int, list[int]] = {}
        for j, v in enumerate(seq):
            pos.setdefault(v, []).append(j)
        for i in range(len(g)):
            if g[i] <= 0:
                continue
            for j in pos.get(g[i], ()):
                r = 0
                while i + r < len(g) and j + r < len(seq) and g[i + r] == seq[j + r]:
                    r += 1
                if r > best["run"]:
                    best.update({"run": r, "catalog": name, "encoding": enc, "start": j, "at": i})
    if best["run"] < KICK_MIN_RUN:
        return empty
    name, enc, j, r = best["catalog"], best["encoding"], best["start"], best["run"]
    vals = dict(KICK_VALUE_CATALOGS)[name]
    if enc == "positions":
        # r gaps span r + 1 catalog values, from value j.
        shown = vals[j : j + min(r + 1, 7)]
        tail = ", ..." if r + 1 > 7 else ""
        best["label"] = f"kicks at {name} [{', '.join(str(v) for v in shown)}{tail}]"
    else:
        shown = vals[j : j + min(r, 6)]
        tail = ", ..." if r > 6 else ""
        best["label"] = f"{name} as kick gaps [{', '.join(str(v) for v in shown)}{tail}]"
    return best


def code_points(
    oc_s: Sequence[float | None], cycle: Sequence[int] | None, model: str = DEFAULT_CODE_MODEL
) -> list[dict]:
    """The series the bits are thresholded from, under either model.

    Under the staircase model the steps are winsorized (`clip_steps`) so a
    single outlying step cannot buy a word of its own.
    """
    if model not in CODE_MODELS:
        raise ValueError(f"unknown code model {model!r}; one of {CODE_MODELS}")
    pts = oc_points(oc_s, cycle)
    if model != "staircase":
        return pts
    inc = increments(pts)
    clipped = clip_steps([p["y"] for p in inc])
    return [{"y": float(v), "c": p["c"]} for v, p in zip(clipped, inc)]


def decode_bit_code(
    oc_s: Sequence[float | None],
    cycle: Sequence[int] | None,
    model: str = DEFAULT_CODE_MODEL,
    fold_at: tuple[int, int] | None = None,
) -> dict | None:
    """Read the O-C as a clock-locked bit word.

    Under the default `staircase` model the bits are the first difference of
    the O-C (see the module docstring); under `square` they are its level;
    under `fold` the levels are folded first and the bits are the steps
    between neighbouring residue means (`_fold_level_word`), for a series
    that has no two consecutive cycles to difference.

    Fold first, threshold second: for every word length L and every bit
    period k in CODE_BIT_CYCLES the series is averaged per residue of
    `(cycle // k) % L` and the *means* are split into two levels, so a word
    repeated r times is read at sqrt(r) the per-cycle signal-to-noise. The
    word that describes the series most cheaply (amplitude-space description
    length, `_fold_word`) wins; the per-cycle agreement, the squareness test
    against a sinusoid and the catalog match then run on it.
    """
    pts = code_points(oc_s, cycle, model)
    if len(pts) < 8:
        return None
    mid = two_level_mid([p["y"] for p in pts])
    if mid is None:
        return None
    fold_word = _fold_level_word if model == "fold" else _fold_word
    best = None
    # `fold_at` = (L, k) skips the search: stage 1 has already run it in
    # compiled form and only wants this word described and named.
    for k in bit_cycles_for(len(pts)) if fold_at is None else (int(fold_at[1]),):
        max_l = min(SEQ_MAX_BITS, len(pts) // (2 * k))
        if max_l < CODE_MIN_WORD_BITS:
            continue
        lengths = range(CODE_MIN_WORD_BITS, max_l + 1) if fold_at is None else (int(fold_at[0]),)
        for L in lengths:
            if L > max_l:
                continue
            fold = fold_word(pts, L, k)
            if fold is None:
                continue
            if (
                best is None
                or fold["cost"] < best["cost"] - 1e-9
                or (abs(fold["cost"] - best["cost"]) <= 1e-9 and L * k < best["n"] * best["k"])
            ):
                best = fold
    if best is None:
        return None
    if model == "fold":
        # Everything below reads steps. The level series has none, so the
        # steps are read against the nearest observation of the previous bit
        # at the fold that won (pseudo_steps); the fold's own word is kept.
        pts = pseudo_steps(pts, best["n"], best["k"])
        if len(pts) < 8:
            return None
        mid = two_level_mid([p["y"] for p in pts])
        if mid is None:
            return None
    desc = _describe_word(pts, best, mid)
    trimmed = _trim_repeating_bits(desc["bits"], desc["score"], desc["conf"], desc["info"])
    if not trimmed:
        return None
    trimmed["mdl"] = best["cost"]
    trimmed["model"] = model
    trimmed["k"] = best["k"]
    trimmed["fold_means"] = best["means"]
    trimmed["fold_counts"] = best["counts"]
    trimmed["sep_z"] = best["sep_z"]
    trimmed["levels"] = list(best["levels"])
    # Two levels against a sinusoid, in whichever series carries the bits: a
    # smooth O-C differenced is still smooth, so the test transfers as it is.
    sq = _squareness(pts, trimmed["bits"], trimmed["n"], best["k"])
    trimmed["squareness"] = None if sq is None else sq["squareness"]
    trimmed["squareness_at_period"] = None if sq is None else sq["squareness_at_period"]
    trimmed["beat_frequency"] = None if sq is None else sq["beat_frequency"]
    # A rotation word is what an irrational period ratio writes, not a message.
    trimmed["rotation_word"] = is_rotation_word(trimmed["bits"])
    # Name the word: the matched filter over every catalog window, at this
    # fold's length and bit period, tolerant of a bit or two. It reads the
    # steps themselves, not the thresholded bits, so a word whose bits are
    # individually below the noise but collectively right still gets a name.
    import matched  # matched imports this module, so it is bound here

    hit = matched.matched_catalog(
        [pt["y"] for pt in pts], [pt["c"] for pt in pts], trimmed["n"], best["k"]
    )
    trimmed["recog"] = recognize_hit(hit, trimmed["n"])
    trimmed["label"] = trimmed["recog"].get("label")
    # The two-level test on the window itself, which a fixed pattern earns at
    # the two repeats the filter needs rather than the three the fold's own
    # word does (window_squareness). Rotation-free, like the squareness above.
    named = bool(trimmed["label"]) and trimmed["label"] != "unrecognized"
    trimmed["window_squareness"] = window_squareness(pts, hit, best["k"]) if named else None
    rotation = trimmed["recog"].get("rotation")
    if rotation:
        _rotate_word(trimmed, rotation)
    _mark_near_misses(trimmed)
    return trimmed


def _mark_near_misses(code: dict) -> None:
    """Which of the recognised values contain a bit the word disagrees on.

    The viewer draws one pip per catalog value (or per four bits of an
    expansion) and dims the ones a near miss falls in. Each span gets `off`,
    and the recognition record an `off_mask` with one character per pip, so
    the pips can be drawn from the stored row without the word.
    """
    recog = code.get("recog") or {}
    window = recog.get("window")
    bits = code.get("bits") or ""
    if not recog.get("label") or recog["label"] == "unrecognized" or not window:
        return
    n = min(len(bits), len(window))

    def off(lo: int, hi: int) -> bool:
        return any(bits[i] != window[i] for i in range(max(0, lo), min(hi, n)))

    if recog.get("encoding") == "bitstream":
        pips = seq_identified_count(recog, bits)
        mask = "".join("1" if off(4 * g, 4 * g + 4) else "0" for g in range(pips))
    else:
        spans = recog.get("spans") or []
        for sp in spans:
            sp["off"] = off(int(sp["lo"]), int(sp["hi"]))
        by_value = {}
        for sp in spans:
            by_value[sp["value"]] = by_value.get(sp["value"], False) or sp["off"]
        # A value in a later repeat of a tiled window has no span; it is
        # never marked, which errs on the side of the exact reading.
        mask = "".join("1" if by_value.get(v, False) else "0" for v in recog.get("nums") or [])
    recog["off_mask"] = mask


def _rotate_word(code: dict, rotation: int) -> None:
    """Turn the decoded word so it reads in catalog order.

    The fold numbers its residues from cycle 0, which is where the observation
    happened to start, not where the transmitter's word does; the matched
    filter searches every rotation and says which one lines up with the
    catalog window. Window bit i is fold residue (i - rotation) mod L, so the
    per-bit arrays roll forward by `rotation`. The squareness and agreement
    score are rotation-free and are left alone.
    """
    n = len(code["bits"])
    if n < 2:
        return
    r = int(rotation) % n
    if r == 0:
        return

    def roll(seq):
        items = list(seq)
        k = r % len(items) if items else 0
        return items[-k:] + items[:-k] if k else items

    code["bits"] = "".join(roll(code["bits"]))
    for key in ("conf", "info", "fold_means", "fold_counts"):
        if code.get(key):
            code[key] = roll(code[key])
    code["rotation"] = r


def recognize_hit(hit: dict | None, word_len: int) -> dict:
    """The viewer's recognition record for a matched-filter hit.

    A hit is only recorded while its p-value is below one: past that the
    label is arithmetic, not a match (at 4,096 bits of expansion almost every
    12-bit word appears somewhere in pi). For a numeric packing the record
    carries the catalog values the window covers, `u` for a value cut off at
    either end, and the span of each value in the rotated word so the viewer
    can group its bits; for an expansion it carries where the match begins.
    `mismatch` is how many of the window's bits the decoded word disagrees
    with, and `window` the catalog bits themselves, so the viewer can mark
    the bits that differ.
    """
    if not hit or not hit.get("name") or float(hit.get("p", 1.0)) >= 1.0:
        return {"label": "unrecognized"}
    name = str(hit["name"])
    enc = str(hit.get("encoding") or "")
    start = int(hit.get("start") or 0)
    L = int(hit.get("len") or word_len)
    inverted = int(hit.get("polarity") or 1) < 0
    window = hit.get("window")
    if inverted and window:
        # The catalog window matched the word with every bit flipped, so in
        # the word's own convention (a set bit is a negative step) the
        # pattern to draw against it is the window's complement.
        window = "".join("1" if ch == "0" else "0" for ch in str(window))
    out = {
        "encoding": enc,
        "start": start,
        "p": float(hit["p"]),
        "p_amplitude": hit.get("p_amplitude"),
        "p_bits": hit.get("p_bits"),
        "rho": hit.get("rho"),
        "z": hit.get("z"),
        "gof": hit.get("gof"),
        "rotation": int(hit.get("rotation") or 0),
        "len": L,
        "mismatch": int(hit.get("mismatch") or 0),
        "window": window,
        "corrected": bool(hit.get("corrected")),
        "inverted": inverted,
    }
    suffix = " (inverted)" if inverted else ""
    if enc == "bitstream":
        out["label"] = (name if start == 0 else f"{name} from bit {start + 1}") + suffix
        out["n"] = L
        return out
    encoder = dict(catalog_encoders()).get(enc)
    values = (KNOWN_BIT_SEQ[int(hit["catalog"])]["values"] or []) if "catalog" in hit else []
    pos = 0
    complete: list[int] = []
    spans: list[dict] = []
    partial_start = False
    partial_end = False
    win_end = start + L
    for v in values:
        b = encoder(v) if encoder else None
        if not b:
            continue
        a, e = pos, pos + len(b)
        pos = e
        if e <= start or a >= win_end:
            continue
        if a < start < e:
            partial_start = True
        if a < win_end < e:
            partial_end = True
        if a >= start and e <= win_end:
            complete.append(v)
            # Spans are in word positions; a tiled window repeats the word,
            # so a value past the first repeat is clipped to it.
            lo, hi = a - start, e - start
            if lo < word_len:
                spans.append({"value": v, "lo": lo, "hi": min(hi, word_len)})
    # A window that covers no whole value -- two fragments of 8-bit numbers,
    # say -- names nothing the viewer can show a digit for, and a fifth of the
    # short numeric windows are like that. It is not a recognition.
    if not complete:
        return {"label": "unrecognized"}
    out.update(
        {
            "label": _format_seq_label(name, complete, partial_start, partial_end) + suffix,
            "nums": complete,
            "spans": spans,
            "n": len(complete),
            "partialStart": partial_start,
            "partialEnd": partial_end,
        }
    )
    return out


def seq_family(label: str | None) -> str:
    lab = str(label or "")
    if not lab or lab == "unrecognized":
        return ""
    cut = lab.find(" [")
    if cut > 0:
        return lab[:cut]
    return lab.split(" ")[0] if lab else ""


def seq_identified_count(recog: dict | None, bits: str | None) -> int:
    if not recog:
        return 0
    if recog.get("encoding") == "bitstream":
        n = recog.get("n")
        if n is None:
            n = len("".join(ch for ch in str(bits or "") if ch in "01"))
        return n // 4
    nums = recog.get("nums")
    if nums:
        return len(nums)
    lab = str(recog.get("label") or "")
    m = lab[lab.find("[") + 1 : lab.find("]")] if "[" in lab and "]" in lab else ""
    if not m:
        return 0
    return len([x.strip() for x in m.split(",") if x.strip() and x.strip() != "u"])


def seq_u_count(recog: dict | None) -> int:
    if not recog or recog.get("encoding") == "bitstream":
        return 0
    u = 0
    if recog.get("partialStart"):
        u += 1
    if recog.get("partialEnd"):
        u += 1
    if u:
        return u
    lab = str(recog.get("label") or "")
    m = lab[lab.find("[") + 1 : lab.find("]")] if "[" in lab and "]" in lab else ""
    if not m:
        return 0
    return len([x.strip() for x in m.split(",") if x.strip() == "u"])


def seq_sort_n(recog: dict | None, bits: str | None) -> float:
    if not recog:
        return 0.0
    return seq_identified_count(recog, bits) + 0.49 * seq_u_count(recog)


def fields_from_code(code: dict | None, spark: str) -> dict:
    """The viewer's sequence fields for an already-decoded word, so a caller
    that has run `decode_bit_code` need not run it again."""
    return _fields_from_code(code, spark)


def _fields_from_code(code: dict | None, spark: str) -> dict:
    out = dict(EMPTY)
    out["code"] = code or None
    out["oc_spark"] = spark or ""
    bits = "".join(ch for ch in str((code or {}).get("bits") or "") if ch in "01")
    L = len(bits)
    out["seq_bit_len"] = L
    # Recorded even when the word matches no catalog sequence: squareness is
    # about the waveform, not about recognising the digits.
    out["seq_squareness"] = (code or {}).get("squareness")
    # The same predicate stage 2 uses to tell a phase code from a smooth wave.
    # Stage 1's triage gate reads it, so the two cannot drift apart.
    out["seq_square"] = is_square_code(code)
    # The catalog route's guard: the same test asked of the named window.
    out["seq_window_squareness"] = (code or {}).get("window_squareness")
    out["seq_window_square"] = is_square_match(code)
    # Why a word was refused, for the viewer: the period-L squareness the
    # test used to stop at, the frequency of the sinusoid that explains the
    # steps, and whether the word is a rotation coding (a two-mode beat).
    out["seq_squareness_at_period"] = (code or {}).get("squareness_at_period")
    out["seq_beat_frequency"] = (code or {}).get("beat_frequency")
    out["seq_rotation_word"] = bool((code or {}).get("rotation_word"))
    if not code or not code.get("label") or code.get("label") == "unrecognized":
        out["seq_p"] = 1.0
        return out
    recog = code.get("recog") or code
    out["seq_matched"] = True
    out["seq_label"] = code["label"]
    out["seq_name"] = seq_family(code["label"])
    out["seq_digit_count"] = seq_identified_count(recog, bits)
    out["seq_sort_n"] = seq_sort_n(recog, bits)
    out["seq_p"] = float(recog.get("p") or 1.0)
    out["seq_off_mask"] = str(recog.get("off_mask") or "")
    return out


def downsample_pairs(pairs: list[tuple[Any, Any]]) -> list[tuple[Any, Any]]:
    n = len(pairs)
    if n <= MAX_OC:
        return pairs
    step = math.ceil(n / MAX_OC)
    return pairs[::step]


def score_oc_pairs(pairs: list[tuple[Any, Any]]) -> dict:
    """pairs: (cycle_index, oc_seconds) in cycle order."""
    if not pairs:
        return dict(EMPTY)
    # The spark is a drawing and can be thinned; the decode cannot. Thinning
    # to every k-th cycle leaves no consecutive pair, and the staircase
    # decoder reads its bits from consecutive steps -- a 3,000-cycle RR Lyrae
    # would decode to nothing at all. So the word is read from every cycle.
    shown = downsample_pairs(pairs)
    spark = spark_path([v for _, v in shown])
    cycle = [c for c, _ in pairs]
    oc_s = [v for _, v in pairs]
    try:
        code = decode_bit_code(oc_s, cycle)
        return _fields_from_code(code, spark)
    except Exception:
        out = dict(EMPTY)
        out["oc_spark"] = spark
        return out


def oc_seconds(det: object, shift: object) -> float | None:
    raw = det if det is not None else shift
    if raw is None:
        return None
    try:
        if raw != raw:  # NaN
            return None
    except (TypeError, ValueError):
        pass
    try:
        sec = float(raw) * 86400.0
    except (TypeError, ValueError):
        return None
    if not math.isfinite(sec):
        return None
    return sec


def score_oc_frame(frame: Any) -> dict:
    """Score a stage-1 O–C DataFrame (cycle_index + oc_detrended_days or shift_days)."""
    if frame is None or getattr(frame, "empty", True):
        return dict(EMPTY)
    cols = list(frame.columns)
    cyc_i = (
        cols.index("cycle_index")
        if "cycle_index" in cols
        else (cols.index("cycle") if "cycle" in cols else None)
    )
    det_i = cols.index("oc_detrended_days") if "oc_detrended_days" in cols else None
    shift_i = cols.index("shift_days") if "shift_days" in cols else None
    pairs = []
    for i, rec in enumerate(frame.itertuples(index=False, name=None)):
        cyc = rec[cyc_i] if cyc_i is not None else i
        det = rec[det_i] if det_i is not None else None
        shift = rec[shift_i] if shift_i is not None else None
        pairs.append((cyc, oc_seconds(det, shift)))
    return score_oc_pairs(pairs)


def acf_maxlag(n: int) -> int:
    """Lag budget for the O-C autocorrelation, matching stage 1's lb_maxlag."""
    return max(5, min(60, n // 3))


def acf_series(values: Iterable[float | None]) -> list[float]:
    """Autocorrelation of a finite O-C series at lags 1..acf_maxlag.

    Used by the viewer's ACF panel; kept here next to spark_path so the page
    never has to recompute anything from the O-C itself.
    """
    xs = [float(v) for v in values if v is not None and math.isfinite(v)]
    n = len(xs)
    if n < 8:
        return []
    maxlag = acf_maxlag(n)
    mean = sum(xs) / n
    c = [x - mean for x in xs]
    denom = sum(v * v for v in c)
    if denom <= 0:
        return []
    return [sum(c[i] * c[i + k] for i in range(n - k)) / denom for k in range(1, maxlag + 1)]
