#!/usr/bin/env python3
"""
stage1_pulse.py -- Stage 1: search TESS, Kepler, ASAS-SN and ZTF light curves of
Cepheids, RR Lyraes and beta Cephei stars for artificial phase modulation, following
Learned et al. (2008, arXiv:0809.0339) and Hippke et al. (2015, ApJ 798, 42).

Star classes (--star-class; 'auto' reads the '# star_class:' header that
stage0_build_targets.py writes, else every class present in the store):
  cepheid  P = 0.05-60 d, single dominant mode, no prewhitening.
  rrlyrae  P = 0.1-2 d; Blazhko shows up as a secondary period. On TESS data
           up to 3 independent modes are prewhitened first, so a double-mode
           star's beat is not read as a two-level code (KIC 5520878).
  bcep     P = 0.03-1 d, low-amplitude and usually multiperiodic: before the
           per-cycle timing fit the light curve is prewhitened of up to
           --prewhiten independent frequencies (not harmonics of the dominant
           mode) so beating between modes is not read as O-C structure.

Two analysis modes, chosen automatically from the cadence:
  DENSE  (TESS 2-min / 20-s, Kepler 1-min):  fit every pulsation cycle with a
         time-shifted template -> per-cycle O-C.
  SPARSE (ASAS-SN / ZTF, ~nightly): convert each epoch's residual into an implied
         timing offset via the template slope, average over blocks of cycles
         -> per-block O-C; plus a residual^2-vs-slope^2 jitter test that
         detects cycle-to-cycle phase scatter block averaging would erase.

Tests: autocorrelation vs shuffle null (incl. prime-lag), Ljung-Box,
permutation entropy, >5-sigma jumps, prewhitening for secondary periods
(multi-mode / Blazhko), sampling-alias check for ground-based data.

Outputs: summary, triage, per-star O-C series, and a compact phase-fold
into the SQLite store (default data/cepheid.sqlite). Light curves are cached
under data/lc_cache. Stage 2 reads current triage from the store. Open
candidates.html against the store (`starseti serve`) to inspect folds and
timing.

A bare --target-file name is looked up in data/starlists/, where stage 0
writes them; a path that exists is used as given.

Usage
-----
  starseti stage1 --synthetic [--inject] [--survey tess|asassn|ztf|all]
  starseti stage1                                    # all store rows (TESS + ASAS-SN + ZTF)
  starseti stage1 --survey tess   --target-file tess_cepheids.txt --jobs 8
  starseti stage1 --survey asassn --targets "RR Lyr" "SU Cas" --band g
  starseti stage1 --survey asassn --target-file asassn_cepheids.txt   # 'ASASSN <id>' lines
  starseti stage1 --survey ztf    --target-file ztf_cepheids.txt      # 'ZTF <oid>' lines
  starseti stage1 --survey ztf    --targets "RR Lyr" "RADEC 291.366 42.784" --band r
  starseti stage1 --survey all    --target-file stars.txt
  starseti stage1 --survey tess   --target-file tess_rrlyrae.txt   # header sets rrlyrae
  starseti stage1 --survey tess   --target-file tess_bcep.txt --cadence fast
  starseti stage1 --star-class bcep --targets "BW Vul" "beta Cep" "nu Eri"
  starseti stage1 --synthetic --inject --star-class cepheid   # injected primes control
  starseti stage1 --synthetic --inject --star-class bcep      # multi-mode + Fibonacci
  starseti stage1 --synthetic --inject --star-class rrlyrae   # RR Lyr + counting
  starseti stage1 --triage          # recompute triage flags and rescore
  starseti stage1 --retry-errors    # re-run download / analysis failures
  python stage1_pulse.py --skip-existing   # resume: skip stars that already have a row
  python stage1_pulse.py --skip-existing --retry-errors   # resume, retrying error rows
  python stage2_vet.py                     # current triage that is not yet vetted

ZTF notes: IRSA's light-curve service (nph_light_curves) is used per target;
targets are 'ZTF <oid>', 'TIC N', 'RADEC ra dec' or a SIMBAD name. All objects
within 1.5" (overlapping fields / quadrants) are merged with per-object
zero-point offsets. Public survey cadence is ~2-3 nights per band, so ZTF
always runs in SPARSE mode. For very large lists the service throttles; keep
--jobs modest (4) or pre-download via the DR HATS catalogs with lsdb.

Requirements: numpy scipy astropy pandas lightkurve astroquery tqdm skypatrol requests
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import os
import re
import socket
import sys
import time
import warnings
from collections.abc import Iterator, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from html import unescape
from itertools import pairwise, permutations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm

try:
    from numba import njit

    HAVE_NUMBA = True
except ImportError:  # pragma: no cover
    HAVE_NUMBA = False

    def njit(*args: object, **kwargs: object) -> Any:
        if args and callable(args[0]) and not kwargs:
            return args[0]

        def _wrap(fn: Any) -> Any:
            return fn

        return _wrap


import layout
import lc_cache
import matched
import score
import sequence
from pools import worker_pool
from store import Store, add_db_arg, json_safe

warnings.filterwarnings("ignore")

# No network call in this program may block forever. MAST, Sky Patrol and
# IRSA are reached through libraries that mostly leave the socket timeout at
# its default of None, and a stalled connection then holds a pool worker
# indefinitely. A process-wide default is the one place that covers all of
# them; requests and urllib3 fall back to it when no timeout is passed.
NETWORK_TIMEOUT_S = 300.0
socket.setdefaulttimeout(NETWORK_TIMEOUT_S)
# Import first so each package can install its AstropyLogger; getLogger(name)
# beforehand leaves a plain Logger and _init_log() crashes on _set_defaults.
import astropy  # noqa: E402, F401
from astropy import log as _astropy_log  # noqa: E402

_astropy_log.setLevel("WARNING")
import astroquery  # noqa: E402

astroquery.log.setLevel("WARNING")
logging.getLogger("lightkurve").setLevel(logging.ERROR)
logging.getLogger("lightkurve.utils").setLevel(logging.ERROR)


@njit(cache=True)
def _fourier_at(
    t: np.ndarray, coef: np.ndarray, period: float, t0: float, nharm: int, dt: float
) -> np.ndarray:
    n = t.shape[0]
    out = np.empty(n)
    c0 = coef[0]
    for i in range(n):
        phase = ((t[i] - dt - t0) / period) % 1.0
        s = c0
        for k in range(1, nharm + 1):
            ang = 2.0 * np.pi * k * phase
            s += coef[2 * k - 1] * np.cos(ang) + coef[2 * k] * np.sin(ang)
        out[i] = s
    return out


@njit(cache=True)
def _fourier_slope_nb(
    t: np.ndarray, coef: np.ndarray, period: float, t0: float, nharm: int
) -> np.ndarray:
    n = t.shape[0]
    out = np.empty(n)
    for i in range(n):
        phase = ((t[i] - t0) / period) % 1.0
        d = 0.0
        for k in range(1, nharm + 1):
            ang = 2.0 * np.pi * k * phase
            a, b = coef[2 * k - 1], coef[2 * k]
            d += 2.0 * np.pi * k / period * (-a * np.sin(ang) + b * np.cos(ang))
        out[i] = d
    return out


@njit(cache=True)
def _chi2_one(
    t: np.ndarray,
    f: np.ndarray,
    invvar: np.ndarray,
    coef: np.ndarray,
    period: float,
    t0: float,
    nharm: int,
    dt: float,
) -> float:
    s = 0.0
    c0 = coef[0]
    for i in range(t.shape[0]):
        phase = ((t[i] - dt - t0) / period) % 1.0
        m = c0
        for k in range(1, nharm + 1):
            ang = 2.0 * np.pi * k * phase
            m += coef[2 * k - 1] * np.cos(ang) + coef[2 * k] * np.sin(ang)
        r = f[i] - m
        s += r * r * invvar[i]
    return s


@njit(cache=True)
def _chi2_grid(
    t: np.ndarray,
    f: np.ndarray,
    invvar: np.ndarray,
    coef: np.ndarray,
    period: float,
    t0: float,
    nharm: int,
    dts: np.ndarray,
) -> np.ndarray:
    """Chi2 vs trial time-shifts: one column per dt (vectorized over the grid)."""
    n = t.shape[0]
    g = dts.shape[0]
    out = np.empty(g)
    c0 = coef[0]
    for j in range(g):
        dt = dts[j]
        s = 0.0
        for i in range(n):
            phase = ((t[i] - dt - t0) / period) % 1.0
            m = c0
            for k in range(1, nharm + 1):
                ang = 2.0 * np.pi * k * phase
                m += coef[2 * k - 1] * np.cos(ang) + coef[2 * k] * np.sin(ang)
            r = f[i] - m
            s += r * r * invvar[i]
        out[j] = s
    return out


@njit(cache=True)
def _dense_timings_nb(
    t: np.ndarray,
    f: np.ndarray,
    invvar: np.ndarray,
    cyc: np.ndarray,
    period: float,
    coef: np.ndarray,
    t0: float,
    nharm: int,
    min_pts: int,
    min_coverage: float,
) -> tuple:
    order = np.argsort(cyc)
    n = t.shape[0]
    ngrid = 121
    dts = np.linspace(-0.15 * period, 0.15 * period, ngrid)
    nfine = 81
    out_c = np.empty(n, np.int64)
    out_tc = np.empty(n)
    out_dt = np.empty(n)
    out_err = np.empty(n)
    out_npts = np.empty(n, np.int64)
    out_chi2 = np.empty(n)
    out_amp = np.empty(n)
    nout = 0
    i = 0
    while i < n:
        c = cyc[order[i]]
        j = i + 1
        while j < n and cyc[order[j]] == c:
            j += 1
        npts = j - i
        if npts >= min_pts:
            tt = np.empty(npts)
            ff = np.empty(npts)
            iv = np.empty(npts)
            tmin = 1e300
            tmax = -1e300
            tsum = 0.0
            for k in range(npts):
                idx = order[i + k]
                tt[k] = t[idx]
                ff[k] = f[idx]
                iv[k] = invvar[idx]
                if tt[k] < tmin:
                    tmin = tt[k]
                if tt[k] > tmax:
                    tmax = tt[k]
                tsum += tt[k]
            if (tmax - tmin) / period >= min_coverage:
                chi2g = _chi2_grid(tt, ff, iv, coef, period, t0, nharm, dts)
                best = 0
                bestv = chi2g[0]
                for g in range(1, ngrid):
                    if chi2g[g] < bestv:
                        bestv = chi2g[g]
                        best = g
                g0 = dts[best]
                half = 0.012 * period
                fine = np.linspace(g0 - half, g0 + half, nfine)
                chi2f = _chi2_grid(tt, ff, iv, coef, period, t0, nharm, fine)
                fb = 0
                fbv = chi2f[0]
                for g in range(1, nfine):
                    if chi2f[g] < fbv:
                        fbv = chi2f[g]
                        fb = g
                dt = fine[fb]
                # Clip outliers against the *shifted* model, per cycle, and
                # refit. The template's own 5-sigma clip is taken against the
                # mean ephemeris, so a cycle that genuinely arrives early has
                # its steep-branch points -- the ones that carry the timing --
                # thrown out, and the shift is read 30-40% low on real hosts
                # (research/injection-sensitivity.md). Here the reference is
                # the cycle's own fine-grid shift, so only points that
                # disagree with the cycle are dropped; never more than half.
                mod0 = _fourier_at(tt, coef, period, t0, nharm, dt)
                res = np.empty(npts)
                for k in range(npts):
                    res[k] = ff[k] - mod0[k]
                scale = 1.4826 * np.median(np.abs(res - np.median(res)))
                if scale > 0.0:
                    nkeep = 0
                    for k in range(npts):
                        if abs(res[k]) < 5.0 * scale:
                            nkeep += 1
                    if nkeep >= min_pts and 2 * nkeep >= npts and nkeep < npts:
                        tt2 = np.empty(nkeep)
                        ff2 = np.empty(nkeep)
                        iv2 = np.empty(nkeep)
                        m = 0
                        tsum = 0.0
                        for k in range(npts):
                            if abs(res[k]) < 5.0 * scale:
                                tt2[m] = tt[k]
                                ff2[m] = ff[k]
                                iv2[m] = iv[k]
                                tsum += tt[k]
                                m += 1
                        tt, ff, iv = tt2, ff2, iv2
                        npts = nkeep
                        chi2f = _chi2_grid(tt, ff, iv, coef, period, t0, nharm, fine)
                        fb = 0
                        fbv = chi2f[0]
                        for g in range(1, nfine):
                            if chi2f[g] < fbv:
                                fbv = chi2f[g]
                                fb = g
                        dt = fine[fb]
                h = 1e-4 * period
                c0 = _chi2_one(tt, ff, iv, coef, period, t0, nharm, dt)
                cp = _chi2_one(tt, ff, iv, coef, period, t0, nharm, dt + h)
                cm = _chi2_one(tt, ff, iv, coef, period, t0, nharm, dt - h)
                curv = (cp - 2.0 * c0 + cm) / (h * h)
                chi2_red = c0 / max(npts - 1, 1)
                if curv > 0.0:
                    err = np.sqrt(2.0 / curv) * np.sqrt(max(chi2_red, 1.0))
                else:
                    err = np.nan
                mod = _fourier_at(tt, coef, period, t0, nharm, dt)
                num = 0.0
                den = 0.0
                cmean = coef[0]
                for k in range(npts):
                    mk = mod[k] - cmean
                    num += mk * (ff[k] - cmean)
                    den += mk * mk
                amp = num / den if den > 1e-12 else 1.0
                out_c[nout] = c
                out_tc[nout] = tsum / npts
                out_dt[nout] = dt
                out_err[nout] = err
                out_npts[nout] = npts
                out_chi2[nout] = chi2_red
                out_amp[nout] = amp
                nout += 1
        i = j
    return (
        nout,
        out_c[:nout],
        out_tc[:nout],
        out_dt[:nout],
        out_err[:nout],
        out_npts[:nout],
        out_chi2[:nout],
        out_amp[:nout],
    )


# How often `_ls_power_grid_nb` re-seeds its rotor from an exact grid
# frequency. Drift is not what limits the kernel -- never re-seeding costs
# under 2e-14 in power, measured against a 60-digit reference -- so this is
# margin, and it buys it for about a tenth of the kernel's time.
LS_RESEED = 32


@njit(cache=True)
def _ls_power_grid_nb(
    t: np.ndarray, y: np.ndarray, w: np.ndarray, freqs: np.ndarray, reseed: int
) -> np.ndarray:
    """Generalized Lomb-Scargle power on a uniform frequency grid.

    The quantity astropy returns for `LombScargle(t, y, dy).power(f,
    method="cython")` with its defaults -- the Zechmeister-Kurster form with a
    floating mean, centred data and standard normalization. Astropy makes two
    passes per frequency, each with its own sin and cos; every sum either pass
    needs can be had from one pass instead, because the time shift tau enters
    only as a single angle and the double-angle terms come from the first
    pass's sine and cosine. Along a uniform grid it needs almost no trig at
    all: exp(i w_(i+1) t) is exp(i w_i t) turned by a per-point rotor that
    does not change.

    Not turning w*t into a sine at every frequency is also what makes it more
    accurate than the direct form on ground-based data. Over an eight-year
    baseline at half a day that product is ~1e5 radians, where one ulp is
    1.5e-11 rad, and the rounding of the product is then the error floor. The
    rotor accumulates phase from one seed and never forms the big product:
    against a 60-digit reference this matches astropy on TESS sampling and
    beats it on ASAS-SN and ZTF.

    Two things it is sensitive to, and both are the caller's to get right by
    handing over the grid itself. The step has to be taken across the grid's
    ends: `freqs[1] - freqs[0]` loses eleven digits to cancellation, and the
    error then multiplies by the step index, which is a 150x accuracy loss by
    the end of a 801-point window. And a re-seed has to land on a grid
    frequency exactly, because one ulp of frequency is 6e-14 rad of phase out
    at these baselines.
    """
    n = t.shape[0]
    nf = freqs.shape[0]
    out = np.empty(nf)
    W = 0.0
    Ysum = 0.0
    YYsum = 0.0
    for j in range(n):
        W += w[j]
        Ysum += w[j] * y[j]
        YYsum += w[j] * y[j] * y[j]
    if not (W > 0.0):
        for i in range(nf):
            out[i] = np.nan
        return out
    Y = Ysum / W
    YY = YYsum / W
    zr = np.empty(n)
    zi = np.empty(n)
    rr = np.empty(n)
    ri = np.empty(n)
    df = (freqs[nf - 1] - freqs[0]) / (nf - 1) if nf > 1 else 0.0
    dom = 2.0 * np.pi * df
    for j in range(n):
        rr[j] = np.cos(dom * t[j])
        ri[j] = np.sin(dom * t[j])
    for i in range(nf):
        if i % reseed == 0:
            om = 2.0 * np.pi * freqs[i]
            for j in range(n):
                zr[j] = np.cos(om * t[j])
                zi[j] = np.sin(om * t[j])
        else:
            for j in range(n):
                a = zr[j] * rr[j] - zi[j] * ri[j]
                b = zr[j] * ri[j] + zi[j] * rr[j]
                zr[j] = a
                zi[j] = b
        # Weighted sums of cos, sin, y cos, y sin, cos 2wt and sin 2wt: every
        # term of the closed form is built from these six.
        A = 0.0
        B = 0.0
        Ay = 0.0
        By = 0.0
        A2 = 0.0
        B2 = 0.0
        for j in range(n):
            c = zr[j]
            s = zi[j]
            wj = w[j]
            A += wj * c
            B += wj * s
            Ay += wj * y[j] * c
            By += wj * y[j] * s
            A2 += wj * (1.0 - 2.0 * s * s)
            B2 += wj * 2.0 * s * c
        Cm = A / W
        Sm = B / W
        C2m = A2 / W - (Cm * Cm - Sm * Sm)
        S2m = B2 / W - 2.0 * Sm * Cm
        # w*tau, the shift that makes the sine and cosine terms orthogonal.
        phi = 0.5 * np.arctan2(S2m, C2m)
        cp = np.cos(phi)
        sp = np.sin(phi)
        c2p = np.cos(2.0 * phi)
        s2p = np.sin(2.0 * phi)
        Ctau = (cp * A + sp * B) / W
        Stau = (cp * B - sp * A) / W
        YCtau = (cp * Ay + sp * By) / W - Y * Ctau
        YStau = (cp * By - sp * Ay) / W - Y * Stau
        half = 0.5 * (c2p * A2 + s2p * B2) / W
        CCtau = 0.5 + half - Ctau * Ctau
        SStau = 0.5 - half - Stau * Stau
        vy = YY - Y * Y
        # A constant series has no periodogram at all: astropy answers NaN and
        # `_refine_freq` reads an all-NaN window as "keep the coarse peak".
        # numba's default error model raises on a zero divide rather than
        # propagating it, so the degenerate cases are named rather than hit.
        if not (vy > 0.0) or not (CCtau > 0.0) or not (SStau > 0.0):
            out[i] = np.nan
        else:
            out[i] = (YCtau * YCtau / CCtau + YStau * YStau / SStau) / vy
    return out


@njit(cache=True)
def _pairwise_sum_nb(a: np.ndarray) -> float:
    """numpy's own summation order, so a numba reduction matches np.sum exactly.

    `np.add.reduce` over a contiguous float64 array sums in pairs -- eight
    running accumulators within a 128-element block, halves above that -- and
    a plain accumulation loop disagrees with it in the last bit or two. The
    null loop below has to reproduce statistics that were computed with numpy
    before it existed, value for value, so it sums the way numpy sums.
    """
    n = a.shape[0]
    if n < 8:
        res = 0.0
        for i in range(n):
            res += a[i]
        return res
    if n <= 128:
        r0 = a[0]
        r1 = a[1]
        r2 = a[2]
        r3 = a[3]
        r4 = a[4]
        r5 = a[5]
        r6 = a[6]
        r7 = a[7]
        i = 8
        stop = n - (n % 8)
        while i < stop:
            r0 += a[i]
            r1 += a[i + 1]
            r2 += a[i + 2]
            r3 += a[i + 3]
            r4 += a[i + 4]
            r5 += a[i + 5]
            r6 += a[i + 6]
            r7 += a[i + 7]
            i += 8
        res = ((r0 + r1) + (r2 + r3)) + ((r4 + r5) + (r6 + r7))
        while i < n:
            res += a[i]
            i += 1
        return res
    n2 = n // 2
    n2 -= n2 % 8
    return _pairwise_sum_nb(a[:n2]) + _pairwise_sum_nb(a[n2:])


@njit(cache=True)
def _acf_nb(x: np.ndarray, maxlag: int) -> np.ndarray:
    x = x - np.mean(x)
    d = np.dot(x, x)
    out = np.empty(maxlag + 1)
    out[0] = 1.0
    if d <= 0.0:
        for k in range(1, maxlag + 1):
            out[k] = 0.0
        return out
    for k in range(1, maxlag + 1):
        out[k] = np.dot(x[:-k], x[k:]) / d
    return out


@njit(cache=True)
def _decode_mdl_nb(bits: np.ndarray, cyc: np.ndarray, maxlen: int) -> tuple:
    """Shortest description of an O-C series read as a clock-locked bit word.

    The same quantity `sequence._score_bit_period` minimises, in the same
    units: for each candidate word length L the series is folded on
    `cycle % L`, each residue is majority-voted into a bit, and the cost is
    L bits for the word plus the binary entropy of the votes that disagreed.
    A clean two-level code locked to the cycle count compresses to about L; a
    series with no clock in it does not compress.

    It lives here rather than in `sequence.py` because it is called once per
    shuffle draw -- `sequence.py` is deliberately stdlib-only so the server can
    import it, and a pure-Python fold inside the null loop costs seconds per
    star. Returns (mdl, L); L is 0 when no length could be folded.
    """
    n = bits.shape[0]
    best = 1e300
    best_l = 0
    acc = np.empty(maxlen + 1, np.int64)
    cnt = np.empty(maxlen + 1, np.int64)
    for length in range(4, maxlen + 1):
        for k in range(length):
            acc[k] = 0
            cnt[k] = 0
        for i in range(n):
            k = int(cyc[i]) % length
            if k < 0:
                k += length
            acc[k] += bits[i]
            cnt[k] += 1
        residual = 0.0
        ok = True
        for k in range(length):
            if cnt[k] < 1:
                ok = False
                break
            ones = acc[k]
            zeros = cnt[k] - ones
            fewer = ones if ones < zeros else zeros
            if fewer > 0:
                e = fewer / cnt[k]
                residual -= cnt[k] * (e * np.log2(e) + (1.0 - e) * np.log2(1.0 - e))
        if ok:
            mdl = length + residual
            if mdl < best:
                best = mdl
                best_l = length
    return best, best_l


@njit(cache=True)
def _pct_sorted(srt: np.ndarray, q: float) -> float:
    """Linear-interpolated percentile of an already sorted array."""
    n = srt.shape[0]
    pos = (n - 1) * q
    lo = int(np.floor(pos))
    hi = int(np.ceil(pos))
    if lo == hi:
        return srt[lo]
    return srt[lo] + (srt[hi] - srt[lo]) * (pos - lo)


@njit(cache=True)
def _fold_index_table(cyc: np.ndarray, n: int, maxlen: int, ks: np.ndarray) -> tuple:
    """Every (bit period, word length) fold the cost search will try.

    Which residue a step falls in is a fact about the cycle numbers, not about
    the values being folded, so the table is the same for the observation and
    for every one of its shuffle draws -- and building it once takes the two
    integer divisions per point per candidate out of the null loop entirely,
    which is where that loop spent most of its time. Candidates with an empty
    residue are dropped here rather than per draw; that is the same `cnt < 1`
    test, and a dropped candidate contributed nothing to the search.

    Returns (tab, cnt, lens, kk): the residue of each point per candidate, how
    many points each residue holds, and the (L, k) the candidate stands for,
    in the order the search walks them -- ties in cost keep the first, so the
    order is part of the answer.
    """
    rows = 0
    for ki in range(ks.shape[0]):
        k = int(ks[ki])
        if k < 1:
            continue
        lmax = min(maxlen, n // (2 * k))
        if lmax >= 4:
            rows += lmax - 3
    tab = np.empty((rows, n), np.int16)
    cnt = np.empty((rows, maxlen + 1), np.int64)
    lens = np.empty(rows, np.int64)
    kk = np.empty(rows, np.int64)
    row = 0
    for ki in range(ks.shape[0]):
        k = int(ks[ki])
        if k < 1:
            continue
        lmax = min(maxlen, n // (2 * k))
        for length in range(4, lmax + 1):
            for j in range(length):
                cnt[row, j] = 0
            for i in range(n):
                j = (int(cyc[i]) // k) % length
                if j < 0:
                    j += length
                tab[row, i] = j
                cnt[row, j] += 1
            ok = True
            for j in range(length):
                if cnt[row, j] < 1:
                    ok = False
                    break
            if ok:
                lens[row] = length
                kk[row] = k
                row += 1
    return tab[:row], cnt[:row], lens[:row], kk[:row]


@njit(cache=True)
def _fold_cost_tab_nb(
    y: np.ndarray, tab: np.ndarray, cnt: np.ndarray, lens: np.ndarray, kk: np.ndarray
) -> tuple:
    """`_fold_cost_nb` against a pre-built fold table. See it for the quantity."""
    n = y.shape[0]
    best = 1e300
    best_l = 0
    best_k = 0
    maxlen = cnt.shape[1] - 1
    acc = np.empty(maxlen + 1)
    means = np.empty(maxlen + 1)
    word = np.empty(maxlen + 1, np.int64)
    # Sorted into a buffer rather than through np.sort: this is the innermost
    # loop of the whole null, and one allocation per candidate per draw was
    # costing more than the sort. A sort is a permutation, so the order it
    # walks is not part of the answer.
    srt = np.empty(maxlen + 1)
    for row in range(tab.shape[0]):
        length = int(lens[row])
        for j in range(length):
            acc[j] = 0.0
        for i in range(n):
            acc[tab[row, i]] += y[i]
        for j in range(length):
            means[j] = acc[j] / cnt[row, j]
            srt[j] = means[j]
        # Otsu split of the L means (sequence.split_means): the cut that
        # explains most of their variance, one mean or more on each side.
        srt[:length].sort()
        if not srt[length - 1] > srt[0]:
            continue
        total = 0.0
        for j in range(length):
            total += srt[j]
        bestb = -1.0
        mid = 0.0
        run = 0.0
        for i in range(1, length):
            run += srt[i - 1]
            m0 = run / i
            m1 = (total - run) / (length - i)
            between = (i / length) * (1.0 - i / length) * (m0 - m1) ** 2
            if between > bestb:
                bestb = between
                mid = (srt[i - 1] + srt[i]) / 2.0
        for j in range(length):
            word[j] = 1 if means[j] < mid else 0
        s0 = 0.0
        s1 = 0.0
        c0 = 0
        c1 = 0
        for i in range(n):
            if word[tab[row, i]] == 1:
                s1 += y[i]
                c1 += 1
            else:
                s0 += y[i]
                c0 += 1
        if c0 == 0 or c1 == 0:
            continue
        l0 = s0 / c0
        l1 = s1 / c1
        rss = 0.0
        for i in range(n):
            d = y[i] - (l1 if word[tab[row, i]] == 1 else l0)
            rss += d * d
        if rss <= 0.0:
            rss = 1e-300
        # ln 2 per bit, ln n per level: see sequence._fold_word.
        cost = n * np.log(rss / n) + length * np.log(2.0) + 2.0 * np.log(n)
        if cost < best:
            best = cost
            best_l = length
            best_k = int(kk[row])
    return best, best_l, best_k


@njit(cache=True)
def _fold_cost_nb(y: np.ndarray, cyc: np.ndarray, maxlen: int, ks: np.ndarray) -> tuple:
    """Fold-first description length of the O-C steps read as a bit word.

    For every bit period k in `ks` and word length L the steps are averaged
    per residue of (cycle // k) % L, the means are split at the 20/80 mid
    (sequence.two_level_mid), and the cost is n ln(RSS/n) for the residuals
    about the two levels, ln 2 per bit and ln n per level --
    the same quantity `sequence._fold_word` minimises. Averaging before
    thresholding is what lets a word repeated r times be read at sqrt(r) the
    per-cycle signal-to-noise. Returns (cost, L, k); L is 0 when nothing
    could be folded.
    """
    tab, cnt, lens, kk = _fold_index_table(cyc, y.shape[0], maxlen, ks)
    return _fold_cost_tab_nb(y, tab, cnt, lens, kk)


@njit(cache=True)
def _otsu_nb(y: np.ndarray) -> float:
    """How bimodal the steps are: Otsu's between-class share of the variance.

    A phase code's steps are quantised -- exactly zero or exactly dt -- so
    their histogram has two peaks with a gap, while the steps of a natural
    random-walk O-C take every size. The statistic is the largest fraction of
    the total variance a single split explains, over splits that leave at
    least a tenth of the points on each side, so it needs no repeat of the
    word at all. Scale-free, so it is calibrated by the same shuffle null as
    everything else.
    """
    n = y.shape[0]
    if n < 12:
        return np.nan
    srt = np.sort(y)
    tot = 0.0
    for i in range(n):
        tot += srt[i]
    mean = tot / n
    var = 0.0
    for i in range(n):
        var += (srt[i] - mean) ** 2
    var /= n
    if var <= 0.0:
        return 0.0
    nmin = max(3, n // 10)
    best = 0.0
    run = 0.0
    for i in range(1, n):
        run += srt[i - 1]
        if i < nmin or n - i < nmin:
            continue
        w0 = i / n
        m0 = run / i
        m1 = (tot - run) / (n - i)
        between = w0 * (1.0 - w0) * (m0 - m1) ** 2
        if between > best:
            best = between
    return best / var


@njit(cache=True)
def _pct_linear_nb(srt: np.ndarray, q: float) -> float:
    """`np.percentile(y, 100*q)` on an already sorted array, to the last bit.

    numpy interpolates from whichever end is nearer -- `a + (b-a)*t` below the
    midpoint and `b - (b-a)*(1-t)` above it -- and the one-sided form disagrees
    with it about twice in a hundred. `_pct_sorted` is the one-sided form and
    stays that way; the clip below is a threshold test, so it wants the same
    limit numpy would have computed rather than one that rounds differently.
    """
    n = srt.shape[0]
    pos = (n - 1) * q
    lo = int(np.floor(pos))
    hi = int(np.ceil(pos))
    if lo == hi:
        return srt[lo]
    t = pos - lo
    d = srt[hi] - srt[lo]
    if t >= 0.5:
        return srt[hi] - d * (1.0 - t)
    return srt[lo] + d * t


@njit(cache=True)
def _clip_steps_nb(
    y: np.ndarray,
    out: np.ndarray,
    sigma: float,
    max_frac: float,
    kick_min_steps: int,
    kick_spread: float,
) -> None:
    """`sequence.clip_steps` for one shuffle draw, written into `out`.

    The observation is winsorized before it is folded, so a draw has to be
    winsorized before it is folded too -- otherwise the null is a distribution
    of a different statistic than the one being tested, which is how the
    decode p-value lost its calibration the first time. Kept here in compiled
    form rather than called out to `sequence`, which is deliberately
    stdlib-only, because it runs once per draw. `kick_min_steps` and
    `kick_spread` are `sequence.one_sided_cluster`'s rule: a cluster of
    same-sized outliers on one side of the median is a sparse code's kicks
    and is left alone.
    """
    n = y.shape[0]
    for i in range(n):
        out[i] = y[i]
    if n < 8:
        return
    srt = np.sort(y)
    q1 = _pct_linear_nb(srt, 0.25)
    med = _pct_linear_nb(srt, 0.5)
    q3 = _pct_linear_nb(srt, 0.75)
    scale = (q3 - q1) / 1.349
    if not scale > 0.0:
        return
    lim = sigma * scale
    beyond = 0
    for i in range(n):
        if abs(y[i] - med) > lim:
            beyond += 1
    # An outlier is rare by definition: a word's minority level, a sparse
    # word's set bits and a square wave's transitions are structure, and
    # pulling them to one value would build the level the fold looks for.
    if beyond == 0 or beyond / n > max_frac:
        return
    if beyond >= kick_min_steps:
        # One-sided and of one size: kicks, not glitches (sequence.one_sided_cluster).
        dist = np.empty(beyond)
        j = 0
        same = True
        first_sign = 0.0
        for i in range(n):
            d = y[i] - med
            if abs(d) > lim:
                sgn = 1.0 if d > 0 else -1.0
                if j == 0:
                    first_sign = sgn
                elif sgn != first_sign:
                    same = False
                dist[j] = abs(d)
                j += 1
        if same:
            dist.sort()
            dmed = _pct_linear_nb(dist, 0.5)
            if dmed > 0.0 and (dist[beyond - 1] - dist[0]) <= kick_spread * dmed:
                return
    at_lo = med - lim
    at_hi = med + lim
    for i in range(n):
        if out[i] < at_lo:
            out[i] = at_lo
        elif out[i] > at_hi:
            out[i] = at_hi


@njit(cache=True)
def _level_fold_cost_tab_nb(
    y: np.ndarray, tab: np.ndarray, cnt: np.ndarray, lens: np.ndarray, kk: np.ndarray
) -> tuple:
    """`sequence._fold_level_word`'s cost search against a pre-built fold table.

    The `fold` model: the *levels* are averaged per residue, the bits are the
    cyclic first differences of the means split at Otsu's cut, and the cost
    is the description length of the levels about the staircase the word
    predicts (residue j at a - s (S_j - W j / L), fitted by weighted least
    squares), n ln(RSS/n) + L ln 2 + 2 ln n. Same currency as
    `_fold_cost_tab_nb`, for a series with no consecutive cycles to
    difference. Returns (cost, L, k); L is 0 when nothing could be folded.
    """
    n = y.shape[0]
    best = 1e300
    best_l = 0
    best_k = 0
    maxlen = cnt.shape[1] - 1
    acc = np.empty(maxlen + 1)
    acc2 = np.empty(maxlen + 1)
    means = np.empty(maxlen + 1)
    diffs = np.empty(maxlen + 1)
    srt = np.empty(maxlen + 1)
    word = np.empty(maxlen + 1, np.int64)
    u = np.empty(maxlen + 1)
    for row in range(tab.shape[0]):
        length = int(lens[row])
        for j in range(length):
            acc[j] = 0.0
            acc2[j] = 0.0
        for i in range(n):
            j = tab[row, i]
            acc[j] += y[i]
            acc2[j] += y[i] * y[i]
        for j in range(length):
            means[j] = acc[j] / cnt[row, j]
        for j in range(length):
            diffs[j] = means[j] - means[j - 1 if j > 0 else length - 1]
            srt[j] = diffs[j]
        srt[:length].sort()
        if not srt[length - 1] > srt[0]:
            continue
        total = 0.0
        for j in range(length):
            total += srt[j]
        bestb = -1.0
        cut = 0.0
        run = 0.0
        for i in range(1, length):
            run += srt[i - 1]
            m0 = run / i
            m1 = (total - run) / (length - i)
            between = (i / length) * (1.0 - i / length) * (m0 - m1) ** 2
            if between > bestb:
                bestb = between
                cut = (srt[i - 1] + srt[i]) / 2.0
        w_total = 0
        for j in range(length):
            word[j] = 1 if diffs[j] < cut else 0
            w_total += word[j]
        if w_total == 0 or w_total == length:
            continue
        # The staircase shape the word predicts, with its ramp removed.
        running = 0.0
        for j in range(length):
            running += word[j]
            u[j] = running - w_total * j / length
        sw = 0.0
        su = 0.0
        sm = 0.0
        suu = 0.0
        sum_ = 0.0
        for j in range(length):
            c = cnt[row, j]
            sw += c
            su += c * u[j]
            sm += c * means[j]
            suu += c * u[j] * u[j]
            sum_ += c * u[j] * means[j]
        det = sw * suu - su * su
        if not det > 0.0:
            continue
        s = -(sw * sum_ - su * sm) / det
        a = (sm + s * su) / sw
        rss = 0.0
        for j in range(length):
            c = cnt[row, j]
            level = a - s * u[j]
            rss += acc2[j] - c * means[j] * means[j] + c * (means[j] - level) ** 2
        if rss <= 0.0:
            rss = 1e-300
        cost = n * np.log(rss / n) + length * np.log(2.0) + 2.0 * np.log(n)
        if cost < best:
            best = cost
            best_l = length
            best_k = int(kk[row])
    return best, best_l, best_k


@njit(cache=True)
def _epoch_null_loop_nb(
    y: np.ndarray,
    perms: np.ndarray,
    tab: np.ndarray,
    cnt: np.ndarray,
    lens: np.ndarray,
    kk: np.ndarray,
) -> np.ndarray:
    """The `fold` model's cost for a block of shuffle draws.

    A draw keeps the observed cycles and permutes the levels among them, so
    only the clock is under test; the (L, k) search runs on every draw, which
    is what puts the look-elsewhere over word length and bit period inside
    the p-value.
    """
    ndraw = perms.shape[0]
    n = y.shape[0]
    out = np.empty(ndraw)
    s = np.empty(n)
    for d in range(ndraw):
        for i in range(n):
            s[i] = y[perms[d, i]]
        out[d] = _level_fold_cost_tab_nb(s, tab, cnt, lens, kk)[0]
    return out


@njit(cache=True)
def _kick_set_stats_nb(
    cyc_obs: np.ndarray,
    sel: np.ndarray,
    seq_flat: np.ndarray,
    seq_off: np.ndarray,
    is_prime: np.ndarray,
) -> tuple:
    """The two kick statistics for one set of kicks.

    `cyc_obs` are the cycles that carry a measured step, sorted; `sel` the
    positions in it of the kicks, sorted. A gap between consecutive kicks
    counts only when every cycle between them was observed (the position
    difference equals the cycle difference); otherwise a kick could hide in
    it and it ends any run. Returns (longest run of consecutive gaps matching
    a window of one of the `seq_flat` sequences, index of that sequence,
    where in it the run begins, fraction of the kick pairs whose separation is
    prime).
    """
    m = sel.shape[0]
    gaps = np.empty(max(m - 1, 0), np.int64)
    for i in range(m - 1):
        dc = cyc_obs[sel[i + 1]] - cyc_obs[sel[i]]
        dp = sel[i + 1] - sel[i]
        gaps[i] = dc if dc == dp else -1
    best_run = 0
    best_seq = -1
    best_start = -1
    nseq = seq_off.shape[0] - 1
    for q in range(nseq):
        lo = seq_off[q]
        hi = seq_off[q + 1]
        for i in range(m - 1):
            if gaps[i] <= 0:
                continue
            for j in range(lo, hi):
                if seq_flat[j] != gaps[i]:
                    continue
                r = 0
                while i + r < m - 1 and j + r < hi and gaps[i + r] == seq_flat[j + r]:
                    r += 1
                if r > best_run:
                    best_run = r
                    best_seq = q
                    best_start = j - lo
    pairs = 0
    prime_pairs = 0
    nprime = is_prime.shape[0]
    for i in range(m):
        for j in range(i + 1, m):
            sep = cyc_obs[sel[j]] - cyc_obs[sel[i]]
            pairs += 1
            if sep < nprime and is_prime[sep]:
                prime_pairs += 1
    frac = prime_pairs / pairs if pairs > 0 else np.nan
    return best_run, best_seq, best_start, frac


@njit(cache=True)
def _kick_null_nb(
    cyc_obs: np.ndarray,
    sels: np.ndarray,
    seq_flat: np.ndarray,
    seq_off: np.ndarray,
    is_prime: np.ndarray,
) -> tuple:
    """`_kick_set_stats_nb` over a block of random kick sets (rows of `sels`,
    each sorted). Returns (runs, prime fractions) per draw."""
    ndraw = sels.shape[0]
    runs = np.empty(ndraw)
    fracs = np.empty(ndraw)
    for d in range(ndraw):
        r, _q, _s, f = _kick_set_stats_nb(cyc_obs, sels[d], seq_flat, seq_off, is_prime)
        runs[d] = r
        fracs[d] = f
    return runs, fracs


@njit(cache=True)
def _null_loop_nb(
    oc: np.ndarray,
    perms: np.ndarray,
    lo_idx: np.ndarray,
    hi_idx: np.ndarray,
    inc_cyc: np.ndarray,
    pr: np.ndarray,
    comp: np.ndarray,
    maxlag: int,
    maxlen: int,
    ks: np.ndarray,
    decodable: bool,
    quantisable: bool,
    clip_sigma: float,
    clip_max_frac: float,
    kick_min_steps: int,
    kick_spread: float,
) -> tuple:
    """Every shuffle statistic, for a block of draws, in one compiled loop.

    The draws arrive as rows of indices rather than being drawn here.
    numpy's `Generator.permutation` is what calibrates every p-value in the
    program, and `oc[rng.permutation(n)]` is exactly the array
    `rng.permutation(oc)` returns -- the shuffle swaps positions, so the index
    form and the value form consume the same stream and land on the same
    permutation. Together with `_pairwise_sum_nb` that makes this null the
    same one the Python loop produced, draw for draw and bit for bit; what
    goes away is the per-draw interpreter overhead, which at the ~30 O-C
    points a typical star has was most of the cost.

    Returns (max|ACF|, Ljung-Box Q, prime-lag contrast, fold cost, Otsu) per
    draw; a statistic this star cannot compute comes back as NaN, which
    `shuffle_pvalue` already drops.
    """
    ndraw = perms.shape[0]
    n = oc.shape[0]
    nstep = lo_idx.shape[0]
    out_max = np.empty(ndraw)
    out_lb = np.empty(ndraw)
    out_prime = np.empty(ndraw)
    out_mdl = np.empty(ndraw)
    out_quant = np.empty(ndraw)
    s = np.empty(n)
    lbterm = np.empty(maxlag)
    inc = np.empty(nstep)
    inc_c = np.empty(nstep)
    do_steps = decodable or quantisable
    do_prime = comp.shape[0] > 0
    # Once for the whole null: the fold table is a function of the cycle
    # numbers, which the shuffle does not touch.
    tab, tcnt, tlens, tkk = _fold_index_table(inc_cyc, nstep, maxlen, ks)
    for d in range(ndraw):
        for i in range(n):
            s[i] = oc[perms[d, i]]
        rs = _acf_nb(s, maxlag)
        # np.max propagates a NaN wherever it sits, and a bare `>` comparison
        # does not, so stop at the first one rather than skipping past it.
        m = 0.0
        for k in range(1, maxlag + 1):
            v = abs(rs[k])
            if np.isnan(v):
                m = v
                break
            if v > m:
                m = v
        out_max[d] = m
        # Q comes free from the ACF this loop already computes (lb_q_from_acf).
        for k in range(1, maxlag + 1):
            lbterm[k - 1] = rs[k] * rs[k] / (n - k)
        out_lb[d] = n * (n + 2) * _pairwise_sum_nb(lbterm)
        if do_prime:
            out_prime[d] = (
                _pairwise_sum_nb(rs[pr]) / pr.shape[0] - _pairwise_sum_nb(rs[comp]) / comp.shape[0]
            )
        else:
            out_prime[d] = np.nan
        out_mdl[d] = np.nan
        out_quant[d] = np.nan
        if do_steps:
            # The draw differenced and split exactly as the observation is
            # (code_increments): only consecutive cycles make a step.
            for p in range(nstep):
                inc[p] = s[hi_idx[p]] - s[lo_idx[p]]
            _clip_steps_nb(inc, inc_c, clip_sigma, clip_max_frac, kick_min_steps, kick_spread)
            if decodable:
                # A draw that folds nowhere returns 1e300: no compression, a
                # miss, recorded as one so the draw count stays honest.
                out_mdl[d] = _fold_cost_tab_nb(inc_c, tab, tcnt, tlens, tkk)[0]
            if quantisable:
                out_quant[d] = _otsu_nb(inc_c)
    return out_max, out_lb, out_prime, out_mdl, out_quant


@njit(cache=True)
def _perm4_index(w: np.ndarray) -> int:
    idx = np.array([0, 1, 2, 3])
    for a in range(3):
        for b in range(3 - a):
            if w[idx[b]] > w[idx[b + 1]] or (w[idx[b]] == w[idx[b + 1]] and idx[b] > idx[b + 1]):
                tmp = idx[b]
                idx[b] = idx[b + 1]
                idx[b + 1] = tmp
    code = 0
    facts = (6, 2, 1, 1)
    for i in range(4):
        smaller = 0
        for j in range(i + 1, 4):
            if idx[j] < idx[i]:
                smaller += 1
        code += smaller * facts[i]
    return code


@njit(cache=True)
def _perm_entropy4(x: np.ndarray, delay: int = 1) -> float:
    n = x.shape[0] - 3 * delay
    if n < 10:
        return np.nan
    counts = np.zeros(24)
    w = np.empty(4)
    for i in range(n):
        w[0] = x[i]
        w[1] = x[i + delay]
        w[2] = x[i + 2 * delay]
        w[3] = x[i + 3 * delay]
        counts[_perm4_index(w)] += 1.0
    s = 0.0
    tot = float(n)
    for c in counts:
        if c > 0.0:
            p = c / tot
            s -= p * np.log2(p)
    return s / np.log2(24.0)


@njit(cache=True)
def _slope_coeff(s2: np.ndarray, r2: np.ndarray) -> float:
    n = s2.shape[0]
    mean_s = 0.0
    mean_r = 0.0
    for i in range(n):
        mean_s += s2[i]
        mean_r += r2[i]
    mean_s /= n
    mean_r /= n
    var_s = 0.0
    cov = 0.0
    for i in range(n):
        ds = s2[i] - mean_s
        var_s += ds * ds
        cov += ds * (r2[i] - mean_r)
    if var_s <= 0.0:
        return 0.0
    return cov / var_s


@njit(cache=True)
def _jitter_boots(s2: np.ndarray, r2: np.ndarray, nboot: int, seed: int) -> np.ndarray:
    np.random.seed(seed)
    n = s2.shape[0]
    boots = np.empty(nboot)
    ss = np.empty(n)
    rr = np.empty(n)
    for b in range(nboot):
        for i in range(n):
            j = np.random.randint(0, n)
            ss[i] = s2[j]
            rr[i] = r2[j]
        boots[b] = _slope_coeff(ss, rr)
    return boots


# Median points per covered cycle needed to route a star to dense timing.
DENSE_MIN_PTS_PER_CYCLE = 30
# Points needed to fit one individual cycle. Lower than the routing threshold
# on purpose: once a star is in dense mode, a short cycle that still covers
# >= min_coverage of the phase is worth timing rather than dropping.
DENSE_MIN_PTS_PER_FIT = 20
SPARSE_MIN_PTS_PER_BLOCK = 8
# Seasonal gaps in ground-based sampling; blocks either side of one are
# independent measurements of the same O-C structure.
SPARSE_SEASON_GAP_DAYS = 60.0
ACF_P_FOLLOWUP = 1e-3
# Side-condition on the jitter test: the fitted timing variance has to stand
# clear of its own bootstrap spread, not just have a small p-value.
JITTER_MIN_SNR = 4.0
# Prewhitening floor: a secondary mode is removed before timing when its
# amplitude exceeds this many times its formal error in the residual spectrum
# (the S/N >= 4 criterion of frequency analysis), rather than a fixed fraction
# of the dominant amplitude, which left the millimagnitude non-radial modes of
# first-overtone RR Lyrae stars in the light curve (review 2, m8).
PREWHITEN_MIN_SNR = 4.0
# ... together with an absolute floor of this fraction of the dominant mode's
# peak-to-peak amplitude, ten times below the 2 per cent used before review 2
# (the non-radial modes of first-overtone RR Lyrae stars sit at ~1 per cent of
# the overtone), so that a long, clean series does not have its formally
# significant noise spikes fitted as modes.
PREWHITEN_MIN_FRACTION = 0.002


def mc_floor(draws: int | None) -> float:
    """Smallest add-one p-value obtainable from this many draws."""
    try:
        n = int(draws or 0)
    except (TypeError, ValueError):
        return 1.0
    return 1.0 / (n + 1) if n > 0 else 1.0


def draws_for_gate(gate: float, margin: float = 10.0) -> int:
    """Draws needed for an add-one p-value to resolve `gate`, with headroom.

    A gate below 1/(draws+1) can only be met by a fully saturated p-value, so
    it stops being a threshold and becomes 'beat every draw'. `margin` is how
    many distinguishable p-values to leave below the gate.
    """
    return math.ceil(margin / float(gate))


# Derived, not chosen: the jitter bootstrap carries the whole sparse-mode case,
# so its draw count has to leave room under ACF_P_FOLLOWUP.
JITTER_NBOOT = draws_for_gate(ACF_P_FOLLOWUP)
# The shuffle null must at least reach the gate.
MIN_SHUFFLE = draws_for_gate(ACF_P_FOLLOWUP, margin=1.0)
_GATE_WARNED: set[str] = set()


def below_gate(p: float | None, draws: int | None, gate: float, what: str) -> bool:
    """Apply a Monte-Carlo gate, complaining if the draws cannot resolve it."""
    floor = mc_floor(draws)
    if draws and floor >= gate and what not in _GATE_WARNED:
        _GATE_WARNED.add(what)
        log(
            f"warning: the {what} gate ({gate:g}) is below the {int(draws)}-draw "
            f"resolution floor ({floor:.2g}); only a saturated p-value can meet it"
        )
    return _p(p) < gate


TIC_NAME_RE = re.compile(r"(?i)^\s*(?:TIC\s*)?(\d+)\s*$")
ASASSN_NAME_RE = re.compile(r"(?i)^\s*ASAS-?SN\s*(\d+)\s*$")
ZTF_NAME_RE = re.compile(r"(?i)^\s*ZTF\s*(\d{9,})\s*$")  # 'ZTF <oid>'
RADEC_NAME_RE = re.compile(r"(?i)^\s*RADEC\s+([-+0-9.]+)\s+([-+0-9.]+)\s*$")
CADENCE_EXPTIME = {"short": (100.0, 150.0), "fast": (15.0, 25.0)}
DEFAULT_TARGETS = {
    "cepheid": [
        "delta Cep",
        "eta Aql",
        "zeta Gem",
        "beta Dor",
        "l Car",
        "X Cyg",
        "T Vul",
        "RT Aur",
        "SU Cas",
        "RR Lyr",
        "XZ Cyg",
        "RS Boo",
    ],
    "rrlyrae": ["RR Lyr", "XZ Cyg", "RS Boo", "SW And", "RR Gem", "TU UMa", "DH Peg", "T Sex"],
    # BW Vul: single-mode, large-amplitude, shock-driven -- the only Cepheid-like
    # relaxation oscillator in the class. The rest are classical multiperiodic
    # beta Cep stars, useful as controls for the prewhitening.
    "bcep": [
        "BW Vul",
        "beta Cep",
        "nu Eri",
        "12 Lac",
        "delta Cet",
        "xi1 CMa",
        "beta CMa",
        "V2052 Oph",
    ],
}
STAR_CLASS_HEADER_RE = re.compile(r"(?i)^\s*#\s*star[_ -]?class\s*:\s*([A-Za-z]+)")


@dataclass(frozen=True)
class StarClass:
    """Per-class analysis presets."""

    name: str
    pmin: float  # period search window (days)
    pmax: float
    max_nharm: int  # BIC-limited Fourier template order
    prewhiten: int  # independent frequencies removed before timing (0 = off)
    min_amplitude: float = 0.0  # skip stars whose template peak-to-peak is below this (rel. flux)
    # Prewhiten ground-based (sparse) light curves too? Nightly sampling puts
    # daily aliases where a second mode would be, so a class that is only
    # occasionally multiperiodic prewhitens its space photometry alone.
    prewhiten_ground: bool = True


STAR_CLASSES = {
    "cepheid": StarClass("cepheid", 0.05, 60.0, 12, 0),
    # RRc/RRd stars carry a second mode at a few per cent of the first often
    # enough that KIC 5520878 (Hippke et al. 2015) is not an oddity: sampled once
    # per cycle its beat reads as a stable two-level "word" at every convergent
    # of the period ratio (research/irrational-beating-trap.md). Up to three
    # independent modes are removed before timing, on TESS data only.
    "rrlyrae": StarClass("rrlyrae", 0.1, 2.0, 12, 3, prewhiten_ground=False),
    # 2-7 h periods; the dominant mode is often <0.05 mag so the template is kept
    # low-order and up to 5 further modes are prewhitened before the cycle fits.
    "bcep": StarClass("bcep", 0.03, 1.0, 6, 5, 0.002),
}
DEFAULT_STAR_CLASS = "cepheid"


def star_class_from_file(path: str | Path) -> str | None:
    """'# star_class: rrlyrae' header written by stage0_build_targets.py."""
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                m = STAR_CLASS_HEADER_RE.match(line)
                if m and m.group(1).lower() in STAR_CLASSES:
                    return m.group(1).lower()
    except OSError:
        pass
    return None


def resolve_star_class(name: str | None, target_file: str | None) -> StarClass:
    if name and name != "auto":
        return STAR_CLASSES[name]
    hint = star_class_from_file(target_file) if target_file else None
    return STAR_CLASSES[hint or DEFAULT_STAR_CLASS]


def resolve_cls(star_class: StarClass | str | None) -> StarClass:
    """Accept a StarClass or a class name.

    Workers are handed the dataclass itself: ProcessPoolExecutor spawns on
    Windows, so a --prewhiten override written into the module-level
    STAR_CLASSES dict does not survive the re-import in a child process and
    the worker silently used the class default. A plain name is still accepted
    so a parent running older code can drive a newer worker.
    """
    if isinstance(star_class, StarClass):
        return star_class
    return STAR_CLASSES.get(str(star_class), STAR_CLASSES[DEFAULT_STAR_CLASS])


def with_prewhiten(cls: StarClass, nmodes: int | None) -> StarClass:
    """The class presets with --prewhiten applied (identity when not given)."""
    if nmodes is None:
        return cls
    return StarClass(**{**asdict(cls), "prewhiten": nmodes})


VERBOSE = False


def log(msg: str) -> None:
    tqdm.write(str(msg))


def vlog(msg: str) -> None:
    if VERBOSE:
        tqdm.write(str(msg))


# ============================================================================
# 1. Data acquisition
# ============================================================================


@dataclass
class LightCurveData:
    t: np.ndarray
    f: np.ndarray
    e: np.ndarray
    cam: np.ndarray | None  # ASAS-SN camera / ZTF object id per point, else None
    survey: str  # "tess" | "kepler" | "asassn" | "ztf"
    time_label: str


# Space-based, near-continuous photometry: one O-C point per cycle, no nightly
# aliases, and the RR Lyrae prewhitening that ground sampling cannot support.
# Kepler is here for named stars only (`--survey kepler --targets "KIC ..."`);
# it has no stage-0 catalog, but KIC 5520878 has to be run through the same
# cuts as everything else (research/irrational-beating-trap.md).
DENSE_SURVEYS = ("tess", "kepler")


def is_dense_survey(survey: str | None) -> bool:
    return str(survey or "").lower() in DENSE_SURVEYS


def parse_tic_id(name: str) -> int | None:
    if (
        ASASSN_NAME_RE.fullmatch(name)
        or ZTF_NAME_RE.fullmatch(name)
        or RADEC_NAME_RE.fullmatch(name)
    ):
        return None
    m = TIC_NAME_RE.fullmatch(name)
    return int(m.group(1)) if m else None


def _chunks(items: Sequence[Any], n: int) -> Iterator[Sequence[Any]]:
    for i in range(0, len(items), n):
        yield items[i : i + n]


# ---- TESS -------------------------------------------------------------------

TESS_AUTHOR_FALLBACKS = ("TESS-SPOC", "QLP", "TGLC", "T16", "GSFC-ELEANOR-LITE")
TESS_PREFERRED_AUTHORS = ("SPOC", *TESS_AUTHOR_FALLBACKS)
_LK_FITS_PATH_RE = re.compile(r"((?:[A-Za-z]:)?[^\s]+?\.fits)")


def _normalize_median(f: np.ndarray, e: np.ndarray | None) -> tuple[np.ndarray, np.ndarray | None]:
    """Put a light curve on relative flux so FFI pipelines can be compared / stitched.

    A non-positive median is not a flux scale -- dividing by a negative one
    silently flips the light curve upside down -- so such a product is left
    alone rather than normalized.
    """
    med = float(np.nanmedian(f))
    if not np.isfinite(med) or med <= 0:
        return f, e
    return f / med, None if e is None else e / med


def arrays_from_lightcurve(lc: Any) -> LightCurveData:
    """Arrays from a Lightkurve product, cleaned the way our own readers clean.

    `remove_nans` drops cadences with no *flux*; it says nothing about the
    time or the uncertainty, and an HLSP is free to carry a zero or negative
    error bar on a cadence whose flux is a finite number. The weights are
    1/e^2, so a single zero is an infinity in the design matrix, and the fit
    comes back as "SVD did not converge in Linear Least Squares" from a place
    that has no idea which cadence caused it -- thirteen RR Lyrae failed this
    way, one of them carrying thirty zero errors and a flux of -9e14.
    `_clip_rel_flux` is the cut the hand-written readers already make.
    """
    lc = lc.remove_nans().remove_outliers(sigma=6)
    t = np.asarray(lc.time.value, float)
    f = np.asarray(lc.flux.value, float)
    e = np.asarray(lc.flux_err.value, float) if lc.flux_err is not None else None
    f, e = _normalize_median(f, e)
    if e is not None:
        # Mark the unusable ones rather than discarding the good ones with
        # them: a product with thirty zero error bars still has a real
        # uncertainty on its other thousand cadences, and `_clip_rel_flux`
        # drops exactly the marked rows below.
        e = np.where(np.isfinite(e) & (e > 0), e, np.nan)
    if e is None or not np.any(np.isfinite(e)):
        scatter = np.nanstd(np.diff(f)) / math.sqrt(2)
        if not np.isfinite(scatter) or scatter <= 0:
            # A perfectly flat product has no scatter to measure. Any positive
            # constant will do -- equal weights are an unweighted fit -- but
            # zero will not: that is an infinite weight, which is the bug this
            # whole path exists to avoid.
            scatter = max(abs(float(np.nanmedian(f))), 1.0) * 1e-6
        e = np.full_like(f, scatter)
    t, f, e = _clip_rel_flux(t, f, e)
    return LightCurveData(t, f, e, None, "tess", "BTJD")


def stitch_lightcurves(lcs: Any) -> LightCurveData | None:
    if not lcs:
        return None
    import lightkurve as lk

    coll = lcs if hasattr(lcs, "stitch") else lk.LightCurveCollection(list(lcs))
    if len(coll) == 0:
        return None
    return arrays_from_lightcurve(coll.stitch(corrector_func=lambda x: x.remove_nans().normalize()))


def stitch_lcds(pieces: list[LightCurveData]) -> LightCurveData | None:
    """Concatenate same-pipeline sectors after putting each on relative flux."""
    pieces = [p for p in pieces if p is not None and len(p.t)]
    if not pieces:
        return None
    flux, err = [], []
    for p in pieces:
        f, e = _normalize_median(np.asarray(p.f, float), np.asarray(p.e, float))
        flux.append(f)
        err.append(e)
    t = np.concatenate([p.t for p in pieces])
    f = np.concatenate(flux)
    e = np.concatenate(err)
    order = np.argsort(t, kind="mergesort")
    t, f, e = t[order], f[order], e[order]
    if len(t) > 1:
        keep = np.ones(len(t), dtype=bool)
        keep[1:] = np.diff(t) > 1e-8
        t, f, e = t[keep], f[keep], e[keep]
    if len(t) < 20:
        return None
    return LightCurveData(t, f, e, None, "tess", "BTJD")


def _clip_rel_flux(
    t: np.ndarray, f: np.ndarray, e: np.ndarray, nsigma: float = 6.0
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    good = np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (e > 0)
    t, f, e = t[good], f[good], e[good]
    if len(f) < 20:
        return t, f, e
    scatter = stats.median_abs_deviation(f, scale="normal")
    if not np.isfinite(scatter) or scatter <= 0:
        return t, f, e
    keep = np.abs(f - np.median(f)) < nsigma * scatter
    return t[keep], f[keep], e[keep]


def _mag_to_rel_flux(mag: Sequence[float], merr: Sequence[float]) -> tuple[np.ndarray, np.ndarray]:
    mag = np.asarray(mag, float)
    merr = np.asarray(merr, float)
    med = np.nanmedian(mag)
    f = 10 ** (-0.4 * (mag - med))
    e = f * merr * np.log(10) / 2.5
    return f, e


def _t16_best_aperture(data: Any) -> int:
    """Pick the IRM aperture with the smallest median uncertainty (1, 2 or 3)."""
    best, best_med = 2, np.inf
    for i in (1, 2, 3):
        if f"IRM{i}" not in data.names or f"IRE{i}" not in data.names:
            continue
        med = np.nanmedian(np.asarray(data[f"IRE{i}"], float))
        if np.isfinite(med) and med < best_med:
            best, best_med = i, med
    return best


def read_t16_fits(path: str | Path) -> LightCurveData | None:
    """T16 Project FFI light curves (Hartman et al. 2025). Lightkurve cannot read these.

    Use the raw IRM magnitudes: EPD/TFA high-pass out variations longer than ~1 day,
    which erases Cepheid pulsation. Times are already BTJD (BJD-2457000).
    """
    from astropy.io import fits

    with fits.open(path, memmap=False) as hdul:
        data = hdul[1].data
        if data is None or "TMID_BJD" not in data.names or "IRM1" not in data.names:
            return None
        k = _t16_best_aperture(data)
        t = np.asarray(data["TMID_BJD"], float)
        f, e = _mag_to_rel_flux(data[f"IRM{k}"], data[f"IRE{k}"])
    t, f, e = _clip_rel_flux(t, f, e)
    if len(t) < 20:
        return None
    return LightCurveData(t, f, e, None, "tess", "BTJD")


def read_tars_fits(path: str | Path) -> LightCurveData | None:
    """TARS FFI products: a LIGHTCURVE table with TIME + FLUX only."""
    from astropy.io import fits

    with fits.open(path, memmap=False) as hdul:
        data = hdul[1].data
        if data is None or "TIME" not in data.names or "FLUX" not in data.names:
            return None
        t = np.asarray(data["TIME"], float)
        f = np.asarray(data["FLUX"], float)
    good = np.isfinite(t) & np.isfinite(f)
    t, f = t[good], f[good]
    if len(t) < 20:
        return None
    med = np.median(f)
    if not np.isfinite(med) or med == 0:
        return None
    f = f / med
    e = np.full_like(f, np.nanstd(np.diff(f)) / math.sqrt(2) if len(f) > 2 else 0.01)
    t, f, e = _clip_rel_flux(t, f, e)
    if len(t) < 20:
        return None
    return LightCurveData(t, f, e, None, "tess", "BTJD")


def read_tess_fits(path: str | Path) -> LightCurveData | None:
    """Read one TESS FITS file, including HLSP formats Lightkurve rejects as 'generic'."""
    path = Path(path)
    if not path.exists() or path.stat().st_size < 2880:
        return None
    try:
        from astropy.io import fits

        with fits.open(path, memmap=False) as hdul:
            names = list(getattr(hdul[1].data, "names", None) or []) if len(hdul) > 1 else []
    except Exception:
        return None
    if "TMID_BJD" in names and "IRM1" in names:
        return read_t16_fits(path)
    if (
        names
        and "TIME" in names
        and "FLUX" in names
        and "PDCSAP_FLUX" not in names
        and "SAP_FLUX" not in names
    ):
        data = read_tars_fits(path)
        if data is not None:
            return data
    try:
        import lightkurve as lk

        return arrays_from_lightcurve(lk.read(path))
    except Exception:
        return None


def _fits_path_from_error(exc: BaseException) -> Path | None:
    text = str(exc)
    hits = [m.group(1) for m in _LK_FITS_PATH_RE.finditer(text)]
    for raw in reversed(hits):
        p = Path(raw)
        if p.suffix.lower() == ".fits":
            return p
    return None


def _lightkurve_cache_dir() -> Path:
    try:
        from lightkurve.config import get_cache_dir

        return Path(get_cache_dir())
    except Exception:
        return Path.home() / ".cache" / "lightkurve"


def _lightkurve_cache_path(row: Any) -> Path | None:
    try:
        obs_id = str(row["obs_id"])
        fname = str(row["productFilename"])
        obs_col = str(row["obs_collection"]) if "obs_collection" in row.colnames else "HLSP"
    except Exception:
        return None
    if not obs_id or not fname or fname == "None":
        return None
    return _lightkurve_cache_dir() / "mastDownload" / obs_col / obs_id / fname


def _unlink_quiet(path: Path) -> None:
    try:
        path.unlink()
    except OSError:
        pass


def _search_tess(target: str, author: str | None = None, cadence: str | None = None) -> Any:
    import lightkurve as lk

    kw = {"mission": "TESS"}
    if author:
        kw["author"] = author
    if cadence:
        kw["cadence"] = cadence
    return lk.search_lightcurve(target, **kw)


def _prefer_tess_authors(sr: Any) -> Any:
    """Keep a single TESS pipeline. Mixing SPOC/QLP/TGLC/T16/Eleanor of the same
    sector concatenates incompatible flux scales and kills the periodogram."""
    if sr is None or len(sr) == 0 or not hasattr(sr, "author"):
        return sr
    authors = np.asarray([str(a).upper() for a in sr.author])
    for pref in TESS_PREFERRED_AUTHORS:
        mask = authors == pref
        if np.any(mask):
            return sr[mask]
    return sr[authors == authors[0]]


def _tess_candidates(sr: Any, author: str, cadence: str) -> Any:
    """Pick one pipeline out of a single search result.

    lightkurve re-queries MAST on every search_lightcurve() call, so the old
    author/cadence fallback chain cost up to seven round trips per star. One
    unfiltered search plus this local preference order is equivalent: the
    requested author at the requested cadence, that author at any cadence,
    then each FFI pipeline in turn.
    """
    if sr is None or len(sr) == 0:
        return None
    try:
        authors = np.asarray([str(a).upper() for a in sr.table["author"]])
        exptime = np.asarray(sr.table["exptime"], dtype=float)
    except Exception:
        return _prefer_tess_authors(sr)
    lo, hi = CADENCE_EXPTIME.get(cadence, CADENCE_EXPTIME["short"])
    at_cadence = (exptime >= lo) & (exptime <= hi)
    want = str(author or "SPOC").upper()
    order = [(authors == want) & at_cadence, authors == want]
    order += [authors == auth.upper() for auth in TESS_AUTHOR_FALLBACKS]
    for mask in order:
        if mask.any():
            return sr[mask]
    return _prefer_tess_authors(sr)


class TessProductsUnavailable(RuntimeError):
    """Every TESS product this star has failed to download.

    Not the same thing as TESS having no product for it: most Gaia variables
    are simply not 2-minute targets, and that is a fact about the mission,
    while this is a fact about the afternoon.
    """


def _download_one_search_row(sr_one: Any) -> LightCurveData | None:
    """Download one MAST product. Read T16/TARS ourselves if Lightkurve cannot.

    None means the product came down and held nothing usable. A product that
    could not be fetched at all raises instead, so `download_tess` can tell
    the two apart -- a star whose every product failed is a retryable failure,
    and recording it as "no data" put it beyond the reach of --retry-errors
    and made a MAST outage look like an absence of data.
    """
    row = (
        sr_one.table[0]
        if getattr(sr_one, "table", None) is not None and len(sr_one.table)
        else None
    )
    cached = _lightkurve_cache_path(row) if row is not None else None
    if cached is not None and cached.exists():
        data = read_tess_fits(cached)
        if data is not None:
            return data
        _unlink_quiet(cached)
    try:
        lc = sr_one.download()
    except Exception as ex:
        path = _fits_path_from_error(ex) or cached
        if path is not None and path.exists():
            data = read_tess_fits(path)
            if data is not None:
                return data
            if "corrupt" in str(ex).lower() or "interrupted" in str(ex).lower():
                _unlink_quiet(path)
                try:
                    lc = sr_one.download()
                except Exception as ex2:
                    path2 = _fits_path_from_error(ex2)
                    if path2 is not None and path2.exists():
                        data = read_tess_fits(path2)
                        if data is not None:
                            return data
                    raise
                return None if lc is None else arrays_from_lightcurve(lc)
        raise
    return None if lc is None else arrays_from_lightcurve(lc)


def load_lightcurve_files(paths: Sequence[str | Path]) -> LightCurveData | None:
    pieces = []
    for p in paths:
        data = read_tess_fits(p)
        if data is not None:
            pieces.append(data)
    return stitch_lcds(pieces)


def download_tess(
    target: str, cadence: str = "short", author: str = "SPOC"
) -> LightCurveData | None:
    """SPOC at the requested cadence, then SPOC any, then one FFI pipeline.

    Preference after SPOC is TESS-SPOC, QLP, TGLC, T16, Eleanor. Never mix
    authors: three Sector-8 reductions of one RRab (TGLC+T16+Eleanor) used to
    be concatenated, which erased the pulsation. Products are fetched one file
    at a time so a single unreadable HLSP does not abort the star.

    One MAST search; the author / cadence preference is applied locally.
    """
    sr = _tess_candidates(_search_tess(target), author=author, cadence=cadence)
    if sr is None or len(sr) == 0:
        return None  # the mission has no product for this star
    pieces, failures = [], []
    for i in range(len(sr)):
        try:
            data = _download_one_search_row(sr[i])
        except Exception as ex:
            failures.append(ex)
            continue
        if data is not None:
            pieces.append(data)
    if pieces:
        return stitch_lcds(pieces)  # one unreadable product does not cost the star
    if failures:
        raise TessProductsUnavailable(
            f"all {len(sr)} TESS product(s) for {target} failed to download; "
            f"first was {type(failures[0]).__name__}: {failures[0]}"
        )
    return None


# stage-1 --cadence -> Kepler cadence. 'short' and 'fast' both mean the 1-min
# short-cadence months (Kepler had nothing faster); 'long' asks for the 30-min
# long cadence outright, which for a P < 0.7 d star routes to sparse timing.
KEPLER_CADENCE = {"short": "short", "fast": "short", "long": "long"}


def download_kepler(target: str, cadence: str = "short") -> LightCurveData | None:
    """Kepler prime-mission light curves of one KIC, stitched across quarters.

    Short cadence where the star had it, else long cadence (or long cadence
    outright with cadence='long'). One author -- the Kepler pipeline's PDCSAP
    flux -- so quarters go on relative flux and concatenate the way TESS
    sectors do. Times stay in BKJD; only differences matter downstream.
    """
    import lightkurve as lk

    want = KEPLER_CADENCE.get(str(cadence).lower(), "short")
    sr = None
    for cad in [want] if want == "long" else ["short", "long"]:
        try:
            sr = lk.search_lightcurve(target, mission="Kepler", author="Kepler", cadence=cad)
        except Exception as ex:
            log(f"  Kepler search failed for {target}: {ex}")
            return None
        if sr is not None and len(sr):
            break
    if sr is None or len(sr) == 0:
        return None
    pieces = []
    for i in range(len(sr)):
        try:
            lc = sr[i].download()
        except Exception as ex:
            vlog(f"  {target}: quarter {i} not read ({ex})")
            continue
        if lc is not None:
            pieces.append(arrays_from_lightcurve(lc))
    data = stitch_lcds(pieces)
    if data is None:
        return None
    return LightCurveData(data.t, data.f, data.e, None, "kepler", "BKJD")


def resolve_target_tics(targets: Sequence[str]) -> dict[str, int]:
    from astroquery.mast import Catalogs

    out, unresolved = {}, []
    for name in targets:
        tic = parse_tic_id(name)
        (out.__setitem__(name, tic) if tic is not None else unresolved.append(name))
    for name in unresolved:
        try:
            tab = Catalogs.query_object(name, catalog="TIC")
            out[name] = int(tab[0]["ID"]) if len(tab) else None
        except Exception:
            out[name] = None
        vlog(f"  resolved {name} -> TIC {out[name]}")
    return out


def _query_tess_obs(tic_ids: Sequence[int], cadence: str, strict: bool) -> Any:
    from astropy.table import vstack
    from astroquery.mast import Observations

    tables = []
    for chunk in tqdm(
        list(_chunks([str(t) for t in tic_ids], 100)), desc="MAST query", unit="chunk", leave=False
    ):
        kw = {"obs_collection": "TESS", "dataproduct_type": "timeseries", "target_name": chunk}
        if strict:
            kw["provenance_name"] = "SPOC"
            kw["t_exptime"] = list(CADENCE_EXPTIME.get(cadence, CADENCE_EXPTIME["short"]))
        try:
            tab = Observations.query_criteria(**kw)
        except Exception as ex:
            log(f"  MAST query failed: {ex}")
            continue
        if tab is not None and len(tab):
            tables.append(tab)
    return vstack(tables, metadata_conflicts="silent") if tables else None


def _lc_products_for_obs(obs: Any) -> Any:
    from astropy.table import vstack
    from astroquery.mast import Observations

    tables = []
    for chunk in tqdm(
        list(_chunks(list(obs["obsid"]), 40)), desc="LC products", unit="chunk", leave=False
    ):
        try:
            tab = Observations.get_product_list(chunk)
        except Exception as ex:
            log(f"  product list failed: {ex}")
            continue
        if tab is None or len(tab) == 0:
            continue
        if "productSubGroupDescription" in tab.colnames:
            tab = tab[tab["productSubGroupDescription"] == "LC"]
        if "productFilename" in tab.colnames:
            tab = tab[[str(fn).lower().endswith(".fits") for fn in tab["productFilename"]]]
        if len(tab):
            tables.append(tab)
    return vstack(tables, metadata_conflicts="silent") if tables else None


def _download_one_product(row: Any, base_dir: Path) -> Path | None:
    from astroquery.mast import Observations

    local_dir = os.path.join(
        base_dir, str(row.get("obs_collection") or "TESS"), str(row.get("obs_id") or "unknown")
    )
    os.makedirs(local_dir, exist_ok=True)
    local_path = os.path.join(local_dir, os.path.basename(str(row["productFilename"])))
    status, msg, _ = Observations.download_file(
        row["dataURI"], local_path=local_path, cache=True, verbose=False
    )
    if status == "ERROR":
        raise RuntimeError(msg or status)
    if read_tess_fits(local_path) is None:
        _unlink_quiet(Path(local_path))
        status, msg, _ = Observations.download_file(
            row["dataURI"], local_path=local_path, cache=False, verbose=False
        )
        if status == "ERROR" or read_tess_fits(local_path) is None:
            raise RuntimeError(msg or "unreadable light-curve FITS after re-download")
    return local_path, status


def prefetch_sparse(
    targets: Sequence[str],
    survey: str,
    jobs: int = 8,
    policy: lc_cache.CachePolicy | None = None,
    hints: dict[str, dict] | None = None,
) -> dict[str, int]:
    """Fill the ASAS-SN / ZTF cache for `targets`, analysing nothing.

    The sparse counterpart of `prefetch_lightcurves`, and the thing that makes
    a sparse injection grid possible: `inject.hosts_from_store` will only draw
    a host whose photometry is already on disk, because re-downloading ~25 s
    of it per injection is not a grid. Downloading is I/O, so threads; the
    analysis this deliberately skips is what needs processes.

    Already-cached targets are not re-fetched unless `policy.refresh` says so.
    A target the survey has nothing for is counted and *not* stored -- see the
    note on negative caching in `lc_cache`.
    """
    if survey not in lc_cache.SURVEYS:
        raise ValueError(f"{survey!r} has no raw-row cache; one of {lc_cache.SURVEYS}")
    hint_of = hints or {}

    def _fetch_ztf(name: str) -> pd.DataFrame | None:
        h = hint_of.get(name) or {}
        return _ztf_raw(name, pos=h.get("pos"), oids=h.get("oids"))

    fetch = _fetch_ztf if survey == "ztf" else _asassn_raw
    params = {"merge_fields": True} if survey == "ztf" else None
    pol = policy or lc_cache.DEFAULT_POLICY
    names = list(dict.fromkeys(str(t).strip() for t in targets if str(t).strip()))
    todo = [n for n in names if not lc_cache.has(survey, n, params, policy=pol)]
    counts = {"targets": len(names), "cached": len(names) - len(todo), "fetched": 0}
    counts.update(empty=0, failed=0)
    log(f"Prefetching {len(todo)} {survey} light curve(s); {counts['cached']} already cached")
    if not todo:
        return counts
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        futs = {pool.submit(fetch, n): n for n in todo}
        for fut in tqdm(as_completed(futs), total=len(futs), desc=f"{survey} prefetch", unit="lc"):
            name = futs[fut]
            try:
                df = fut.result()
            except Exception as ex:
                vlog(f"  {name}: {ex}")
                counts["failed"] += 1
                continue
            if df is None or not len(df):
                counts["empty"] += 1
                continue
            if lc_cache.save(survey, name, df, params, policy=pol) is None:
                counts["failed"] += 1
            else:
                counts["fetched"] += 1
    log(f"  {counts['fetched']} cached, {counts['empty']} with no data, {counts['failed']} failed")
    return counts


def prefetch_lightcurves(
    targets: Sequence[str], cadence: str, download_dir: str | Path, download_jobs: int
) -> dict[str, list[str]]:
    """One MAST search + parallel LC download. Returns {target: [fits paths]}."""
    from astroquery.mast import Observations

    download_dir = Path(download_dir)
    download_dir.mkdir(parents=True, exist_ok=True)
    base_dir = str(download_dir / "mastDownload")
    log("Resolving target names to TIC IDs ...")
    name_to_tic = resolve_target_tics(targets)
    tic_ids = sorted({t for t in name_to_tic.values() if t is not None})
    if not tic_ids:
        return {n: [] for n in targets}
    try:
        Observations.enable_cloud_dataset(verbose=False)
    except Exception:
        pass
    log(f"Querying MAST for {len(tic_ids)} TIC IDs ...")
    obs = _query_tess_obs(tic_ids, cadence, strict=True)
    found = (
        {int(x) for x in obs["target_name"] if str(x).isdigit()}
        if obs is not None and "target_name" in obs.colnames
        else set()
    )
    missing = [t for t in tic_ids if t not in found]
    if missing:
        extra = _query_tess_obs(missing, cadence, strict=False)
        if extra is not None:
            from astropy.table import vstack

            obs = extra if obs is None else vstack([obs, extra], metadata_conflicts="silent")
    if obs is None or len(obs) == 0:
        return {n: [] for n in targets}
    if "obsid" in obs.colnames:
        _, uniq = np.unique(np.asarray(obs["obsid"], dtype=str), return_index=True)
        obs = obs[np.sort(uniq)]
    products = _lc_products_for_obs(obs)
    if products is None or len(products) == 0:
        return {n: [] for n in targets}
    rows = [
        {
            "obs_collection": str(r["obs_collection"])
            if "obs_collection" in products.colnames
            else "TESS",
            "obs_id": str(r["obs_id"]),
            "productFilename": str(r["productFilename"]),
            "dataURI": str(r["dataURI"]),
        }
        for r in products
    ]
    log(f"  downloading {len(rows)} light-curve files with {download_jobs} thread(s) ...")
    tic_paths = {t: [] for t in tic_ids}
    fname_tic = re.compile(r"-(\d{16})-")
    with ThreadPoolExecutor(max_workers=max(1, download_jobs)) as pool:
        futs = [pool.submit(_download_one_product, row, base_dir) for row in rows]
        for fut in tqdm(as_completed(futs), total=len(futs), desc="download", unit="file"):
            try:
                path, _ = fut.result()
            except Exception as ex:
                log(f"  download error: {ex}")
                continue
            m = fname_tic.search(Path(path).name)
            if m:
                tic_paths.setdefault(int(m.group(1)), []).append(path)
    index = {
        n: (
            sorted(set(tic_paths.get(name_to_tic.get(n), [])))
            if name_to_tic.get(n) is not None
            else []
        )
        for n in targets
    }
    log(f"  local light curves for {sum(1 for p in index.values() if p)} / {len(targets)} targets")
    (download_dir / "index.json").write_text(json.dumps(index, indent=1), encoding="utf-8")
    return index


# One --band flag serves two surveys that want different answers, so its
# default is a sentinel resolved per survey rather than a literal band.
# ASAS-SN g and V are different *epochs* of the survey (V 2012-18, g 2018+)
# and mixing them is meaningless; ZTF g and r are simultaneous coverage of the
# same star, worth ~2.5x the epochs merged once camera_scales handles the
# amplitude difference.
AUTO_BAND = "auto"
DEFAULT_BAND = {"asassn": "g", "ztf": "gr"}


def resolve_band(survey: str, band: str | None) -> str:
    """The effective band for a survey. 'auto' / None take the survey default."""
    if band in (None, AUTO_BAND):
        return DEFAULT_BAND.get(survey, "g")
    return band


# ---- ASAS-SN ----------------------------------------------------------------


def _asassn_raw(target: str) -> pd.DataFrame | None:
    """Sky Patrol photometry for one target, exactly as the client returns it.

    The network half of `download_asassn`, split out so `lc_cache` has
    something to store. Nothing here cuts rows: band, quality and the outlier
    clip are read-side concerns, and a cached blob has to serve every band.
    """
    from pyasassn.client import SkyPatrolClient

    client = SkyPatrolClient(verbose=False)
    m_as = ASASSN_NAME_RE.fullmatch(target)
    tic = parse_tic_id(target)
    if m_as:
        lcc = client.query_list(
            [int(m_as.group(1))], catalog="stellar_main", id_col="asas_sn_id", download=True
        )
    elif tic is not None:
        # query the stellar_main catalog by TIC id, then pull curves by asas_sn_id
        df_idx = client.adql_query(
            f"SELECT asas_sn_id, ra_deg, dec_deg FROM stellar_main WHERE tic_id = {tic}"
        )
        if df_idx is None or len(df_idx) == 0:
            return None
        lcc = client.query_list(
            list(df_idx["asas_sn_id"]), catalog="stellar_main", id_col="asas_sn_id", download=True
        )
    else:
        lcc = client.simbad_lookup(target, download=True)
    if lcc is None or len(lcc) == 0:
        return None
    # len(lcc) counts catalog rows, so a collection carrying no photometry at
    # all still reaches here with an empty (column-less) frame. That is no
    # data, not a KeyError on mag_err.
    df = lcc.data
    if df is None or df.empty:
        return None
    return df


def download_asassn(
    target: str,
    band: str = AUTO_BAND,
    include_poor: bool = False,
    *,
    cache: lc_cache.CachePolicy | None = None,
) -> LightCurveData | None:
    """ASAS-SN Sky Patrol light curve for a SIMBAD-resolvable name (or 'TIC N').

    Raw rows come from `lc_cache` when it has a fresh entry and from Sky
    Patrol otherwise; every cut below runs either way, so what the cache
    changes is the network bill and never the light curve.
    """
    band = resolve_band("asassn", band)
    df = lc_cache.load("asassn", target, policy=cache)
    if df is None:
        df = _asassn_raw(target)
        if df is None:
            return None
        lc_cache.save("asassn", target, df, policy=cache)
    df = df.copy()
    missing = [c for c in ("jd", "mag", "mag_err") if c not in df.columns]
    if missing:
        raise RuntimeError(
            "ASAS-SN light curve missing "
            + ", ".join(missing)
            + "; got "
            + ", ".join(str(c) for c in df.columns)
        )
    if not include_poor and "quality" in df:
        df = df[df["quality"] == "G"]
    df = df[
        np.isfinite(df["jd"])  # see the note on hjd in download_ztf
        & (df["mag_err"] < 99)
        & np.isfinite(df["mag"])
        & (df["mag_err"] > 0)
    ]
    if band in ("g", "V") and "phot_filter" in df:
        df = df[df["phot_filter"] == band]
    if len(df) < 50:
        return None
    df = df.sort_values("jd")
    t = df["jd"].to_numpy(float) - 2450000.0
    mag, merr = df["mag"].to_numpy(float), df["mag_err"].to_numpy(float)
    cam = df["camera"].astype(str).to_numpy() if "camera" in df else None
    f = 10 ** (-0.4 * (mag - np.median(mag)))
    e = f * merr * np.log(10) / 2.5
    good = np.abs(f - np.median(f)) < 8 * stats.median_abs_deviation(f, scale="normal")
    return LightCurveData(
        t[good],
        f[good],
        e[good],
        None if cam is None else cam[good],
        "asassn",
        f"JD-2450000 ({band})",
    )


ZTF_LC_URL = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"
# Which filtercodes survive the local cut -- and, since `_ztf_raw` stopped
# sending BANDNAME so that one cached download serves every band, the only
# thing that decides the band at all. 'all' keeps whatever came back.
ZTF_KEEP_BANDS = {"g": ("g",), "r": ("r",), "i": ("i",), "V": ("g",), "gr": ("g", "r"), "all": ()}
ZTF_MATCH_RADIUS_DEG = 1.5 / 3600.0


def parse_ztf_oid(name: str) -> int | None:
    m = ZTF_NAME_RE.fullmatch(name)
    return int(m.group(1)) if m else None


_ZTF_FAULT_RE = re.compile(r'value="ERROR">(.*?)</INFO>', re.S)


def _ztf_fault(text: str) -> str:
    """The IRSA UsageFault message out of an error VOTable, for the error row."""
    m = _ZTF_FAULT_RE.search(text or "")
    raw = m.group(1) if m else (text or "")
    return re.sub(r"\s+", " ", unescape(raw)).strip()[:200]


class ZtfUnavailable(RuntimeError):
    """IRSA's light-curve service could not be reached.

    Deliberately not the same thing as the service answering that it has
    nothing for this target: one is a star with no ZTF photometry, the other
    is a star whose photometry we failed to ask for.
    """


# Waits between attempts on a ZTF call that could not be answered. Long enough
# that a short throttle is waited out rather than spent, and short enough that
# a service which is really down does not hold a worker for the afternoon.
ZTF_BACKOFF_S = (5.0, 15.0, 45.0, 120.0)
# Once a call has used its whole budget, the service is not having a blip.
# Every ZTF fetch in this process fails immediately for this long instead of
# each one paying the full budget again -- which both ends the run promptly
# and stops eight workers hammering something that just asked them to stop.
ZTF_COOLDOWN_S = 300.0
_ztf_cooldown_until = 0.0


def _ztf_unavailable(reason: str) -> ZtfUnavailable:
    """Start the cool-down and return the error to raise."""
    global _ztf_cooldown_until
    _ztf_cooldown_until = time.monotonic() + ZTF_COOLDOWN_S
    return ZtfUnavailable(
        f"IRSA ZTF light-curve service unreachable ({reason}); "
        f"pausing ZTF fetches for {ZTF_COOLDOWN_S:.0f}s"
    )


def _ztf_wait(attempt: int, retry_after: str | None = None) -> None:
    """Sleep before the next attempt, honouring Retry-After when it is sane."""
    delay = ZTF_BACKOFF_S[min(attempt, len(ZTF_BACKOFF_S) - 1)]
    if retry_after:
        try:
            delay = max(delay, min(float(retry_after), ZTF_BACKOFF_S[-1]))
        except (TypeError, ValueError):
            pass
    time.sleep(delay)


def _ztf_get(params: dict) -> pd.DataFrame | None:
    """One call to the IRSA ZTF light-curve API; a DataFrame, or None for
    "the survey has nothing here".

    A 4xx is raised, not swallowed: an ID / NOBS_MIN parameter clash used to
    come back as None and get stored as "no data", which hid the fact that
    every by-oid query was failing.

    Being unable to reach the service at all is raised for the same reason,
    and it is the more dangerous of the two because of how a throttle
    presents. IRSA does not answer a rate-limited caller with 429; it accepts
    the request and never replies, so it arrives as a read timeout. Returning
    None for that wrote "no data" onto stars with years of photometry, one
    quiet row at a time, and nothing in the run said the survey had stopped
    talking. Now the attempt budget waits out a blip, and anything longer
    raises `ZtfUnavailable` -- which `process_target` records as a download
    failure, so `--retry-errors` finds it again.
    """
    import io

    import requests

    if time.monotonic() < _ztf_cooldown_until:
        raise ZtfUnavailable("ZTF fetches are in cool-down after an earlier failure")
    params = {"FORMAT": "csv", **params}
    last = "no attempt made"
    for attempt in range(len(ZTF_BACKOFF_S) + 1):
        final = attempt == len(ZTF_BACKOFF_S)
        try:
            r = requests.get(ZTF_LC_URL, params=params, timeout=NETWORK_TIMEOUT_S)
            if r.status_code == 200:
                if not r.text.strip():
                    return None  # answered, and the answer is "nothing here"
                df = pd.read_csv(io.StringIO(r.text))
                return df if len(df) else None
            if 400 <= r.status_code < 500 and r.status_code != 429:
                raise RuntimeError(f"ZTF API {r.status_code}: {_ztf_fault(r.text)}")
            # 429, 503, and every other 5xx: the service is refusing rather
            # than answering, which is exactly what waiting is for.
            last = f"HTTP {r.status_code}"
            if not final:
                _ztf_wait(attempt, r.headers.get("Retry-After"))
        except RuntimeError:
            raise
        except Exception as ex:
            last = f"{type(ex).__name__}: {ex}"
            if not final:
                _ztf_wait(attempt)
    raise _ztf_unavailable(last)


def _resolve_position(target: str) -> tuple[float, float] | None:
    """(ra, dec) in degrees for 'RADEC ra dec', 'TIC N', or a SIMBAD name."""
    m = RADEC_NAME_RE.fullmatch(target)
    if m:
        return float(m.group(1)), float(m.group(2))
    tic = parse_tic_id(target)
    if tic is not None:
        from astroquery.mast import Catalogs

        tab = Catalogs.query_criteria(catalog="TIC", ID=tic)
        if len(tab):
            return float(tab[0]["ra"]), float(tab[0]["dec"])
        return None
    from astroquery.simbad import Simbad

    tab = Simbad.query_object(target)
    if tab is None or len(tab) == 0:
        return None
    ra_col = "ra" if "ra" in tab.colnames else "RA"
    dec_col = "dec" if "dec" in tab.colnames else "DEC"
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    ra, dec = str(tab[0][ra_col]), str(tab[0][dec_col])
    if ":" in ra or " " in ra.strip():
        c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    else:
        c = SkyCoord(float(ra), float(dec), unit="deg")
    return float(c.ra.deg), float(c.dec.deg)


def _ztf_raw(
    target: str,
    merge_fields: bool = True,
    pos: tuple[float, float] | None = None,
    oids: Sequence[int] | None = None,
) -> pd.DataFrame | None:
    """IRSA epochs for one target, exactly as the API returned them.

    The network half of `download_ztf`. `BANDNAME` is deliberately *not* sent:
    it is the only fetch parameter the requested band used to change, and
    sending it would fork the cache into a copy of the same star per band.
    Dropping it can only widen the response -- `keep_bands` below cuts to the
    requested filters either way -- and it loosens `NOBS_MIN` from 20 epochs
    in one band to 20 across all of them, which the local `min_epochs` cut of
    50 more than covers.

    So `merge_fields` is the only thing here the caller can vary, and the only
    thing in the cache key besides the target -- `pos` changes only how the
    cone's centre is arrived at, not what is asked for.

    Three ways in, cheapest first, and the two hints come from
    `store.ztf_fetch_hints`.

    `oids` is `stage0_targets.ztf_oids_all`: the overlapping-field objects
    stage 0 already resolved. Asking for each by id is what the cone search
    is *for*, and the cone is the expensive call -- 8 to 60 s against about a
    second per id, because IRSA resolves a position against the whole object
    table. On the stars this was checked against, the union of the per-id
    fetches was the cone's answer exactly: same oid set, same row count.

    `pos` is the catalog position, for a star with no stored oid list. It
    skips only the by-oid query that the old path ran to take the median of
    its own astrometry and then threw away with every row it came with.

    With neither, both calls run as before. A cone centred on Gaia's position
    is not centred where ZTF's astrometry puts the star, so a hinted result
    is checked for the object that was asked for and anything else falls back
    rather than quietly returning a different star's oids.

    The oid list is stage 0's view of the data release it was built against.
    A cone search sees whatever IRSA serves now, so an object added to a
    later release is one this path will not pick up; `--lc-cache-refresh`
    with no hints is the way to re-derive it.
    """
    epoch_only = {"BAD_CATFLAGS_MASK": 32768}
    pos_base = {**epoch_only, "NOBS_MIN": 20}
    oid = parse_ztf_oid(target)
    if oid is not None:
        if merge_fields and oids:
            want = list(dict.fromkeys([*oids, oid]))
            parts = [_ztf_get({"ID": int(o), **epoch_only}) for o in want]
            parts = [p for p in parts if p is not None and len(p)]
            if parts:
                return pd.concat(parts, ignore_index=True)
        if merge_fields and pos is not None:
            more = _ztf_get(
                {"POS": f"CIRCLE {pos[0]:.6f} {pos[1]:.6f} {ZTF_MATCH_RADIUS_DEG:.6f}", **pos_base}
            )
            if more is not None and "oid" in more and (more["oid"] == oid).any():
                return more
        df = _ztf_get({"ID": oid, **epoch_only})
        if df is None:
            return None
        if merge_fields:
            ra, dec = float(df["ra"].median()), float(df["dec"].median())
            more = _ztf_get(
                {"POS": f"CIRCLE {ra:.6f} {dec:.6f} {ZTF_MATCH_RADIUS_DEG:.6f}", **pos_base}
            )
            if more is not None:
                df = more
        return df
    pos = _resolve_position(target)
    if pos is None:
        return None
    return _ztf_get(
        {"POS": f"CIRCLE {pos[0]:.6f} {pos[1]:.6f} {ZTF_MATCH_RADIUS_DEG:.6f}", **pos_base}
    )


def download_ztf(
    target: str,
    band: str = AUTO_BAND,
    merge_fields: bool = True,
    min_epochs: int = 50,
    *,
    cache: lc_cache.CachePolicy | None = None,
    hint: dict | None = None,
) -> LightCurveData | None:
    """ZTF DR light curve from IRSA for 'ZTF <oid>', 'TIC N', 'RADEC ra dec' or a
    SIMBAD name. Every ZTF object id is one field/quadrant/filter, so a star that
    falls in overlapping fields has several oids; with merge_fields we pull all of
    them within 1.5" and treat each oid as a 'camera' so its zero-point offset is
    fit jointly with the template (same machinery as ASAS-SN cameras).

    band='gr' -- the ZTF default -- keeps both filters. An oid is single-band
    by construction, so the
    per-oid 'camera' already separates them, and each oid is median-normalised
    here so the two filters arrive on a common level rather than ~1 mag apart.
    What the offsets cannot absorb is that a pulsator's amplitude is smaller in
    r than in g -- see camera_scales, which is what makes merging safe.

    Raw epochs come from `lc_cache` when it has a fresh entry; the band cut and
    everything after it run either way.
    """
    band = resolve_band("ztf", band)
    keep_bands = ZTF_KEEP_BANDS.get(band, ("g",))
    params = {"merge_fields": bool(merge_fields)}
    df = lc_cache.load("ztf", target, params, policy=cache)
    if df is None:
        df = _ztf_raw(
            target,
            merge_fields=merge_fields,
            pos=(hint or {}).get("pos"),
            oids=(hint or {}).get("oids"),
        )
        if df is None:
            return None
        lc_cache.save("ztf", target, df, params, policy=cache)
    # clean: catflags==0 (best quality), finite photometry, single band
    # hjd as well as the photometry: one DR24 object carries a single NaN
    # epoch time, and one of those makes the baseline NaN, which makes
    # astropy's frequency grid NaN, which fails the star deep inside
    # `autofrequency` with "cannot convert float NaN to integer".
    df = df[
        np.isfinite(df["hjd"])
        & np.isfinite(df["mag"])
        & np.isfinite(df["magerr"])
        & (df["magerr"] > 0)
    ]
    if "catflags" in df:
        df = df[df["catflags"] == 0]
    if keep_bands and "filtercode" in df:
        df = df[df["filtercode"].isin([f"z{b}" for b in keep_bands])]
    if len(df) < min_epochs:
        return None
    df = df.sort_values("hjd")
    t = df["hjd"].to_numpy(float) - 2450000.0  # HJD keeps timing consistent with ASAS-SN
    mag, merr = df["mag"].to_numpy(float), df["magerr"].to_numpy(float)
    cam = df["oid"].astype(str).to_numpy()  # one 'camera' per field/quadrant object
    # Per-oid, not global: g and r sit ~1 mag apart, and a global median would
    # put a whole filter's worth of points at a flux level the template has to
    # absorb as a huge offset before it can see the pulsation at all.
    ref = np.zeros_like(mag)
    for c in np.unique(cam):
        m = cam == c
        ref[m] = np.median(mag[m])
    f = 10 ** (-0.4 * (mag - ref))
    e = f * merr * np.log(10) / 2.5
    good = np.abs(f - np.median(f)) < 8 * stats.median_abs_deviation(f, scale="normal")
    label = "+".join(keep_bands) if keep_bands else "all"
    return LightCurveData(t[good], f[good], e[good], cam[good], "ztf", f"HJD-2450000 (ZTF {label})")


# ============================================================================
# 2. Period + Fourier template (with per-camera offsets)
# ============================================================================

# Press & Rybicki 1989; astropy 8's default 'lra' NUFFT is more accurate but
# much slower on long TESS series, and we always refine the peak exactly.
_FAST_LS_KWDS = {"algorithm": "fasper"}


def _parabolic_argmax(freq: np.ndarray, power: np.ndarray, i: int) -> float:
    """Sub-bin frequency of a uniformly spaced 3-sample peak."""
    if i <= 0 or i >= len(power) - 1:
        return float(freq[i])
    y0, y1, y2 = (float(power[i + d]) for d in (-1, 0, 1))
    denom = y0 - 2.0 * y1 + y2
    if denom == 0.0 or not np.isfinite(denom):
        return float(freq[i])
    delta = 0.5 * (y0 - y2) / denom
    if abs(delta) > 1.0:
        return float(freq[i])
    return float(freq[i] + delta * float(freq[i] - freq[i - 1]))


def _ls_spectrum(
    ls: Any, fmin: float, fmax: float, samples_per_peak: int = 8
) -> tuple[np.ndarray, np.ndarray]:
    """O(N log N) periodogram on an automatic frequency grid."""
    kw = {
        "minimum_frequency": float(fmin),
        "maximum_frequency": float(fmax),
        "samples_per_peak": samples_per_peak,
    }
    try:
        freq, power = ls.autopower(method="fast", method_kwds=_FAST_LS_KWDS, **kw)
    except Exception:
        freq, power = ls.autopower(method="auto", **kw)
    return np.asarray(freq, dtype=float), np.asarray(power, dtype=float)


def _ls_power_grid(ls: Any, fg: np.ndarray) -> np.ndarray:
    """LS power at `fg`: the compiled grid kernel where it applies, astropy else.

    `_ls_power_grid_nb` reproduces exactly one configuration -- one term, a
    floating mean, centred data, standard normalization, plain arrays, a
    uniform grid -- so anything else goes back to astropy rather than quietly
    returning a different periodogram. Without numba the kernel is a Python
    loop over every (frequency, point) pair and astropy is the faster of the
    two by orders of magnitude, so that falls back as well.
    """
    step = np.diff(fg) if fg.size >= 2 else np.zeros(0)
    if not (
        HAVE_NUMBA
        and fg.size >= 2
        # The rotor is one fixed turn per step, so a grid that is not evenly
        # spaced would silently be read as one that is.
        and bool(np.all(np.abs(step - step[0]) <= 1e-6 * abs(step[0])))
        and getattr(ls, "fit_mean", False)
        and getattr(ls, "center_data", False)
        and getattr(ls, "normalization", "") == "standard"
        and getattr(ls, "nterms", 1) == 1
        and getattr(ls.t, "unit", None) is None
        and getattr(ls.y, "unit", None) is None
    ):
        return np.asarray(ls.power(fg, method="cython"), dtype=float)
    y = np.asarray(ls.y, dtype=float)
    w = np.ones_like(y) if ls.dy is None else np.asarray(ls.dy, dtype=float) ** -2.0
    if not np.all(np.isfinite(w)):  # a zero error bar is astropy's to refuse
        return np.asarray(ls.power(fg, method="cython"), dtype=float)
    # astropy centres the data first; with fit_mean that is a simplification
    # rather than a change of answer, but the kernel assumes it has happened.
    y = np.ascontiguousarray(y - np.dot(w, y) / w.sum())
    return _ls_power_grid_nb(
        np.ascontiguousarray(ls.t, dtype=float),
        y,
        np.ascontiguousarray(w),
        np.ascontiguousarray(fg, dtype=float),
        LS_RESEED,
    )


def _refine_freq(ls: Any, fpk: float, rel: float = 2e-3, n: int = 801) -> float:
    """Exact Lomb-Scargle on a tight window around a coarse peak."""
    lo, hi = fpk * (1.0 - rel), fpk * (1.0 + rel)
    if not (0.0 < lo < hi):
        return float(fpk)
    fg = np.linspace(lo, hi, int(n))
    p = _ls_power_grid(ls, fg)
    if not np.any(np.isfinite(p)):
        return float(fpk)  # a zero-variance series has no periodogram at all
    i = int(np.nanargmax(p))
    return _parabolic_argmax(fg, p, i)


def _peaks_desc(
    freq: np.ndarray, power: np.ndarray, min_sep: float, nkeep: int = 24
) -> list[tuple[float, float]]:
    """Strongest local maxima, greedily spaced by min_sep (one bin per peak)."""
    p = np.asarray(power, dtype=float)
    if p.size < 3:
        if p.size == 0 or not np.isfinite(p).any():
            return []
        i = int(np.nanargmax(p))
        return [(float(freq[i]), float(p[i]))]
    interior = (p[1:-1] >= p[:-2]) & (p[1:-1] >= p[2:]) & np.isfinite(p[1:-1])
    idx = np.flatnonzero(interior) + 1
    if idx.size == 0:
        i = int(np.nanargmax(p))
        return [(float(freq[i]), float(p[i]))]
    order = idx[np.argsort(p[idx])[::-1]]
    out, chosen = [], []
    for i in order:
        fi, pi = float(freq[i]), float(p[i])
        if any(abs(fi - fc) < min_sep for fc in chosen):
            continue
        chosen.append(fi)
        out.append((fi, pi))
        if len(out) >= nkeep:
            break
    return out


def find_period(
    t: np.ndarray,
    f: np.ndarray,
    pmin: float = 0.05,
    pmax: float = 60.0,
    hint: float | None = None,
    hint_frac: float = 0.02,
) -> float:
    """Dominant period: fast Lomb-Scargle peak + local exact refine.

    If `hint` (Gaia P) is finite, search only ±hint_frac. Coarse sampling is
    ~8 points per 1/T peak (Press-Rybicki); refine_period later tightens P
    in the time domain, so this only has to land in the right basin."""
    from astropy.timeseries import LombScargle

    baseline = t.max() - t.min()
    pmax = min(pmax, baseline / 3)
    if hint is not None and np.isfinite(hint) and hint > 0:
        lo = max(pmin, hint * (1.0 - hint_frac))
        hi = min(pmax, hint * (1.0 + hint_frac))
        if not (lo < hi):
            lo, hi = hint * (1.0 - hint_frac), hint * (1.0 + hint_frac)
        fmin, fmax = 1.0 / hi, 1.0 / lo
        spp = 12
    else:
        fmin, fmax = 1.0 / pmax, 1.0 / pmin
        spp = 8
    if not (fmax > fmin > 0):
        return float(hint) if (hint is not None and np.isfinite(hint) and hint > 0) else float(pmin)
    ls = LombScargle(t, f)
    freq, power = _ls_spectrum(ls, fmin, fmax, samples_per_peak=spp)
    if not np.any(np.isfinite(power)):
        # A hint longer than the baseline can put the whole window below the
        # frequency resolution, and the periodogram comes back all-NaN. Fall
        # back rather than raising out of the middle of a survey run.
        return float(hint) if (hint is not None and np.isfinite(hint) and hint > 0) else float(pmin)
    i = int(np.nanargmax(power))
    return 1.0 / _refine_freq(ls, _parabolic_argmax(freq, power, i))


def fourier_design(phase: np.ndarray, nharm: int) -> np.ndarray:
    cols = [np.ones_like(phase)]
    for k in range(1, nharm + 1):
        cols += [np.cos(2 * np.pi * k * phase), np.sin(2 * np.pi * k * phase)]
    return np.vstack(cols).T


def camera_design(cam: np.ndarray | None) -> tuple[np.ndarray | None, list[str]]:
    if cam is None:
        return None, []
    names = sorted(set(cam))
    cols = [(cam == c).astype(float) for c in names[1:]]
    return (np.vstack(cols).T if cols else None), names[1:]


def fit_fourier(
    t: np.ndarray,
    f: np.ndarray,
    period: float,
    nharm: int,
    cam: np.ndarray | None = None,
    w: np.ndarray | None = None,
    t0: float | None = None,
) -> tuple[np.ndarray, dict[str, float], float]:
    if t0 is None:
        t0 = t.min()
    A = fourier_design(((t - t0) / period) % 1.0, nharm)
    C, camnames = camera_design(cam)
    if C is not None:
        A = np.hstack([A, C])
    sw = np.sqrt(np.ones_like(f) if w is None else w)
    coef, *_ = np.linalg.lstsq(A * sw[:, None], f * sw, rcond=None)
    nf = 2 * nharm + 1
    return coef[:nf], dict(zip(camnames, coef[nf:])), t0


def cam_offset_array(cam: np.ndarray | None, offsets: dict[str, float]) -> np.ndarray | float:
    return 0.0 if cam is None else np.array([offsets.get(c, 0.0) for c in cam])


def fourier_model(
    t: np.ndarray, coef: np.ndarray, period: float, t0: float, nharm: int
) -> np.ndarray:
    t = np.ascontiguousarray(t, dtype=np.float64)
    coef = np.ascontiguousarray(coef, dtype=np.float64)
    return _fourier_at(t, coef, float(period), float(t0), int(nharm), 0.0)


def fourier_slope(
    t: np.ndarray, coef: np.ndarray, period: float, t0: float, nharm: int
) -> np.ndarray:
    t = np.ascontiguousarray(t, dtype=np.float64)
    coef = np.ascontiguousarray(coef, dtype=np.float64)
    return _fourier_slope_nb(t, coef, float(period), float(t0), int(nharm))


def choose_nharm(
    t: np.ndarray, f: np.ndarray, period: float, cam: np.ndarray | None = None, max_h: int = 12
) -> int:
    best, best_bic, n = 2, np.inf, len(f)
    for h in range(2, max_h + 1):
        coef, off, t0 = fit_fourier(t, f, period, h, cam)
        r = f - cam_offset_array(cam, off) - fourier_model(t, coef, period, t0, h)
        bic = n * np.log(np.sum(r**2) / n) + (2 * h + 1 + len(off)) * np.log(n)
        if bic < best_bic:
            best, best_bic = h, bic
    return best


def robust_template(
    t: np.ndarray,
    f: np.ndarray,
    e: np.ndarray,
    period: float,
    cam: np.ndarray | None = None,
    nharm: int | None = None,
    t0: float | None = None,
    max_nharm: int = 12,
) -> tuple[np.ndarray, dict[str, float], float, int, np.ndarray]:
    """Iteratively fit template + camera offsets with 5-sigma clipping."""
    nharm = nharm or choose_nharm(t, f, period, cam, max_h=max_nharm)
    keep = np.ones(len(t), bool)
    for _ in range(4):
        coef, off, t0 = fit_fourier(
            t[keep],
            f[keep],
            period,
            nharm,
            None if cam is None else cam[keep],
            1 / e[keep] ** 2,
            t0=t0 if t0 is not None else t.min(),
        )
        r = f - cam_offset_array(cam, off) - fourier_model(t, coef, period, t0, nharm)
        scale = stats.median_abs_deviation(r[keep], scale="normal")
        new = np.abs(r) < 5 * scale
        # A perfectly modelled or quantised series leaves a residual whose
        # spread is numerical dust, and 5x dust would discard the whole star.
        # A real 5-sigma clip never removes a meaningful share of the points,
        # so treat a mask that does as a degenerate scale and stop clipping.
        if not np.isfinite(scale) or scale <= 0 or new.sum() * 2 < keep.sum():
            break
        if new.sum() == keep.sum():
            break
        keep = new
    return coef, off, t0, nharm, keep


def template_peak_to_peak(coef: np.ndarray, nharm: int) -> float:
    pp = np.linspace(0, 1, 2000)
    m = fourier_design(pp, nharm) @ coef
    return float(m.max() - m.min())


# A camera whose scale is measured from too few points, or whose points cover
# too little of the template's range, gets 1.0 instead of a noisy number.
CAM_SCALE_MIN_PTS = 20
CAM_SCALE_MIN_RANGE = 0.3  # of the template peak-to-peak
CAM_SCALE_LIMITS = (0.2, 5.0)
CAM_SCALE_TOL = 0.02  # below this the scales are not worth a re-fit


def camera_scales(
    t: np.ndarray,
    f: np.ndarray,
    e: np.ndarray,
    period: float,
    coef: np.ndarray,
    off: dict[str, float],
    t0: float,
    nharm: int,
    cam: np.ndarray | None,
) -> dict[str, float]:
    """Per-camera amplitude scale against the shared template.

    Camera offsets put every camera on the same mean level. They do not put it
    on the same *amplitude* scale, and three things make that a real difference
    rather than a refinement: a pulsator's amplitude in ZTF r is ~0.7-0.8x its
    amplitude in g; two ZTF object ids of one star sit in different
    fields/quadrants with different blending, so each is diluted differently;
    and ASAS-SN cameras differ in aperture and focus. Fit one template across a
    mix of those and the residual becomes a function of phase -- which
    `jitter_test` reads as cycle-to-cycle phase jitter, and `sparse_timings`
    reads as a timing offset, because it divides by a template slope that is
    wrong by the same factor. Merging bands without this would manufacture
    exactly the signal the search is looking for.

    Solves f = b_c + a_c * (T - T0) per camera by weighted least squares, where
    T0 is the template mean, and returns {camera: a_c} normalised so the
    point-weighted mean scale is 1 (the template keeps its own amplitude).
    """
    if cam is None:
        return {}
    names = sorted(set(cam))
    if len(names) < 2:
        return {}
    model = fourier_model(t, coef, period, t0, nharm)
    t0_level = float(coef[0])
    span = template_peak_to_peak(coef, nharm)
    out, npts = {}, {}
    for c in names:
        m = cam == c
        n = int(m.sum())
        x = model[m] - t0_level
        if n < CAM_SCALE_MIN_PTS or span <= 0 or (x.max() - x.min()) < CAM_SCALE_MIN_RANGE * span:
            out[c], npts[c] = 1.0, n
            continue
        y = f[m] - off.get(c, 0.0)
        w = 1.0 / np.maximum(e[m], 1e-30) ** 2
        # [x, 1] design: slope is the amplitude scale, intercept re-absorbs the
        # level so a wrong offset does not leak into the scale.
        A = np.vstack([x, np.ones(n)]).T
        sw = np.sqrt(w)
        try:
            sol, *_ = np.linalg.lstsq(A * sw[:, None], y * sw, rcond=None)
        except np.linalg.LinAlgError:
            out[c], npts[c] = 1.0, n
            continue
        a = float(sol[0])
        lo, hi = CAM_SCALE_LIMITS
        out[c], npts[c] = (a if np.isfinite(a) and lo <= a <= hi else 1.0), n
    total = sum(npts.values())
    mean = sum(out[c] * npts[c] for c in names) / total if total else 1.0
    if not np.isfinite(mean) or mean <= 0:
        return {}
    return {c: out[c] / mean for c in names}


def apply_camera_scales(
    f: np.ndarray,
    cam: np.ndarray | None,
    off: dict[str, float],
    scales: dict[str, float],
    level: float,
) -> np.ndarray:
    """Put every camera on the template's amplitude scale.

    f -> level + (f - b_c - level)/a_c: the excursion about the template mean
    is what gets rescaled, not the flux itself. The offset is removed here and
    re-fit afterwards, because dividing by a_c changes it.
    """
    if cam is None or not scales:
        return f
    a = np.array([scales.get(c, 1.0) for c in cam])
    b = np.array([off.get(c, 0.0) for c in cam])
    return level + (f - b - level) / np.where(a == 0, 1.0, a)


def prewhiten_modes(
    t: np.ndarray,
    r: np.ndarray,
    e: np.ndarray,
    period: float,
    nmax: int = 5,
    fap: float = 1e-3,
    fmax: float | None = None,
    harm_tol: float = 0.02,
    min_amplitude: float = 0.0,
    sideband_tol: float = 0.15,
    min_snr: float = 0.0,
) -> tuple[np.ndarray, list[dict]]:
    """Iteratively remove independent pulsation modes from a residual series.

    Multiperiodic pulsators (beta Cep, RRd, double-mode Cepheids) leave their
    secondary modes in the residual after the dominant-mode template is
    subtracted. Fitting a time-shifted template cycle by cycle then reads the
    beating as O-C variation with the beat period. Each residual spectrum is a
    fast Lomb-Scargle; harmonic leftover of the template is skipped by walking
    down other significant peaks in the *same* spectrum rather than notching a
    million-point grid and recomputing. A peak is fitted as a sinusoid and
    subtracted, up to nmax times or until nothing exceeds the false-alarm
    level. The floor is set by the residual spectrum's own noise: a peak is
    removed only when its fitted amplitude exceeds min_snr times the formal
    amplitude error, scaled by the residual's reduced chi2 (the usual S/N >= 4
    criterion of frequency analysis). A fixed fraction of the dominant
    amplitude (min_amplitude) can be given instead; the 2 per cent floor used
    before review 2 sat above the millimagnitude 0.61 and 0.68 modes of
    first-overtone RR Lyrae stars, which then contaminated the top of the
    ranking (review 2, m8). Peaks within sideband_tol of a harmonic are removed
    but tagged 'sideband_suspect' so stage 2 can see whether the prewhitening
    may have eaten part of a signal (a beacon's own phase-modulation sidebands
    sit there).

    Returns (model, modes): model is the summed sinusoids evaluated on t."""
    from astropy.timeseries import LombScargle

    baseline = t.max() - t.min()
    if len(t) < 100 or baseline <= 0 or nmax <= 0:
        return np.zeros_like(r), []
    dt = np.median(np.diff(np.sort(t)))
    nyq = 0.5 / max(dt, 1e-6)
    f0 = 1.0 / period
    # Independent modes sit near the class period window, not at the photometric
    # Nyquist. Default: up to 6 f0, but at least the bcep short-period edge.
    if fmax is None:
        fmax = min(nyq, max(6.0 * f0, 1.0 / STAR_CLASSES["bcep"].pmin))
    else:
        fmax = min(float(fmax), nyq)
    fmin = 2.0 / baseline
    if not (fmax > fmin):
        return np.zeros_like(r), []
    model = np.zeros_like(r)
    modes = []
    resid = r.copy()
    min_sep = 1.5 / baseline
    w = 1.0 / e
    for _ in range(nmax):
        ls = LombScargle(t, resid, e)
        freq, power = _ls_spectrum(ls, fmin, fmax)
        thresh = ls.false_alarm_level(fap, minimum_frequency=fmin, maximum_frequency=fmax)
        picked = False
        for f_rough, pwr in _peaks_desc(freq, power, min_sep, nkeep=nmax + 16):
            if pwr < thresh:
                break
            fpk = _refine_freq(ls, f_rough)
            ratio = fpk / f0
            k = round(ratio)
            is_harm = k >= 1 and abs(ratio - k) < harm_tol
            is_subharm = any(abs(ratio - 1.0 / m) < harm_tol / m for m in (2, 3, 4))
            if is_harm or is_subharm:
                continue
            if any(abs(fpk / m["frequency_per_day"] - 1.0) < harm_tol for m in modes):
                continue
            A = np.vstack([np.sin(2 * np.pi * fpk * t), np.cos(2 * np.pi * fpk * t)]).T
            c, *_ = np.linalg.lstsq(A * w[:, None], resid * w, rcond=None)
            amp = float(math.hypot(*c))
            if amp < min_amplitude:
                return model, modes
            if min_snr > 0:
                # Formal amplitude error from the normal matrix, inflated by
                # the residual's reduced chi2 so that unmodelled variability
                # counts as noise; the S/N of frequency analysis.
                Aw = A * w[:, None]
                cov = np.linalg.pinv(Aw.T @ Aw)
                chi2_nu = float(np.sum(((resid - A @ c) * w) ** 2)) / max(len(resid) - 2, 1)
                sig_amp = math.sqrt(max(0.5 * float(np.trace(cov)) * chi2_nu, 1e-300))
                if amp / sig_amp < min_snr:
                    return model, modes
            comp = A @ c
            modes.append(
                {
                    "frequency_per_day": fpk,
                    "period_days": 1.0 / fpk,
                    "amplitude": amp,
                    "ratio_to_dominant": ratio,
                    "power": float(pwr),
                    "sideband_suspect": bool(k >= 1 and abs(ratio - k) < sideband_tol),
                    "amplitude_snr": (float(amp / sig_amp) if min_snr > 0 else None),
                }
            )
            resid = resid - comp
            model = model + comp
            picked = True
            break
        if not picked:
            break
    return model, modes


# ============================================================================
# 3. Timing residuals
# ============================================================================


@dataclass
class Timing:
    index: int
    t_center: float
    shift_days: float
    shift_err_days: float
    npts: int
    chi2_red: float
    amp_ratio: float = 1.0  # per-cycle amplitude scale relative to template (dense mode)


def dense_timings(
    t: np.ndarray,
    f: np.ndarray,
    e: np.ndarray,
    period: float,
    coef: np.ndarray,
    t0: float,
    nharm: int,
    min_coverage: float = 0.6,
    boundary_t0: float | None = None,
) -> list[Timing]:
    """Per-cycle template time-shift fit. boundary_t0 lets stage 2 move the cycle cuts.

    The shift is a coarse-then-fine chi2 grid (not Brent). A 3-point Brent
    bracket around the coarse minimum often is not a local min when the
    template is a poor match, and used to abort the whole star.
    """
    if len(t) == 0:
        return []
    t = np.ascontiguousarray(t, dtype=np.float64)
    f = np.ascontiguousarray(f, dtype=np.float64)
    e = np.ascontiguousarray(e, dtype=np.float64)
    coef = np.ascontiguousarray(coef, dtype=np.float64)
    period, t0, nharm = float(period), float(t0), int(nharm)
    bt0 = t0 if boundary_t0 is None else float(boundary_t0)
    invvar = 1.0 / np.maximum(e, 1e-30) ** 2
    cyc = np.floor((t - bt0) / period).astype(np.int64)
    nout, cs, tcs, dts, errs, npts, chi2s, amps = _dense_timings_nb(
        t, f, invvar, cyc, period, coef, t0, nharm, DENSE_MIN_PTS_PER_FIT, float(min_coverage)
    )
    return [
        Timing(
            int(cs[i]),
            float(tcs[i]),
            float(dts[i]),
            float(errs[i]),
            int(npts[i]),
            float(chi2s[i]),
            float(amps[i]),
        )
        for i in range(int(nout))
    ]


def epoch_offsets(
    t: np.ndarray,
    f: np.ndarray,
    e: np.ndarray,
    period: float,
    coef: np.ndarray,
    t0: float,
    nharm: int,
    slope_frac: float = 0.3,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Each epoch's implied timing offset, its error, its cycle, and whether
    the template was steep enough there to read one.

    The residual from the template divided by the template slope at that
    phase is the time shift that would have produced it; epochs near an
    extremum (slope below `slope_frac` of the maximum) carry no timing and
    are masked. This is the series both the block O-C (`sparse_timings`) and
    the per-cycle series (`cycle_timings`) are built from.
    """
    resid = f - fourier_model(t, coef, period, t0, nharm)
    slope = fourier_slope(t, coef, period, t0, nharm)
    smax = np.max(
        np.abs(fourier_slope(np.linspace(t0, t0 + period, 2000), coef, period, t0, nharm))
    )
    use = np.abs(slope) > slope_frac * smax
    safe_slope = np.where(use, slope, 1.0)
    dt = np.where(use, -resid / safe_slope, np.nan)
    dterr = np.where(use, e / np.abs(safe_slope), np.nan)
    cyc = np.floor((t - t0) / period).astype(int)
    return dt, dterr, cyc, use


def cycle_timings(
    t: np.ndarray, dt: np.ndarray, dterr: np.ndarray, cyc: np.ndarray, use: np.ndarray
) -> list[Timing]:
    """One timing per *observed cycle*, from the epochs that fall in it.

    Ground-based surveys sample a short-period star once or twice a night,
    so most cycles are never seen and the ones that are hold one or two
    epochs. Block averaging (`sparse_timings`) recovers a smooth O-C from
    that at the price of the bits: a block is K ~ 100 cycles wide and the
    step between blocks is the sum of ~100 bits, Gaussian by the central
    limit theorem. This keeps the cycle: the inverse-variance mean of the
    epochs in each observed cycle, with its error inflated by sqrt(chi2_red)
    where there is more than one, indexed by cycle number. Folding these on
    the word (`analyse_epoch_fold`) integrates over every repeat the baseline
    holds without ever needing two adjacent cycles to have been observed.
    """
    m = use & np.isfinite(dt) & np.isfinite(dterr) & (dterr > 0)
    if not m.any():
        return []
    out = []
    order = np.argsort(cyc[m], kind="stable")
    cs = cyc[m][order]
    ts = t[m][order]
    ds = dt[m][order]
    es = dterr[m][order]
    starts = np.flatnonzero(np.r_[True, cs[1:] != cs[:-1]])
    ends = np.r_[starts[1:], cs.size]
    for a, b in zip(starts, ends):
        w = 1.0 / es[a:b] ** 2
        mean = float(np.sum(w * ds[a:b]) / np.sum(w))
        npts = int(b - a)
        chi2_red = float(np.sum(w * (ds[a:b] - mean) ** 2) / (npts - 1)) if npts > 1 else 1.0
        err = math.sqrt(1.0 / float(np.sum(w))) * math.sqrt(max(chi2_red, 1.0))
        out.append(Timing(int(cs[a]), float(ts[a:b].mean()), mean, err, npts, chi2_red))
    return out


def sparse_timings(
    t: np.ndarray,
    f: np.ndarray,
    e: np.ndarray,
    period: float,
    coef: np.ndarray,
    t0: float,
    nharm: int,
    slope_frac: float = 0.3,
    target_block_pts: int = SPARSE_MIN_PTS_PER_BLOCK,
    block_offset: int = 0,
    blocks_k: int | None = None,
) -> tuple[list[Timing], int]:
    """Epoch-level implied timing offsets averaged over blocks of K cycles.

    Returns (timings, K). `block_offset` shifts where the block boundaries
    fall, which is the sparse analogue of moving the per-cycle cuts in dense
    mode: real timing structure survives a re-blocking, an artifact of where
    the averaging happened to land does not. `blocks_k` pins K so a shifted
    run re-blocks the same data rather than re-deriving a different width.
    """
    dt, dterr, cyc, use = epoch_offsets(t, f, e, period, coef, t0, nharm, slope_frac)
    pts_per_cycle = len(t) / max(cyc.max() + 1, 1)
    if blocks_k:
        K = max(1, int(blocks_k))
    else:
        K = max(1, math.ceil(target_block_pts / max(pts_per_cycle * use.mean(), 1e-6)))
    blk = (cyc + int(block_offset)) // K
    out = []
    for b in np.unique(blk):
        m = (blk == b) & use
        if m.sum() < max(3, target_block_pts // 2):
            continue
        # Clip within the block against the block's own timing, for the same
        # reason dense_timings clips per cycle: the template mask is taken
        # against the mean ephemeris and would remove a genuinely shifted
        # epoch. Never drop more than half the block.
        med = np.median(dt[m])
        mad = 1.4826 * np.median(np.abs(dt[m] - med))
        if mad > 0:
            inl = m & (np.abs(dt - med) < 5 * mad)
            if inl.sum() >= max(3, target_block_pts // 2) and 2 * inl.sum() >= m.sum():
                m = inl
        w = 1 / dterr[m] ** 2
        mean = np.sum(w * dt[m]) / np.sum(w)
        chi2_red = np.sum(w * (dt[m] - mean) ** 2) / max(m.sum() - 1, 1)
        err = math.sqrt(1 / np.sum(w)) * math.sqrt(max(chi2_red, 1.0))
        out.append(
            Timing(
                int(b), float(t[m].mean()), float(mean), float(err), int(m.sum()), float(chi2_red)
            )
        )
    return out, K


def timing_drift_design(
    t: np.ndarray, slope: np.ndarray, period: float, t0: float, mod_period: float | None = None
) -> np.ndarray:
    """Design matrix of a slow timing drift, linearised into flux.

    A timing offset dt(t) small against the scale over which the template
    slope changes shows up in the flux residual as r = -F'(phi) dt(t). The
    columns are -F' times [1, x, x^2] with x = (t - t0)/P, i.e. a quadratic
    ephemeris (period error and a constant Pdot, the terms the dense path's
    detrending removes), plus [sin, cos] at `mod_period` when one is given
    (a Blazhko-range modulation, coherent over tens to hundreds of days)."""
    x = (np.asarray(t, dtype=float) - float(t0)) / float(period)
    x = x - x.mean()
    cols = [np.ones_like(x), x, x * x]
    if mod_period:
        w = 2.0 * np.pi * np.asarray(t, dtype=float) / float(mod_period)
        cols += [np.sin(w), np.cos(w)]
    return -np.asarray(slope, dtype=float)[:, None] * np.vstack(cols).T


def fit_timing_drift(
    r: np.ndarray, e: np.ndarray, A: np.ndarray
) -> tuple[np.ndarray, np.ndarray, float]:
    """Weighted least squares of the residual on a drift design.

    Returns (coefficients, fitted flux component, chi2 of the fit)."""
    w = 1.0 / np.asarray(e, dtype=float)
    c, *_ = np.linalg.lstsq(A * w[:, None], np.asarray(r, dtype=float) * w, rcond=None)
    comp = A @ c
    chi2 = float(np.sum(((r - comp) * w) ** 2))
    return c, comp, chi2


def jitter_test(
    t: np.ndarray,
    f: np.ndarray,
    e: np.ndarray,
    period: float,
    coef: np.ndarray,
    t0: float,
    nharm: int,
    nboot: int = JITTER_NBOOT,
    rng: np.random.Generator | None = None,
    detrend: bool = True,
) -> dict:
    """Regress resid^2 on template slope^2: slope coefficient = sigma_t^2 (phase jitter).

    With `detrend` (the default) a quadratic timing drift is removed from the
    residuals first, so that a secular period change over a multi-year
    ground-based baseline is not read as jitter: a typical RR Lyrae Pdot of a
    few 1e-10 d/d drifts the phase by tens of seconds over five years,
    comparable to the codes the sparse arm detects (review 2, m5). The dense
    path removes the same quadratic from its O-C.

    p_jitter_bootstrap is add-one (see shuffle_pvalue), so it is floored at
    1/(nboot+1) rather than collapsing to 0; jitter_nboot is reported next to
    it so a consumer can tell saturation from significance.
    """
    r = f - fourier_model(t, coef, period, t0, nharm)
    slope = fourier_slope(t, coef, period, t0, nharm)
    if detrend and len(t) > 12:
        _c, comp, _chi2 = fit_timing_drift(r, e, timing_drift_design(t, slope, period, t0))
        r = r - comp
    r2 = r**2
    s2 = slope**2
    r2 = np.ascontiguousarray(r2, dtype=np.float64)
    s2 = np.ascontiguousarray(s2, dtype=np.float64)
    var_t = float(_slope_coeff(s2, r2))
    seed = 7 if rng is None else int(rng.integers(0, 2**31 - 1))
    boots = _jitter_boots(s2, r2, int(nboot), seed)
    err = float(boots.std())
    return {
        "jitter_sigma_t_seconds": math.sqrt(max(var_t, 0.0)) * 86400,
        "jitter_sigma_t_err_seconds": float(math.sqrt(err) * 86400) if err > 0 else None,
        "jitter_var_t_over_err": float(var_t / err) if err > 0 else None,
        "jitter_nboot": int(nboot),
        "p_jitter_bootstrap": shuffle_pvalue(boots, 0.0, upper=False),
    }


# Sampling periods (days) a coherent-jitter fit may not land on: daily and
# half-daily aliases, the synodic month and its half, the year and its
# harmonics. A modulation within COHERENT_ALIAS_TOL of any of these is the
# survey's window function, not the star.
COHERENT_ALIAS_PERIODS = (0.5, 1.0, 14.77, 29.53, 91.31, 121.75, 182.62, 365.25)
COHERENT_ALIAS_TOL = 0.05
# Modulation periods searched, days: the Blazhko range and long-period
# binaries, capped at half the baseline so the sinusoid is separable from the
# quadratic ephemeris the design already carries.
COHERENT_PMIN = 10.0
COHERENT_PMAX = 1000.0
COHERENT_NPERIODS = 300
# Minimum improvement in chi2 for the sinusoid to count: 2 free parameters, a
# look-elsewhere over COHERENT_NPERIODS trial periods (2 ln 300 ~ 11), and a
# margin. A sinusoid that clears this on a ~1000-epoch light curve is real.
COHERENT_MIN_DCHI2 = 30.0
# The fraction of the (quadratic-detrended) jitter variance the coherent
# sinusoid has to explain for the jitter to be called modulation.
COHERENT_MIN_EXPLAINED = 0.5
# A repeated word is itself coherent, at the word period: L bits at one bit
# per cycle repeat every L pulsation cycles, and a sinusoid at that period
# explained 75% of the synthetic primes control's jitter. So a fitted period
# that is a whole number of cycles in the decoder's word-length range is
# folded there and its steps are asked to be squarer than a sinusoid, as the
# dense squareness gate asks (Section 4.7 of the paper); a code-like fold is
# not vetoed.
COHERENT_WORD_MIN, COHERENT_WORD_MAX = 4, 64
COHERENT_BIT_CYCLES = (1, 2, 4, 8, 16, 32)
COHERENT_MIN_EPOCHS_PER_RESIDUE = 3
# The code model has to beat the wave by this many nats (odds of ~150:1)
# before the fold is called code-like: on white noise and on a sinusoid that
# happens to sit at a whole number of cycles the best split wins by 0-2 nats,
# an injected slow word by 75-200.
COHERENT_CODE_MIN_GAIN_NATS = 5.0


def _fold_code_evidence(t: np.ndarray, dt: np.ndarray, w: np.ndarray, pm: float, n: int) -> dict:
    """Are per-epoch offsets folded on the modulation period two levels or a wave?

    The offsets are binned into n phase bins of the fitted period `pm` (one
    bin per bit when pm is the word period of a one-bit-per-cycle code),
    averaged with their errors, and differenced into steps. Two models of the
    steps are priced in nats, BIC-style: the best sinusoid at any frequency
    (three parameters and a look-elsewhere over the grid) against the best
    two-level split (two levels plus ln 2 per bit, the fold decoder's price).
    The bit price is what stops a two-level split of noisy steps winning by
    default: an optimal split removes 2/pi of Gaussian noise variance on its
    own, which is worth less than ln 2 per step. A positive `cost_gain_nats`
    means the steps are cheaper as a code than as a wave.
    """
    out = {"n_bins": n, "squareness": None, "cost_gain_nats": None, "code_like": False}
    phase = np.mod((t - t.min()) / float(pm), 1.0)
    res = np.minimum((phase * n).astype(np.int64), n - 1)
    m, s2 = np.zeros(n), np.zeros(n)
    for r in range(n):
        sel = (res == r) & (w > 0)
        if sel.sum() < COHERENT_MIN_EPOCHS_PER_RESIDUE:
            return out
        m[r] = np.average(dt[sel], weights=w[sel])
        s2[r] = float(np.average((dt[sel] - m[r]) ** 2, weights=w[sel])) / sel.sum()
    steps = np.roll(m, -1) - m
    var = np.roll(s2, -1) + s2
    var = np.maximum(var, 1e-3 * var.mean() if var.mean() > 0 else 1e-300)
    wt = 1.0 / var
    x = np.arange(n, dtype=float)
    # the wave: best sinusoid plus a constant at any frequency of the fold
    grid = np.linspace(1.0 / n, 0.5, 4 * n)
    chi_sin = np.inf
    for f in grid:
        A = np.vstack([np.ones(n), np.sin(2 * np.pi * f * x), np.cos(2 * np.pi * f * x)]).T
        c, *_ = np.linalg.lstsq(A * np.sqrt(wt)[:, None], steps * np.sqrt(wt), rcond=None)
        chi_sin = min(chi_sin, float(np.sum(wt * (steps - A @ c) ** 2)))
    # the code: the split of the sorted steps with the smallest chi2
    order = np.argsort(steps)
    chi_word = np.inf
    for i in range(1, n):
        lo, hi = order[:i], order[i:]
        mlo = np.average(steps[lo], weights=wt[lo])
        mhi = np.average(steps[hi], weights=wt[hi])
        chi = float(
            np.sum(wt[lo] * (steps[lo] - mlo) ** 2) + np.sum(wt[hi] * (steps[hi] - mhi) ** 2)
        )
        chi_word = min(chi_word, chi)
    cost_sin = 0.5 * chi_sin + 1.5 * math.log(n) + math.log(grid.size)
    cost_word = 0.5 * chi_word + math.log(n) + n * math.log(2.0)
    out.update(
        {
            "squareness": float(max(0.0, 1.0 - chi_word / chi_sin)) if chi_sin > 0 else 0.0,
            "cost_gain_nats": float(cost_sin - cost_word),
            "code_like": bool(cost_sin - cost_word >= COHERENT_CODE_MIN_GAIN_NATS),
        }
    )
    return out


def coherent_jitter_test(
    t: np.ndarray,
    f: np.ndarray,
    e: np.ndarray,
    period: float,
    coef: np.ndarray,
    t0: float,
    nharm: int,
) -> dict:
    """Is the per-epoch phase jitter coherent over a Blazhko-range period?

    A one-bit-per-cycle code's implied per-epoch timing offsets are white; a
    Blazhko star's are coherent over the modulation period, tens to hundreds
    of days. This is the sparse-mode analogue of the dense path's beat veto
    (review 2, M3): fit one slow sinusoid plus the quadratic ephemeris to the
    linearised offsets, re-measure the jitter variance on what is left, and
    report how much of it the sinusoid took. The caller vetoes when the
    explained fraction clears COHERENT_MIN_EXPLAINED, the fit clears
    COHERENT_MIN_DCHI2 and the period is not a sampling alias.
    """
    t = np.asarray(t, dtype=float)
    r = f - fourier_model(t, coef, period, t0, nharm)
    slope = fourier_slope(t, coef, period, t0, nharm)
    baseline = float(t.max() - t.min()) if t.size else 0.0
    pmax = min(COHERENT_PMAX, baseline / 2.0)
    out = {
        "period_days": None,
        "semi_amplitude_seconds": None,
        "delta_chi2": None,
        "sigma_t_before_seconds": None,
        "sigma_t_after_seconds": None,
        "explained_fraction": None,
        "is_alias": None,
        "coherent": False,
    }
    if t.size < 50 or not (pmax > COHERENT_PMIN):
        return out
    A0 = timing_drift_design(t, slope, period, t0)
    _c0, comp0, chi0 = fit_timing_drift(r, e, A0)
    s2 = np.ascontiguousarray(slope**2, dtype=np.float64)
    var_before = float(_slope_coeff(s2, np.ascontiguousarray((r - comp0) ** 2, dtype=np.float64)))
    best = None
    for pm in np.geomspace(COHERENT_PMIN, pmax, COHERENT_NPERIODS):
        c, comp, chi = fit_timing_drift(r, e, timing_drift_design(t, slope, period, t0, pm))
        if best is None or chi < best[0]:
            best = (chi, float(pm), c, comp)
    chi1, pm, c, comp = best
    # The residual fit's own reduced chi2 scales the improvement, as it does
    # the jitter bootstrap: the errors are the survey's, not the star's.
    dof = max(t.size - A0.shape[1] - 2, 1)
    dchi2 = (chi0 - chi1) / max(chi1 / dof, 1e-300)
    var_after = float(_slope_coeff(s2, np.ascontiguousarray((r - comp) ** 2, dtype=np.float64)))
    # Clamped: the slope estimate can come out negative on white noise.
    explained = min(1.0, max(0.0, 1.0 - var_after / var_before)) if var_before > 0 else 0.0
    alias = any(abs(pm - a) / a < COHERENT_ALIAS_TOL for a in COHERENT_ALIAS_PERIODS)
    # A repeated word is coherent at the word period; fold there, one bin per
    # pulsation cycle, and ask whether the steps are two levels or a wave.
    # The fold is in time on the fitted period, not by cycle index, because a
    # per-mille period error smears a cycle-index fold over a long baseline.
    # One bin per bit at every bit period the ground-based arm can hear (the
    # slow codes of Table 3 hold a bit for 8-32 cycles); the best fold decides.
    ev = {"n_bins": None, "squareness": None, "cost_gain_nats": None, "code_like": False}
    ok = np.abs(slope) > 0
    dt = np.zeros(t.size)
    dt[ok] = -(r - comp0)[ok] / slope[ok]
    w = np.where(ok, slope**2 / np.asarray(e, dtype=float) ** 2, 0.0)
    for k in COHERENT_BIT_CYCLES:
        n_word = round(pm / (k * float(period)))
        if COHERENT_WORD_MIN <= n_word <= COHERENT_WORD_MAX:
            got = _fold_code_evidence(t, dt, w, pm, n_word)
            if got["cost_gain_nats"] is not None and (
                ev["cost_gain_nats"] is None or got["cost_gain_nats"] > ev["cost_gain_nats"]
            ):
                ev = {**got, "bit_cycles": k}
    code_like = bool(ev["code_like"])
    out.update(
        {
            "period_days": pm,
            # The slope is per day, so the drift coefficients are in days.
            "semi_amplitude_seconds": float(math.hypot(c[3], c[4]) * 86400),
            "delta_chi2": float(dchi2),
            "sigma_t_before_seconds": math.sqrt(max(var_before, 0.0)) * 86400,
            "sigma_t_after_seconds": math.sqrt(max(var_after, 0.0)) * 86400,
            "explained_fraction": float(explained),
            "is_alias": bool(alias),
            "word_cycles": ev["n_bins"],
            "word_bit_cycles": ev.get("bit_cycles"),
            "fold_squareness": ev["squareness"],
            "fold_cost_gain_nats": ev["cost_gain_nats"],
            "code_like": code_like,
            "coherent": bool(
                dchi2 >= COHERENT_MIN_DCHI2
                and explained >= COHERENT_MIN_EXPLAINED
                and not alias
                and not code_like
            ),
        }
    )
    return out


# ============================================================================
# 4/5. Statistics on the O-C series
# ============================================================================


def shuffle_pvalue(null: Any, obs: float, upper: bool = True) -> float:
    """Add-one Monte-Carlo p-value: (#{null at least as extreme} + 1) / (N + 1).

    The plain fraction returns exactly 0 as soon as the observed statistic
    beats every draw, and the log-evidence curve in score.py reads p = 0 as
    ~1e-12 -- six orders of magnitude below anything N draws can resolve. The
    add-one form floors the p-value at 1/(N+1), which is the honest limit, and
    also keeps a non-finite statistic from scoring as p = 0.
    """
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    n = null.size
    if n == 0 or obs is None or not np.isfinite(obs):
        return 1.0
    hits = int(np.count_nonzero(null >= obs if upper else null <= obs))
    return float((hits + 1) / (n + 1))


def primes_upto(n: int) -> list[int]:
    s = np.ones(n + 1, bool)
    s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]:
            s[i * i :: i] = False
    return np.flatnonzero(s)


def acf(x: Any, maxlag: int) -> np.ndarray:
    x = np.ascontiguousarray(x, dtype=np.float64)
    return _acf_nb(x, int(maxlag))


def lb_maxlag(n: int) -> int:
    """Lag budget for the Ljung-Box / ACF tests.

    Deliberately a large fraction of n: a beacon repeating a multi-byte word
    puts its structure at the word length, so capping at the textbook ~10 lags
    would look straight past the signal. The price is that the chi2
    approximation for Q is not trustworthy this far out, which is why the
    reported p-value is calibrated by permutation instead (see ljung_box_p).
    """
    return max(5, min(60, n // 3))


def lb_q_from_acf(r: Any, n: int, maxlag: int) -> float:
    """Ljung-Box Q from an already-computed ACF (r[0] is lag 1)."""
    return float(n * (n + 2) * np.sum(np.asarray(r) ** 2 / (n - np.arange(1, maxlag + 1))))


_LN10 = math.log(10.0)


def log10_sf(dist: Any, stat: Any, *args: Any) -> float:
    """log10 of an upper tail probability.

    scipy's .sf() underflows to exactly 0.0 for an extreme statistic, and a
    stored p-value of 0 is indistinguishable from one that genuinely sat at the
    double-precision floor. .logsf() keeps the magnitude, so a tail of 1e-320
    is reported as -320 instead of thrown away.
    """
    try:
        val = float(dist.logsf(stat, *args)) / _LN10
    except (ValueError, TypeError, ZeroDivisionError):
        return float("nan")
    return val if math.isfinite(val) else float("nan")


def _chi2_log10_sf_asymptotic(q: float, df: int) -> float:
    """Leading term of the upper incomplete gamma, valid for Q >> df.

    log Q(a, x) ~ -x + (a - 1) ln x - lgamma(a), with a = df/2, x = Q/2.
    Agrees with chi2.logsf to ~0.02 dex where that still works, and carries on
    past ~1e-273 where it underflows.
    """
    a, x = df / 2.0, q / 2.0
    if x <= 0.0 or a <= 0.0:
        return float("nan")
    return (-x + (a - 1.0) * math.log(x) - math.lgamma(a)) / _LN10


def chi2_log10_sf(q: Any, df: Any) -> float:
    """log10 of the chi2 upper tail, without collapsing to an exact zero."""
    try:
        q, df = float(q), int(df)
    except (TypeError, ValueError):
        return float("nan")
    if not (math.isfinite(q) and q > 0.0 and df > 0):
        return float("nan")
    val = log10_sf(stats.chi2, q, df)
    return val if math.isfinite(val) else _chi2_log10_sf_asymptotic(q, df)


def _f_log10_sf_asymptotic(stat: float, d1: int, d2: int) -> float:
    """Leading term of the F upper tail, valid for stat >> 1.

    sf(x; d1, d2) = I_z(d2/2, d1/2) with z = d2 / (d2 + d1 x), and for small z
    the regularised incomplete beta is I_z(a, b) ~ z**a / (a B(a, b)). scipy's
    own .logsf() underflows to -inf above stat ~ 1e6 -- far sooner than the
    chi2 one -- so a strong F-test would otherwise lose its magnitude entirely.
    """
    a, b = d2 / 2.0, d1 / 2.0
    if not (stat > 0.0 and a > 0.0 and b > 0.0):
        return float("nan")
    z = d2 / (d2 + d1 * stat)
    if z <= 0.0:
        return float("-inf")
    log_beta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    return (a * math.log(z) - math.log(a) - log_beta) / _LN10


def f_log10_sf(stat: Any, d1: Any, d2: Any) -> float:
    """log10 of the F upper tail, without collapsing to an exact zero."""
    try:
        stat, d1, d2 = float(stat), int(d1), int(d2)
    except (TypeError, ValueError):
        return float("nan")
    if not (math.isfinite(stat) and stat > 0.0 and d1 > 0 and d2 > 0):
        return float("nan")
    val = log10_sf(stats.f, stat, d1, d2)
    return val if math.isfinite(val) else _f_log10_sf_asymptotic(stat, d1, d2)


def ljung_box(x: Any, maxlag: int) -> tuple[float, float]:
    """(Q, chi2 p). The chi2 tail is unreliable at maxlag ~ n/3; kept for
    reference and for callers with no permutation null to hand."""
    n = len(x)
    r = acf(x, maxlag)[1:]
    q = lb_q_from_acf(r, n, maxlag)
    return q, float(stats.chi2.sf(q, maxlag))


def ljung_box_p(
    x: Any, maxlag: int | None = None, nperm: int = 1000, rng: np.random.Generator | None = None
) -> dict:
    """Ljung-Box Q with a permutation p-value.

    Under the null the series is exchangeable, so permuting it and recomputing
    Q gives an exact reference distribution -- no chi2 assumption, which is
    what breaks down at the lag budget this pipeline needs (see lb_maxlag).
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    maxlag = lb_maxlag(n) if maxlag is None else int(maxlag)
    q, p_chi2 = ljung_box(x, maxlag)
    rng = rng or np.random.default_rng(17)
    null = np.empty(int(nperm))
    for i in range(int(nperm)):
        null[i] = lb_q_from_acf(acf(rng.permutation(x), maxlag)[1:], n, maxlag)
    return {
        "ljung_box_Q": q,
        "p_ljung_box": shuffle_pvalue(null, q),
        "p_ljung_box_chi2": p_chi2,
        "log10_p_ljung_box_chi2": chi2_log10_sf(q, maxlag),
        "ljung_box_maxlag": maxlag,
        "ljung_box_draws": int(nperm),
    }


def permutation_entropy(x: Any, order: int = 4, delay: int = 1) -> float:
    x = np.ascontiguousarray(x, dtype=np.float64)
    if order == 4:
        pe = _perm_entropy4(x, int(delay))
        return None if not np.isfinite(pe) else float(pe)
    n = len(x) - (order - 1) * delay
    if n < 10:
        return None
    pats = {p: i for i, p in enumerate(permutations(range(order)))}
    counts = np.zeros(len(pats))
    for i in range(n):
        counts[pats[tuple(np.argsort(x[i : i + order * delay : delay]))]] += 1
    p = counts[counts > 0] / counts.sum()
    return float(-np.sum(p * np.log2(p)) / np.log2(len(pats)))


_PE_NULL_CACHE: dict[tuple[int, int], np.ndarray] = {}
_PE_NULL_CACHE_MAX = 1024


def permutation_entropy_null(n: int, ndraw: int, order: int = 4, delay: int = 1) -> np.ndarray:
    """Null distribution of the permutation entropy for a series of length n.

    Permutation entropy reads only the *ordering* inside each window, so its
    null depends on n and nothing else -- not on the values, and not on which
    star they came from. Tabulating it once per (n, ndraw) takes one of the
    three statistics out of the per-star shuffle loop, which otherwise redrew
    the same distribution nshuffle times for every one of tens of thousands of
    stars.

    Ranks are drawn from arange(n), so the draws are ties-free; an O-C series
    of floats effectively never ties, and _perm_entropy4 breaks any tie it does
    meet by index.
    """
    key = (int(n), int(ndraw))
    hit = _PE_NULL_CACHE.get(key)
    if hit is not None:
        return hit
    # Deterministic per n: the same table whichever star asks for it first.
    rng = np.random.default_rng(0x9E3779B9 ^ (key[0] * 1_000_003) ^ key[1])
    ranks = np.arange(key[0], dtype=np.float64)
    draws = []
    for _ in range(key[1]):
        pe = permutation_entropy(rng.permutation(ranks), order=order, delay=delay)
        if pe is not None:
            draws.append(pe)
    table = np.asarray(draws, dtype=float)
    if len(_PE_NULL_CACHE) >= _PE_NULL_CACHE_MAX:
        _PE_NULL_CACHE.clear()
    _PE_NULL_CACHE[key] = table
    return table


def is_alias(
    period_days: float, pulsation_period: float, ground_based: bool, tol: float = 0.02
) -> bool:
    if not ground_based:
        return False
    fsec, fpul = 1 / period_days, 1 / pulsation_period
    suspects = [1.0, 2.0, 1 / 29.53, 1 / 365.25]
    suspects += [abs(fpul - k) for k in (1, 2)] + [fpul + 1, fpul + 2]
    suspects += [abs(fpul - 1 / 365.25), fpul + 1 / 365.25]
    return any(abs(fsec - s) / s < tol for s in suspects if s > 0)


def secondary_frequencies(
    tc: np.ndarray,
    oc: np.ndarray,
    sample_period: float,
    pulsation_period: float,
    ground_based: bool,
    nmax: int = 3,
    fap: float = 1e-3,
) -> list[dict]:
    from astropy.timeseries import LombScargle

    baseline = tc.max() - tc.min()
    if baseline < 5 * sample_period or len(tc) < 12:
        return []
    fmin, fmax = 1 / baseline, 0.5 / sample_period
    freq = np.linspace(fmin, fmax, 20_000)
    ls = LombScargle(tc, oc)
    thresh = ls.false_alarm_level(fap, minimum_frequency=fmin, maximum_frequency=fmax)
    found, resid = [], oc.copy()
    for _ in range(nmax):
        ls = LombScargle(tc, resid)
        p = ls.power(freq)
        i = int(np.argmax(p))
        if p[i] < thresh:
            break
        P = float(1 / freq[i])
        found.append(
            {
                "period_days": P,
                "power": float(p[i]),
                "period_ratio_to_pulsation": P / pulsation_period,
                "alias_suspect": bool(is_alias(P, pulsation_period, ground_based)),
            }
        )
        resid = resid - ls.model(tc, freq[i]) + resid.mean()
    return found


def _two_level_mid_np(x: np.ndarray) -> float | None:
    """`sequence.two_level_mid` on an array: midway between the 20th and 80th
    percentiles, None without a spread. numpy's linear percentile is the same
    interpolation `sequence._percentile` uses, so the two agree to rounding."""
    if x.size == 0:
        return None
    lo, hi = np.percentile(x, [20.0, 80.0])
    if not hi > lo:
        return None
    return float((lo + hi) / 2.0)


def code_increments(oc: np.ndarray, cyc: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The first difference of an O-C series in cycle order, per cycle.

    Returns the increments, the cycle each one belongs to (the later cycle of
    the pair, where the step lands), and the positions of the pairs that made
    them -- which is how the shuffle null differences each permuted draw the
    same way (`_null_loop_nb`). Mirrors `sequence.increments`: only
    consecutive cycles make a step; a step across a gap (unseen bits plus
    whatever a sector boundary did to the level) is dropped, and duplicates
    contribute nothing. The steps are winsorized the way
    `sequence.code_points` winsorizes them (`sequence.clip_steps`), so one
    outlying step cannot buy a word on its own; the null loop clips every
    draw the same way (`_clip_steps_nb`), so the null is drawn from the same
    series the observation is read from.
    """
    cyc = np.asarray(cyc, dtype=np.int64)
    pairs = np.flatnonzero(np.diff(cyc) == 1)
    oc = np.asarray(oc, dtype=float)
    return sequence.clip_steps(oc[pairs + 1] - oc[pairs]), cyc[pairs + 1], pairs


# How many shuffle draws are handed to the null kernel at once. The block
# only bounds the memory the index matrix takes (draws x n_series int64); the
# draws themselves come off the generator in the same order whatever it is.
NULL_BLOCK_DRAWS = 4096
NULL_BLOCK_CELLS = 1_000_000

# ---- kicks -----------------------------------------------------------------
# A step this many robust sigmas from the median of the steps is a kick. Four,
# not the clip's three: the Gaussian tail past 4 sigma is 6e-5 per step, so a
# thousand-cycle series produces one by chance about once in sixteen, and the
# route below asks for KICK_MIN of them on one side before it reads anything.
KICK_SIGMA = 4.0
KICK_MIN = 4
# Kicks on the other side of the median from the majority: a bad cycle is one
# step down and one back up, so a code with a glitch in it has one. More than
# a quarter of the kicks pointing the other way is not one-sided.
KICK_COUNTER_MAX_FRACTION = 0.25

# The integer catalogs the kick gaps are read against, flattened for the
# compiled kernel: sequence.KICK_GAP_SEQUENCES in its own order, so the
# kernel's (sequence index, start) names the same window the Python matcher
# labels.
_KICK_SEQ_FLAT = np.concatenate(
    [np.asarray(seq, dtype=np.int64) for _n, _e, seq in sequence.KICK_GAP_SEQUENCES]
)
_KICK_SEQ_OFF = np.cumsum(
    [0] + [len(seq) for _n, _e, seq in sequence.KICK_GAP_SEQUENCES], dtype=np.int64
)


def _is_prime_table(n: int) -> np.ndarray:
    table = np.zeros(max(int(n), 2) + 1, dtype=np.bool_)
    table[primes_upto(table.size - 1)] = True
    return table


def find_kicks(inc: np.ndarray, inc_cyc: np.ndarray) -> dict:
    """The rare, large, one-sided steps of an O-C: the set bits of a code
    whose transmitter saves up and kicks.

    `inc` are the raw (unclipped) steps in cycle order and `inc_cyc` the
    cycle each lands on. A kick is a step more than KICK_SIGMA robust sigmas
    from the median; the majority side is the code's side and the minority
    are counted against it. Returns the kick cycles, their positions in the
    step series, the sign, the counter-kick count and whether the set is
    one-sided enough to read.
    """
    y = np.asarray(inc, dtype=float)
    cyc = np.asarray(inc_cyc, dtype=np.int64)
    out = {
        "kick_n": 0,
        "kick_counter_n": 0,
        "kick_sign": 0,
        "kick_cycles": [],
        "kick_positions": np.zeros(0, dtype=np.int64),
        "kick_one_sided": False,
        "kick_sigma_step": None,
    }
    if y.size < 8:
        return out
    q1, med, q3 = np.percentile(y, [25, 50, 75])
    scale = (q3 - q1) / 1.349
    if not scale > 0:
        return out
    dev = y - med
    lim = KICK_SIGMA * scale
    up = np.flatnonzero(dev > lim)
    down = np.flatnonzero(dev < -lim)
    if up.size == down.size == 0:
        return out
    if up.size == down.size:
        return {**out, "kick_counter_n": int(up.size)}
    sign = 1 if up.size > down.size else -1
    pos = up if sign > 0 else down
    counter = down.size if sign > 0 else up.size
    out.update(
        {
            "kick_n": int(pos.size),
            "kick_counter_n": int(counter),
            "kick_sign": sign,
            "kick_cycles": [int(c) for c in cyc[pos]],
            "kick_positions": pos,
            "kick_one_sided": bool(counter <= KICK_COUNTER_MAX_FRACTION * pos.size),
            "kick_sigma_step": float(np.median(np.abs(dev[pos])) / scale),
        }
    )
    return out


def kick_word(
    inc: np.ndarray, inc_cyc: np.ndarray, nshuffle: int, rng: np.random.Generator
) -> dict:
    """Read the kicks of an O-C as a message written in their spacing.

    Two statistics on the set of kick cycles, each calibrated against
    `nshuffle` random sets of as many kicks drawn from the same observed
    cycles: the longest run of consecutive gaps that reads as a window of an
    integer catalog in either encoding (`sequence.kick_gap_match`), and the
    fraction of kick pairs whose separation is a prime number -- the
    statistic the one previous search of this kind used, asked of the kicks
    rather than of the autocorrelation. Both are add-one Monte-Carlo
    p-values, and the null is a null over *which* cycles kicked, so the
    number of kicks, the gaps in the data and the catalog's own trials are
    all inside it.
    """
    y = np.asarray(inc, dtype=float)
    cyc = np.asarray(inc_cyc, dtype=np.int64)
    found = find_kicks(y, cyc)
    out = {k: v for k, v in found.items() if k != "kick_positions"}
    out.update(
        {
            "kick_gaps": None,
            "kick_label": None,
            "kick_run": 0,
            "p_kick_catalog": None,
            "kick_prime_frac": None,
            "p_kick_prime": None,
            "kick_draws": 0,
        }
    )
    m = found["kick_n"]
    if m < KICK_MIN or not found["kick_one_sided"]:
        return out
    order = np.argsort(cyc, kind="stable")
    cyc_obs = cyc[order]
    pos_sorted = np.sort(np.searchsorted(cyc_obs, cyc[found["kick_positions"]]))
    span = int(cyc_obs[-1] - cyc_obs[0]) + 1
    is_prime = _is_prime_table(span)
    run, seq_i, _start, frac = _kick_set_stats_nb(
        cyc_obs, pos_sorted, _KICK_SEQ_FLAT, _KICK_SEQ_OFF, is_prime
    )
    gaps = []
    for a, b in pairwise(pos_sorted):
        dc, dp = int(cyc_obs[b] - cyc_obs[a]), int(b - a)
        gaps.append(dc if dc == dp else -1)
    ndraw = int(nshuffle)
    nstep = cyc_obs.size
    sels = np.empty((ndraw, m), dtype=np.int64)
    for d in range(ndraw):
        sels[d] = np.sort(rng.permutation(nstep)[:m])
    runs, fracs = _kick_null_nb(cyc_obs, sels, _KICK_SEQ_FLAT, _KICK_SEQ_OFF, is_prime)
    label = None
    if run >= sequence.KICK_MIN_RUN and seq_i >= 0:
        label = sequence.kick_gap_match(gaps).get("label")
    out.update(
        {
            "kick_gaps": gaps,
            "kick_label": label,
            "kick_run": int(run),
            "p_kick_catalog": shuffle_pvalue(runs, float(run)),
            "kick_prime_frac": None if not np.isfinite(frac) else float(frac),
            "p_kick_prime": None if not np.isfinite(frac) else shuffle_pvalue(fracs, float(frac)),
            "kick_draws": ndraw,
        }
    )
    return out


KICK_CATALOG_VERDICT = "sparse kick code: gaps match"
KICK_PRIME_VERDICT = "kick spacings favour prime lags"


def kick_flags(res: dict) -> list[str]:
    """The verdict fragments the kick route writes, from its stored fields."""
    flags = []
    if not res.get("kick_one_sided") or int(res.get("kick_n") or 0) < KICK_MIN:
        return flags
    draws = res.get("kick_draws")
    if (
        below_gate(res.get("p_kick_catalog"), draws, ACF_P_FOLLOWUP, "kick gaps")
        and int(res.get("kick_run") or 0) >= sequence.KICK_MIN_RUN
    ):
        flags.append(f"{KICK_CATALOG_VERDICT} {res.get('kick_label') or 'a catalog'}")
    if below_gate(res.get("p_kick_prime"), draws, ACF_P_FOLLOWUP, "kick prime lags"):
        flags.append(KICK_PRIME_VERDICT)
    return flags


def analyse_oc(
    timings: Sequence[Timing],
    sample_period: float,
    pulsation_period: float,
    ground_based: bool,
    nshuffle: int = 2000,
    rng: np.random.Generator | None = None,
) -> tuple[dict, np.ndarray]:
    rng = rng or np.random.default_rng(42)
    tc = np.array([c.t_center for c in timings])
    oc = np.array([c.shift_days for c in timings])
    err = np.array([c.shift_err_days for c in timings])
    n = len(oc)
    if n < 12:
        return {"n_series": n, "note": "too few O-C points for statistics"}, oc
    idx = np.array([c.index for c in timings], float)
    A = np.vstack([np.ones(n), idx, idx**2]).T
    cc, *_ = np.linalg.lstsq(A, oc, rcond=None)
    oc_d = oc - A @ cc

    maxlag = lb_maxlag(n)
    r = acf(oc_d, maxlag)
    pr = primes_upto(maxlag)
    comp = np.setdiff1d(np.arange(2, maxlag + 1), pr)
    prime_stat = r[pr].mean() - r[comp].mean() if len(comp) else np.nan
    max_abs = np.max(np.abs(r[1:]))
    q, p_lb_chi2 = ljung_box(oc_d, maxlag)
    pe_obs = permutation_entropy(oc_d)
    # Rank-based, so its null is a function of n alone and is tabulated once
    # per length rather than redrawn inside the loop below.
    null_pe = permutation_entropy_null(n, nshuffle)
    # The O-C read as a clock-locked bit word, and the null for it. The bits
    # are the *steps* of the O-C, one per cycle -- a triggered pulse resets the
    # clock, so a phase code is a staircase and its word is the first
    # difference (sequence.py, research/modulation-physics.md section 2). A
    # shuffle keeps the values and destroys the cycle order; each draw is then
    # differenced and split exactly the way the observation is, so the draws
    # search the same statistic and only the clock is under test.
    cyc_idx = np.array([c.index for c in timings], dtype=np.int64)
    order = np.argsort(cyc_idx, kind="stable")
    inc_obs, inc_cyc, pairs = code_increments(oc_d[order], cyc_idx[order])
    # The steps before the clip: a kick is exactly the step the clip would
    # have pulled in, so the kick route reads the raw series.
    inc_raw = oc_d[order][pairs + 1] - oc_d[order][pairs]
    # Where each step's two points sit in the unsorted series, so a draw can
    # be differenced without being re-sorted.
    lo_idx, hi_idx = order[pairs], order[pairs + 1]
    # Fold first: the statistic averages the steps per residue over every
    # word length and bit period before splitting them, so a word repeated r
    # times is read at sqrt(r) the per-cycle signal-to-noise. The draws run
    # the same search, which is what puts the look-elsewhere over L and k
    # inside the p-value. `_otsu_nb` asks the content-free question -- are
    # the steps quantised? -- which needs no repeat at all.
    maxlen = min(sequence.SEQ_MAX_BITS, inc_obs.size // 2)
    # The bit-period grid grows with the series (sequence.bit_cycles_for), so
    # a long series is searched for slow codes too; the draws below search
    # the same grid, which is what charges it.
    ks = np.array(sequence.bit_cycles_for(inc_obs.size), dtype=np.int64)
    decodable = inc_obs.size >= 8 and maxlen >= 4
    if decodable:
        mdl_obs, mdl_len, mdl_k = _fold_cost_nb(inc_obs, inc_cyc, maxlen, ks)
        decodable = mdl_len > 0
    else:
        mdl_obs, mdl_len, mdl_k = float("nan"), 0, 0
    quant_obs = _otsu_nb(inc_obs) if inc_obs.size >= 12 else float("nan")
    quantisable = bool(np.isfinite(quant_obs))
    # The draws run in blocks so the index matrix stays small on a star with a
    # long O-C series; the permutations come off `rng` in the same order a
    # draw-at-a-time loop took them, so the null does not depend on the block
    # size. The Ljung-Box Q comes free from the ACF each draw already
    # computes, which is what calibrates it by permutation rather than by a
    # chi2 tail that is not valid at maxlag ~ n/3.
    parts: dict[str, list[np.ndarray]] = {k: [] for k in ("max", "lb", "prime", "mdl", "quant")}
    block = max(1, min(int(nshuffle), NULL_BLOCK_DRAWS, max(1, NULL_BLOCK_CELLS // max(n, 1))))
    drawn = 0
    while drawn < int(nshuffle):
        nb = min(block, int(nshuffle) - drawn)
        perms = np.empty((nb, n), dtype=np.int64)
        for b in range(nb):
            perms[b] = rng.permutation(n)
        for key, arr in zip(
            ("max", "lb", "prime", "mdl", "quant"),
            _null_loop_nb(
                oc_d,
                perms,
                lo_idx,
                hi_idx,
                inc_cyc,
                pr,
                comp,
                int(maxlag),
                int(maxlen),
                ks,
                bool(decodable),
                bool(quantisable),
                float(sequence.STEP_CLIP_SIGMA),
                float(sequence.STEP_CLIP_MAX_FRACTION),
                int(sequence.KICK_MIN_STEPS),
                float(sequence.KICK_SPREAD),
            ),
        ):
            parts[key].append(arr)
        drawn += nb
    null = {k: (np.concatenate(v) if v else np.zeros(0)) for k, v in parts.items()}
    null_max, null_prime, null_lb = null["max"], null["prime"], null["lb"]
    null_mdl, null_quant = null["mdl"], null["quant"]
    # Cheaper is more structured, so the tail runs the other way. The draws
    # search word length and bit period exactly as the observation does,
    # which is what puts the look-elsewhere over L and k inside the p-value
    # instead of beside it.
    p_decode = shuffle_pvalue(np.array(null_mdl), mdl_obs, upper=False) if decodable else None
    p_quant = shuffle_pvalue(np.array(null_quant), quant_obs) if quantisable else None
    # The receiver a beacon designer expects: the folded steps correlated
    # against every catalog window, charged analytically for every window,
    # rotation, length and bit period tried (matched.py).
    hit = (
        matched.matched_catalog(inc_obs, inc_cyc, int(mdl_len), int(mdl_k))
        if decodable and inc_obs.size >= 2 * matched.MIN_LEN
        else {}
    )
    seq = sequence.score_oc_pairs(
        [(int(c.index), float(v) * 86400.0) for c, v in zip(timings, oc_d)]
    )
    p_max = shuffle_pvalue(null_max, max_abs)
    p_prime = shuffle_pvalue(null_prime, prime_stat) if len(comp) else None
    p_pe = (
        shuffle_pvalue(null_pe, pe_obs, upper=False)
        if pe_obs is not None and len(null_pe)
        else None
    )
    p_lb = shuffle_pvalue(null_lb, q)
    mad = stats.median_abs_deviation(oc_d, scale="normal")
    jumps = [int(timings[i].index) for i in np.flatnonzero(np.abs(oc_d) > 5 * mad)]
    sec = secondary_frequencies(tc, oc_d, sample_period, pulsation_period, ground_based)
    # The kicks: rare one-sided steps, read as a message in their spacing.
    kicks = kick_word(inc_raw, inc_cyc, nshuffle, rng)

    flags = []
    if p_lb < ACF_P_FOLLOWUP or p_max < ACF_P_FOLLOWUP:
        flags.append("significant autocorrelation")
    if p_pe is not None and p_pe < ACF_P_FOLLOWUP:
        flags.append("low permutation entropy")
    real_sec = [s for s in sec if not s["alias_suspect"]]
    if real_sec:
        flags.append("coherent secondary period(s)")
    elif sec:
        flags.append("secondary period(s) coincide with sampling aliases -> ignore")
    # A code is two levels on a fixed clock. The p-value says the clock is
    # real; `seq_square` -- stage 2's own predicate -- says the levels are
    # levels and not a smooth wave, which is the whole difference between a
    # phase code and Blazhko modulation.
    if p_decode is not None and p_decode < ACF_P_FOLLOWUP and seq.get("seq_square"):
        flags.append("coherent two-level phase code")
    # Steps of two sizes and nothing in between: what a code's clock does and
    # a random walk's does not.
    # Both need the two levels to beat a sinusoid, as the decode route does:
    # smooth modulation has an arcsine step distribution that Otsu reads as
    # two peaks, and folded at some length it correlates at 0.9 with a square
    # catalog window. On real RR Lyrae 7-9% fired without this guard.
    if p_quant is not None and p_quant < ACF_P_FOLLOWUP and seq.get("seq_square"):
        flags.append("quantised steps")
    # One catalog route. `seq_p` is the matched filter's own look-elsewhere
    # p-value, and the window it names has to beat a sinusoid on the steps
    # (sequence.is_square_match): a smooth wave that happens to correlate
    # with a catalog window is Blazhko with a name, not a code. The window
    # test needs only the two repeats the filter itself needs, so every named
    # word has one, where `seq_square` has none below three.
    if float(seq.get("seq_p") or 1.0) < score.SEQ_P_GATE and seq.get("seq_window_square"):
        flags.append(CATALOG_VERDICT)
    flags.extend(kick_flags(kicks))
    return {
        "n_series": n,
        **kicks,
        "oc_rms_seconds": float(np.std(oc_d) * 86400),
        "median_timing_err_seconds": float(np.nanmedian(err) * 86400),
        "secular_fit_coeffs": [float(v) for v in cc],
        "shuffle_draws": int(null_max.size),
        "max_abs_acf": float(max_abs),
        "p_max_acf_shuffle": p_max,
        "prime_lag_stat": None if not np.isfinite(prime_stat) else float(prime_stat),
        "p_prime_lag_shuffle": p_prime,
        "ljung_box_Q": q,
        "p_ljung_box": p_lb,
        "p_ljung_box_chi2": p_lb_chi2,
        "log10_p_ljung_box_chi2": chi2_log10_sf(q, maxlag),
        "ljung_box_maxlag": int(maxlag),
        "permutation_entropy": pe_obs,
        "p_entropy_shuffle": p_pe,
        "decode_mdl": None if not decodable else float(mdl_obs),
        "decode_bit_len": int(mdl_len),
        "decode_bit_cycles": int(mdl_k),
        "decode_model": sequence.DEFAULT_CODE_MODEL,
        "p_decode_shuffle": p_decode,
        "quantised_stat": None if not quantisable else float(quant_obs),
        "p_quantised_shuffle": p_quant,
        "p_catalog_matched": None if not hit else float(hit["p"]),
        "catalog_matched_label": None if not hit else hit.get("label"),
        "catalog_matched_rho": None if not hit else hit.get("rho"),
        "catalog_matched_len": None if not hit else hit.get("len"),
        "catalog_matched_k": None if not hit else hit.get("k"),
        # The gate reads these out of the payload, so a --triage rescore can
        # re-apply it. The full row goes to its own table, and rides along here
        # under `sequence_row` -- like `phase_fold`, a derived field the store
        # writes elsewhere and `store.payload_fields` keeps out of the payload.
        # The decode is the most expensive thing either writer does, and both
        # used to repeat it on the frame after this function had already done
        # it; handing the row over means the word is read exactly once.
        "sequence_row": seq,
        "seq_p": seq.get("seq_p"),
        "seq_bit_len": seq.get("seq_bit_len"),
        "seq_squareness": seq.get("seq_squareness"),
        "seq_square": bool(seq.get("seq_square")),
        "seq_window_squareness": seq.get("seq_window_squareness"),
        "seq_window_square": bool(seq.get("seq_window_square")),
        "seq_label": seq.get("seq_label"),
        "jump_indices_gt5sigma": jumps,
        "secondary_periods": sec,
        "verdict": "; ".join(flags) if flags else "no evidence of non-random structure",
    }, oc_d


# ============================================================================
# The per-cycle series of a sparse survey, folded on the word
# ============================================================================

# Fewer observed cycles than this and the fold has nothing to integrate: a
# 4-bit word at one cycle per bit needs each residue seen several times.
EPOCH_MIN_CYCLES = 24
EPOCH_CATALOG_VERDICT = (
    "epoch fold: decoded word matches a catalog sequence, squarer than a sinusoid"
)
EPOCH_CODE_VERDICT = "coherent two-level phase code (epoch fold)"

# The fields `analyse_epoch_fold` writes, so a reader can tell "not run" from
# "run and found nothing".
EPOCH_FIELDS = (
    "epoch_n_cycles",
    "epoch_n_epochs",
    "epoch_oc_rms_seconds",
    "epoch_shuffle_draws",
    "epoch_decode_mdl",
    "epoch_decode_bit_len",
    "epoch_decode_bit_cycles",
    "p_epoch_decode_shuffle",
    "epoch_seq_p",
    "epoch_seq_label",
    "epoch_seq_bit_len",
    "epoch_seq_squareness",
    "epoch_seq_square",
    "epoch_seq_window_squareness",
    "epoch_seq_window_square",
    "epoch_code",
)


def analyse_epoch_fold(
    timings: Sequence[Timing], nshuffle: int = 2000, rng: np.random.Generator | None = None
) -> tuple[dict, np.ndarray]:
    """Read a per-cycle timing series with gaps as a repeated bit word.

    `timings` are one per observed cycle (`cycle_timings`); most cycles are
    missing and almost no two adjacent ones are present, so the staircase
    cannot be differenced. Instead the detrended levels are folded on every
    word length and bit period the series can hold (`_level_fold_cost_tab_nb`,
    the `fold` model) and the bits are read from the differences of the
    residue means; the cheapest description over (L, k) is compared against
    the same search on `nshuffle` permutations of the levels among the
    observed cycles, so the choice of word length and bit period is charged
    inside the p-value. The winning fold is then described, tested for
    squareness and named through `sequence.decode_bit_code(model="fold")`.

    Returns the `epoch_*` fields and the detrended series in days.
    """
    rng = rng or np.random.default_rng(43)
    n = len(timings)
    out: dict = dict.fromkeys(EPOCH_FIELDS)
    out["epoch_n_cycles"] = n
    out["epoch_n_epochs"] = int(sum(int(c.npts) for c in timings))
    if n < EPOCH_MIN_CYCLES:
        out["epoch_note"] = "too few observed cycles for the epoch fold"
        return out, np.array([c.shift_days for c in timings])
    cyc = np.array([c.index for c in timings], dtype=np.int64)
    order = np.argsort(cyc, kind="stable")
    cyc = cyc[order]
    oc = np.array([c.shift_days for c in timings])[order]
    idx = cyc.astype(float)
    A = np.vstack([np.ones(n), idx, idx**2]).T
    cc, *_ = np.linalg.lstsq(A, oc, rcond=None)
    oc_d = oc - A @ cc
    out["epoch_oc_rms_seconds"] = float(np.std(oc_d) * 86400)
    maxlen = min(sequence.SEQ_MAX_BITS, n // 2)
    ks = np.array(sequence.bit_cycles_for(n), dtype=np.int64)
    tab, tcnt, tlens, tkk = _fold_index_table(cyc, n, maxlen, ks)
    if tab.shape[0] == 0:
        out["epoch_note"] = "no word length folds on these cycles"
        return out, oc_d[np.argsort(order)]
    mdl_obs, mdl_len, mdl_k = _level_fold_cost_tab_nb(oc_d, tab, tcnt, tlens, tkk)
    if mdl_len == 0:
        out["epoch_note"] = "no two-level word on these cycles"
        return out, oc_d[np.argsort(order)]
    null: list[np.ndarray] = []
    block = max(1, min(int(nshuffle), NULL_BLOCK_DRAWS, max(1, NULL_BLOCK_CELLS // max(n, 1))))
    drawn = 0
    while drawn < int(nshuffle):
        nb = min(block, int(nshuffle) - drawn)
        perms = np.empty((nb, n), dtype=np.int64)
        for b in range(nb):
            perms[b] = rng.permutation(n)
        null.append(_epoch_null_loop_nb(oc_d, perms, tab, tcnt, tlens, tkk))
        drawn += nb
    null_mdl = np.concatenate(null) if null else np.zeros(0)
    out.update(
        {
            "epoch_shuffle_draws": int(null_mdl.size),
            "epoch_decode_mdl": float(mdl_obs),
            "epoch_decode_bit_len": int(mdl_len),
            "epoch_decode_bit_cycles": int(mdl_k),
            "p_epoch_decode_shuffle": shuffle_pvalue(null_mdl, float(mdl_obs), upper=False),
        }
    )
    code = None
    try:
        code = sequence.decode_bit_code(
            [float(v) * 86400.0 for v in oc_d],
            [int(c) for c in cyc],
            model="fold",
            fold_at=(int(mdl_len), int(mdl_k)),
        )
    except Exception:  # a word that cannot be described is no word
        code = None
    fields = sequence.fields_from_code(code, "")
    out.update(
        {
            "epoch_seq_p": fields.get("seq_p"),
            "epoch_seq_label": fields.get("seq_label"),
            "epoch_seq_bit_len": fields.get("seq_bit_len"),
            "epoch_seq_squareness": fields.get("seq_squareness"),
            "epoch_seq_square": bool(fields.get("seq_square")),
            "epoch_seq_window_squareness": fields.get("seq_window_squareness"),
            "epoch_seq_window_square": bool(fields.get("seq_window_square")),
            "epoch_code": None if not code else json_safe(code),
        }
    )
    # Back in the caller's order.
    inv = np.empty(n, dtype=np.int64)
    inv[order] = np.arange(n)
    return out, oc_d[inv]


def epoch_flags(res: dict) -> list[str]:
    """The verdict fragments the epoch-fold routes write, from their fields.

    The same two gates the dense decode carries, asked of the per-cycle fold:
    the description-length p-value with the fold's own word squarer than a
    sinusoid, and a named catalog word at the sequence gate with its window
    squarer than a sinusoid.
    """
    flags = []
    draws = res.get("epoch_shuffle_draws")
    if below_gate(res.get("p_epoch_decode_shuffle"), draws, ACF_P_FOLLOWUP, "epoch fold") and bool(
        res.get("epoch_seq_square")
    ):
        flags.append(EPOCH_CODE_VERDICT)
    p_seq = res.get("epoch_seq_p")
    if (
        p_seq is not None
        and float(p_seq) < score.SEQ_P_GATE
        and bool(res.get("epoch_seq_window_square"))
    ):
        flags.append(EPOCH_CATALOG_VERDICT)
    return flags


def epoch_frame(timings: Sequence[Timing], oc_d: np.ndarray) -> pd.DataFrame:
    """The per-cycle series as the store's `stage1_epoch_oc` rows."""
    df = pd.DataFrame([asdict(c) for c in timings]).rename(columns={"index": "cycle_index"})
    df["oc_detrended_days"] = np.asarray(oc_d, dtype=float)
    return df.drop(columns=["amp_ratio"], errors="ignore")


def _add_flags(verdict: str, flags: Sequence[str]) -> str:
    if not flags:
        return verdict
    joined = "; ".join(flags)
    return joined if verdict.startswith("no evidence") else verdict + "; " + joined


# ============================================================================
# Compact phase fold for the HTML viewer
# ============================================================================


def pack_phase_fold(
    t: np.ndarray,
    f: np.ndarray,
    period: float,
    coef: np.ndarray,
    off: dict[str, float],
    t0: float,
    nharm: int,
    keep: np.ndarray,
    cam: np.ndarray | None,
    max_pts: int = 2500,
    n_template: int = 400,
) -> dict:
    """Downsampled phase-fold + Fourier template, stored on the stage-1 payload."""
    flux = np.asarray(f, float) - cam_offset_array(cam, off)
    ph = ((np.asarray(t, float) - t0) / period) % 1.0
    ph, flux = ph[keep], flux[keep]
    finite = np.isfinite(ph) & np.isfinite(flux)
    ph, flux = ph[finite], flux[finite]
    if ph.size == 0:
        return None
    order = np.argsort(ph)
    if ph.size > max_pts:
        order = order[np.linspace(0, ph.size - 1, max_pts).astype(int)]
    ph, flux = ph[order], flux[order]
    pp = np.linspace(0.0, 1.0, n_template)
    tf = fourier_design(pp, nharm) @ coef
    return {
        "phase": [round(float(x), 5) for x in ph],
        "flux": [round(float(x), 6) for x in flux],
        "template_phase": [round(float(x), 5) for x in pp],
        "template_flux": [round(float(x), 6) for x in tf],
    }


# ============================================================================
# Synthetic self-test
# ============================================================================

# Minimal-binary encodings of eight successive catalog numbers, one sequence per class.
# Timing advance per set bit is dt/P = 1e-3. That is a working choice between the two
# numbers Learned et al. (2008) actually give -- a 0.05-cycle advance in their worked
# example, and an energy budget that implies ~1e-6 at unit gain -- not a figure from
# the paper (research/modulation-physics.md section 1). The bit word is tiled
# SYNTH_INJECT_REPS times so the O-C carries that many copies of it.
#
# The waveform is a staircase, not a square wave: a set bit triggers the pulse early
# and the clock stays reset, so the shift at cycle c is dt times the number of set
# bits up to and including c (section 2 of the same note, `code_shift`).
SYNTH_INJECT_FRAC = 1e-3
# Ground-based sampling cannot see a 1e-3 shift, and pretending otherwise makes
# the sparse control a guaranteed failure rather than a test: ~0.015 mag per
# epoch on a 0.57 d RRab is ~470 s of per-epoch timing noise (~165 s per block)
# against a 49 s injection. The sparse control therefore injects 10x more. It
# exercises the sparse machinery; it is not a claim that the sparse arm reaches
# the Learned et al. 1e-3 trigger scaling.
SYNTH_INJECT_FRAC_SPARSE = 1e-2
SYNTH_INJECT_REPS = 5
# Per-point photometric scatter for the synthetic TESS light curves, chosen so
# the injected dt/P = 1e-3 shift lands well above the per-cycle timing error it
# produces. This is a positive control: at the old scatter the injection sat at
# ~0.8 sigma per cycle and the smoke test came back "no evidence", which cannot
# tell a working detector from a broken one. Bright TESS targets reach ~100 ppm
# per 2-min point, so these are still conservative.
#   sigma_t ~ sigma_f / (2 pi A / P * sqrt(N_per_cycle))
SYNTH_TESS_NOISE = {"cepheid": 0.002, "rrlyrae": 0.0015, "bcep": 0.0001}
# Ground-based per-epoch scatter for a bright target (ASAS-SN V~11 reaches
# ~0.015 mag; ZTF a little better). The old 0.15 / 0.09 were faint-end values.
SYNTH_GROUND_NOISE = {"asassn": 0.015, "ztf": 0.012}


def synth_inject_frac(survey: str) -> float:
    return SYNTH_INJECT_FRAC if is_dense_survey(survey) else SYNTH_INJECT_FRAC_SPARSE


TESS_SECTOR_DAYS = 13.7
TESS_SCIENCE_DAYS = 12.9
SYNTH_INJECT_NUMS = {
    "cepheid": ("primes", (2, 3, 5, 7, 11, 13, 17, 19)),
    "bcep": ("Fibonacci", (1, 1, 2, 3, 5, 8, 13, 21)),
    "rrlyrae": ("counting", (1, 2, 3, 4, 5, 6, 7, 8)),
}


def synth_inject_bits(star_class: str) -> np.ndarray:
    _name, nums = SYNTH_INJECT_NUMS.get(star_class, SYNTH_INJECT_NUMS["cepheid"])
    bits = "".join(bin(n)[2:] for n in nums)
    return np.array([int(b) for b in bits], dtype=int)


def synth_inject_code(cyc: np.ndarray, star_class: str) -> np.ndarray:
    bits = synth_inject_bits(star_class)
    n = int(np.max(cyc)) + 1
    return np.tile(bits, n // len(bits) + 2)[:n]


def code_shift(
    code: np.ndarray,
    cyc: np.ndarray,
    frac: float,
    period: float,
    model: str = sequence.DEFAULT_CODE_MODEL,
) -> np.ndarray:
    """Time shift of every sample for a bit code, in days; negative is early.

    `code[c]` is the bit of cycle c. Under the `staircase` model the shift at
    cycle c is dt times the number of set bits up to and including c -- each
    triggered pulse resets the clock and nothing restores it -- so the O-C is
    the running sum of the word. Under `square` it is dt times the current bit
    alone, the waveform the pipeline assumed before the physics was worked
    through; kept so the two can be compared on the same star.
    """
    if model not in sequence.CODE_MODELS:
        raise ValueError(f"unknown code model {model!r}; one of {sequence.CODE_MODELS}")
    level = np.cumsum(code) if model == "staircase" else np.asarray(code)
    return -float(frac) * float(period) * level[cyc]


# TESS 2-min science cadence: what the real controls are generated at, and
# what --synthetic runs use.
SYNTH_TESS_CADENCE_MIN = 2.0
# Optional ceiling on how many samples one pulsation cycle gets. None in
# production. The suite sets it -- see the note in tests/conftest.py -- because
# at long periods 2-min sampling puts thousands of points in a cycle and every
# per-cycle fit pays for all of them. Expressed per cycle rather than as a
# cadence so it binds on the 5.4-day Cepheid and leaves the 0.2-day beta Cep,
# which has no points to spare, exactly as it was.
SYNTH_MAX_PTS_PER_CYCLE: float | None = None


def synth_tess_times(
    t0: float, period: float, star_class: str, cadence_min: float | None = None
) -> np.ndarray:
    """TESS 2-min sampling long enough for SYNTH_INJECT_REPS copies of the bit word."""
    cadence_min = SYNTH_TESS_CADENCE_MIN if cadence_min is None else cadence_min
    if SYNTH_MAX_PTS_PER_CYCLE:
        cadence_min = max(cadence_min, float(period) * 1440.0 / SYNTH_MAX_PTS_PER_CYCLE)
    nbits = max(1, len(synth_inject_bits(star_class)))
    n_cycles = SYNTH_INJECT_REPS * nbits
    span = n_cycles * float(period) * (TESS_SECTOR_DAYS / TESS_SCIENCE_DAYS) + float(period)
    t = np.arange(0.0, span, cadence_min / 1440.0) + float(t0)
    return t[((t - t.min()) % TESS_SECTOR_DAYS) < TESS_SCIENCE_DAYS]


def synthetic(
    survey: str = "tess",
    inject: bool = False,
    seed: int = 1,
    period: float | None = None,
    star_class: str = "cepheid",
) -> LightCurveData:
    rng = np.random.default_rng(seed)
    if star_class == "bcep":
        return synthetic_bcep(inject=inject, seed=seed, period=period)
    if star_class == "rrlyrae":
        return synthetic_rrlyrae(survey=survey, inject=inject, seed=seed, period=period)
    if survey == "tess":
        period = period or 5.366
        t = synth_tess_times(1000.0, period, star_class)
        noise, cam, label = SYNTH_TESS_NOISE["cepheid"], None, "BTJD"
    elif survey == "ztf":
        period = period or 0.5668
        days = np.arange(0, 7 * 365.25)
        nights = days[
            ((days % 365.25) < 220) & (rng.random(days.size) < 0.45)
        ]  # ~2-3 night cadence
        t = np.sort(nights + 0.8 + 0.3 * rng.random(nights.size)) + 8200.0
        noise, label = SYNTH_GROUND_NOISE["ztf"], "HJD-2450000"
        cam = rng.choice(["oidA", "oidB"], size=t.size, p=[0.7, 0.3])  # overlapping fields
    else:
        period = period or 0.5668
        days = np.arange(0, 8 * 365.25)
        nights = days[((days % 365.25) < 240) & (rng.random(days.size) < 0.7)]
        t = (
            np.sort(
                np.concatenate([n + 0.9 + 0.25 * rng.random(rng.integers(1, 4)) for n in nights])
            )
            + 8000.0
        )
        noise, label = SYNTH_GROUND_NOISE["asassn"], "JD-2450000"
        cam = rng.choice(["bd", "be", "bf"], size=t.size)
    amps, phs = [0.25, 0.10, 0.05, 0.025, 0.012], [0.0, 0.9, 1.7, 2.5, 3.2]
    if survey == "tess":
        amps = [a * 0.4 for a in amps]
    cyc = np.floor((t - t.min()) / period).astype(int)
    shifts = np.zeros_like(t)
    if inject:
        code = synth_inject_code(cyc, star_class)
        shifts = code_shift(code, cyc, synth_inject_frac(survey), period)
    phase = (t - shifts - t.min()) / period
    f = 1.0 + sum(
        a * np.sin(2 * np.pi * (k + 1) * phase + p) for k, (a, p) in enumerate(zip(amps, phs))
    )
    f += rng.normal(0, noise, t.size)
    if cam is not None:
        offs = {"bd": 0.0, "be": 0.03, "bf": -0.02, "oidA": 0.0, "oidB": 0.015}
        f += np.array([offs[c] for c in cam])
    return LightCurveData(t, f, np.full_like(f, noise), cam, survey, label)


def synthetic_rrlyrae(
    survey: str = "tess", inject: bool = False, seed: int = 1, period: float | None = None
) -> LightCurveData:
    """Asymmetric RRab analogue: P~0.57 d, steep rise, TESS 2-min covering five code repeats."""
    rng = np.random.default_rng(seed)
    period = period or 0.5668
    if survey == "tess":
        t = synth_tess_times(3000.0, period, "rrlyrae")
        noise, cam, label = SYNTH_TESS_NOISE["rrlyrae"], None, "BTJD"
    elif survey == "ztf":
        days = np.arange(0, 7 * 365.25)
        nights = days[((days % 365.25) < 220) & (rng.random(days.size) < 0.45)]
        t = np.sort(nights + 0.8 + 0.3 * rng.random(nights.size)) + 8200.0
        noise, label = SYNTH_GROUND_NOISE["ztf"], "HJD-2450000"
        cam = rng.choice(["oidA", "oidB"], size=t.size, p=[0.7, 0.3])
    else:
        days = np.arange(0, 8 * 365.25)
        nights = days[((days % 365.25) < 240) & (rng.random(days.size) < 0.7)]
        t = (
            np.sort(
                np.concatenate([n + 0.9 + 0.25 * rng.random(rng.integers(1, 4)) for n in nights])
            )
            + 8000.0
        )
        noise, label = SYNTH_GROUND_NOISE["asassn"], "JD-2450000"
        cam = rng.choice(["bd", "be", "bf"], size=t.size)
    amps, phs = [0.22, 0.11, 0.055, 0.028, 0.014], [0.0, 0.45, 1.0, 1.6, 2.2]
    cyc = np.floor((t - t.min()) / period).astype(int)
    shifts = np.zeros_like(t)
    if inject:
        code = synth_inject_code(cyc, "rrlyrae")
        shifts = code_shift(code, cyc, synth_inject_frac(survey), period)
    phase = (t - shifts - t.min()) / period
    f = 1.0 + sum(
        a * np.sin(2 * np.pi * (k + 1) * phase + p) for k, (a, p) in enumerate(zip(amps, phs))
    )
    f += rng.normal(0, noise, t.size)
    if cam is not None:
        offs = {"bd": 0.0, "be": 0.03, "bf": -0.02, "oidA": 0.0, "oidB": 0.015}
        f += np.array([offs[c] for c in cam])
    return LightCurveData(t, f, np.full_like(f, noise), cam, survey, label)


def synthetic_bcep(
    inject: bool = False,
    seed: int = 1,
    period: float | None = None,
    cadence_min: float | None = None,
) -> LightCurveData:
    """Multiperiodic beta Cep analogue: dominant radial mode P~0.2 d (4.8 h, BW Vul-like
    but low amplitude), two independent non-radial modes at incommensurate
    frequencies, TESS 2-min cadence covering five code repeats. The beat between the modes
    is exactly what would masquerade as O-C structure without prewhitening."""
    rng = np.random.default_rng(seed)
    period = period or 0.2010
    t = synth_tess_times(2000.0, period, "bcep", cadence_min=cadence_min)
    f0 = 1.0 / period
    cyc = np.floor((t - t.min()) / period).astype(int)
    shifts = np.zeros_like(t)
    if inject:
        code = synth_inject_code(cyc, "bcep")
        shifts = code_shift(code, cyc, SYNTH_INJECT_FRAC, period)
    phase = (t - shifts - t.min()) / period
    # dominant mode with a mildly non-sinusoidal shape
    f = 1.0 + 0.020 * np.sin(2 * np.pi * phase) + 0.004 * np.sin(4 * np.pi * phase + 0.8)
    # independent modes (frequency ratios well away from integers)
    f += 0.006 * np.sin(2 * np.pi * (1.3137 * f0) * (t - t.min()) + 0.3)
    f += 0.003 * np.sin(2 * np.pi * (0.7261 * f0) * (t - t.min()) + 2.1)
    noise = SYNTH_TESS_NOISE["bcep"]
    f += rng.normal(0, noise, t.size)
    return LightCurveData(t, f, np.full_like(f, noise), None, "tess", "BTJD")


def synthetic_from_name(name: str) -> LightCurveData | None:
    """'SYNTH cepheid' / 'SYNTH bcep' / 'SYNTH rrlyrae' (optional 'null', survey tag)."""
    low = name.lower()
    sv = "asassn" if "asassn" in low else "ztf" if "ztf" in low else "tess"
    if "rrlyrae" in low or "rr lyr" in low:
        cls = "rrlyrae"
    elif "bcep" in low:
        cls = "bcep"
    else:
        cls = "cepheid"
    inject = "null" not in low
    return synthetic(sv, inject, star_class=cls)


# ============================================================================
# Driver
# ============================================================================


def points_per_covered_cycle(t: np.ndarray, period: float) -> float:
    """Median number of points in cycles that actually contain data (ignores gaps),
    so long multi-sector baselines don't get misrouted to sparse mode."""
    cyc = np.floor((t - t.min()) / period).astype(int)
    # A non-finite time casts to a large negative int, and `np.bincount`
    # answers that with "'list' argument must have no negative elements" --
    # true, and no help at all in finding the cadence that did it. The loaders
    # drop those now; this is the backstop that names the invariant.
    cyc = cyc[cyc >= 0]
    if cyc.size == 0:
        return 0.0
    counts = np.bincount(cyc)
    counts = counts[counts > 0]
    return float(np.median(counts)) if counts.size else 0.0


def refine_period(
    t: np.ndarray,
    f: np.ndarray,
    e: np.ndarray,
    period: float,
    coef: np.ndarray,
    off: dict[str, float],
    t0: float,
    nharm: int,
    cam: np.ndarray | None = None,
    niter: int = 4,
) -> tuple[float, np.ndarray, dict[str, float]]:
    """Linearised period refinement. A period error dP makes the residual
    r(t) ~ -slope(t) * dP * (t - t0)/P; solve for dP by weighted least squares and
    re-fit the template. Without this, the residual phase drift over a multi-year
    ground-based baseline masquerades as cycle-to-cycle jitter."""
    for _ in range(niter):
        fk = f - cam_offset_array(cam, off)
        x = (t - t0) / period
        r = fk - fourier_model(t, coef, period, t0, nharm)
        g = -fourier_slope(t, coef, period, t0, nharm) * x
        w = 1 / e**2
        dP = np.sum(w * g * r) / max(np.sum(w * g * g), 1e-300)
        if not np.isfinite(dP) or abs(dP) < 1e-12 * period:
            break
        period = period + dP
        coef, off, t0 = fit_fourier(t, f, period, nharm, cam, w, t0=t0)
    return period, coef, off


def prepare(
    lcd: LightCurveData, cls: StarClass | None = None, period_hint: float | None = None
) -> dict:
    """Period, template, clipped arrays. Shared with stage 2.

    For classes with cls.prewhiten > 0 the flux is cleaned of that many
    independent modes (see prewhiten_modes) and the template re-fitted; 'f' and
    'f_full' are then the cleaned flux (clipped / full length).

    period_hint is the Gaia SOS period from the store when one exists; the
    Lomb-Scargle search is then a ±2% window around that value."""
    cls = cls or STAR_CLASSES[DEFAULT_STAR_CLASS]
    period = find_period(lcd.t, lcd.f, pmin=cls.pmin, pmax=cls.pmax, hint=period_hint)
    coef, off, t0, nharm, keep = robust_template(
        lcd.t, lcd.f, lcd.e, period, lcd.cam, max_nharm=cls.max_nharm
    )
    period, coef, off = refine_period(
        lcd.t[keep],
        lcd.f[keep],
        lcd.e[keep],
        period,
        coef,
        off,
        t0,
        nharm,
        None if lcd.cam is None else lcd.cam[keep],
    )
    coef, off, t0, nharm, keep = robust_template(
        lcd.t, lcd.f, lcd.e, period, lcd.cam, nharm=nharm, t0=t0
    )
    f_full = lcd.f.copy()
    modes = []
    if cls.prewhiten > 0 and (cls.prewhiten_ground or is_dense_survey(lcd.survey)):
        r = lcd.f - cam_offset_array(lcd.cam, off) - fourier_model(lcd.t, coef, period, t0, nharm)
        # amplitude floor: S/N >= PREWHITEN_MIN_SNR against the residual
        # spectrum's noise (not a fixed fraction of the dominant mode, which
        # left the millimagnitude 0.61/0.68 modes in); search the class period
        # window (independent modes, not the photometric Nyquist).
        nyq = 0.5 / max(float(np.median(np.diff(np.sort(lcd.t[keep])))), 1e-6)
        _mode_model, modes = prewhiten_modes(
            lcd.t[keep],
            r[keep],
            lcd.e[keep],
            period,
            nmax=cls.prewhiten,
            min_snr=PREWHITEN_MIN_SNR,
            min_amplitude=PREWHITEN_MIN_FRACTION * template_peak_to_peak(coef, nharm),
            fmax=min(nyq, 1.0 / cls.pmin),
        )
        if modes:
            full_model = np.zeros_like(lcd.f)
            for m in modes:
                A = np.vstack(
                    [
                        np.sin(2 * np.pi * m["frequency_per_day"] * lcd.t),
                        np.cos(2 * np.pi * m["frequency_per_day"] * lcd.t),
                    ]
                ).T
                w = 1.0 / lcd.e[keep]
                c, *_ = np.linalg.lstsq(A[keep] * w[:, None], r[keep] * w, rcond=None)
                full_model += A @ c
                r = r - A @ c
            f_full = lcd.f - full_model
            coef, off, t0, nharm, keep = robust_template(
                lcd.t, f_full, lcd.e, period, lcd.cam, nharm=nharm, t0=t0
            )
            period, coef, off = refine_period(
                lcd.t[keep],
                f_full[keep],
                lcd.e[keep],
                period,
                coef,
                off,
                t0,
                nharm,
                None if lcd.cam is None else lcd.cam[keep],
            )
            coef, off, t0, nharm, keep = robust_template(
                lcd.t, f_full, lcd.e, period, lcd.cam, nharm=nharm, t0=t0
            )
    # Amplitude, not just level: see camera_scales. Only worth the two extra
    # fits when the cameras actually disagree, so a single-camera TESS star and
    # a well-matched pair both skip it.
    scales = camera_scales(lcd.t, f_full, lcd.e, period, coef, off, t0, nharm, lcd.cam)
    if scales and max(abs(s - 1.0) for s in scales.values()) > CAM_SCALE_TOL:
        f_full = apply_camera_scales(f_full, lcd.cam, off, scales, float(coef[0]))
        coef, off, t0, nharm, keep = robust_template(
            lcd.t, f_full, lcd.e, period, lcd.cam, nharm=nharm, t0=t0
        )
        period, coef, off = refine_period(
            lcd.t[keep],
            f_full[keep],
            lcd.e[keep],
            period,
            coef,
            off,
            t0,
            nharm,
            None if lcd.cam is None else lcd.cam[keep],
        )
        coef, off, t0, nharm, keep = robust_template(
            lcd.t, f_full, lcd.e, period, lcd.cam, nharm=nharm, t0=t0
        )
    else:
        scales = {}
    f_all = f_full - cam_offset_array(lcd.cam, off)
    fk = f_all[keep]
    return {
        "period": period,
        "coef": coef,
        "off": off,
        "scales": scales,
        "t0": t0,
        "nharm": nharm,
        "keep": keep,
        "t": lcd.t[keep],
        "f": fk,
        "e": lcd.e[keep],
        # Every point, offsets removed. The timing fits take these: the 5-sigma
        # mask above is against the mean ephemeris, and a cycle that really
        # arrives early is exactly what it clips. The timers clip per cycle
        # against their own shifted model instead.
        "t_all": lcd.t,
        "f_all": f_all,
        "e_all": lcd.e,
        "f_full": f_full,
        "modes": modes,
        "star_class": cls.name,
        "peak_to_peak": template_peak_to_peak(coef, nharm),
        "gaia_period_days": None if period_hint is None else float(period_hint),
    }


def run_target(
    name: str,
    lcd: LightCurveData,
    nshuffle: int,
    cls: StarClass | None = None,
    period_hint: float | None = None,
) -> tuple[dict, pd.DataFrame | None]:
    cls = cls or STAR_CLASSES[DEFAULT_STAR_CLASS]
    P = prepare(lcd, cls, period_hint=period_hint)
    tk, fk, ek, period = P["t"], P["f"], P["e"], P["period"]
    pts_per_cycle = points_per_covered_cycle(tk, period)
    mode = "dense" if pts_per_cycle >= DENSE_MIN_PTS_PER_CYCLE else "sparse"
    # All points, not the template's kept mask: see prepare()'s note. The
    # timers do their own clipping against the shifted model.
    ta, fa, ea = P["t_all"], P["f_all"], P["e_all"]
    if mode == "dense":
        timings, K = dense_timings(ta, fa, ea, period, P["coef"], P["t0"], P["nharm"]), 1
    else:
        timings, K = sparse_timings(ta, fa, ea, period, P["coef"], P["t0"], P["nharm"])
    base = {
        "target": name,
        "survey": lcd.survey,
        "star_class": cls.name,
        "mode": mode,
        "cycles_per_block": K,
        "period_days": float(period),
        "n_points": len(tk),
        "points_per_cycle": float(pts_per_cycle),
        "baseline_days": float(tk.max() - tk.min()),
        "nharm": P["nharm"],
        "template_peak_to_peak": float(P["peak_to_peak"]),
        "prewhitened_modes": P["modes"] or None,
        "gaia_period_days": P.get("gaia_period_days"),
        "camera_offsets": {k: float(v) for k, v in P["off"].items()} or None,
        "camera_scales": {k: float(v) for k, v in (P.get("scales") or {}).items()} or None,
    }
    if P["peak_to_peak"] < cls.min_amplitude:
        return {
            **base,
            "n_series": 0,
            "note": (
                f"template peak-to-peak {P['peak_to_peak']:.4f} below "
                f"{cls.name} floor {cls.min_amplitude}; skip"
            ),
        }, None
    if len(timings) < 12:
        return {
            **base,
            "n_series": len(timings),
            "note": "fewer than 12 O-C points; skip statistics",
        }, None
    stats_d, oc_d = analyse_oc(
        timings, K * period, period, not is_dense_survey(lcd.survey), nshuffle=nshuffle
    )
    extra: dict = {}
    if mode == "sparse":
        jt = jitter_test(tk, fk, ek, period, P["coef"], P["t0"], P["nharm"])
        stats_d.update(jt)
        if (
            below_gate(
                jt["p_jitter_bootstrap"], jt.get("jitter_nboot"), ACF_P_FOLLOWUP, "jitter bootstrap"
            )
            and (jt["jitter_var_t_over_err"] or 0) > JITTER_MIN_SNR
        ):
            stats_d["verdict"] = _add_flags(
                stats_d["verdict"],
                [f"excess cycle-to-cycle phase jitter ({jt['jitter_sigma_t_seconds']:.0f} s rms)"],
            )
        # The per-cycle series, folded on the word: the content the blocks
        # average away (cycle_timings, analyse_epoch_fold).
        dt, dterr, cyc, use = epoch_offsets(ta, fa, ea, period, P["coef"], P["t0"], P["nharm"])
        ctim = cycle_timings(ta, dt, dterr, cyc, use)
        ep, ep_d = analyse_epoch_fold(ctim, nshuffle=nshuffle)
        stats_d.update(ep)
        stats_d["verdict"] = _add_flags(stats_d["verdict"], epoch_flags(ep))
        if ctim:
            extra["epoch_oc"] = epoch_frame(ctim, ep_d)
    df = pd.DataFrame([asdict(c) for c in timings])
    df = df.rename(columns={"index": "cycle_index"})
    df["oc_detrended_days"] = oc_d
    fold = pack_phase_fold(
        lcd.t, P["f_full"], period, P["coef"], P["off"], P["t0"], P["nharm"], P["keep"], lcd.cam
    )
    return {**base, **stats_d, "phase_fold": fold, **extra}, df


def load_data(
    name: str,
    survey: str,
    paths: list[str],
    cadence: str,
    band: str,
    cache: lc_cache.CachePolicy | None = None,
    hint: dict | None = None,
) -> LightCurveData | None:
    """The light curve for one (target, survey), from disk where there is one.

    TESS reads `paths` -- the FITS the bulk prefetch already wrote -- and only
    downloads what is missing. Kepler goes through lightkurve's own cache, one
    quarter at a time. ASAS-SN and ZTF have no files to be handed, so
    the equivalent lookup happens one level down, inside the downloaders,
    against `lc_cache`; `cache` is how a run says how fresh an entry has to be.
    """
    if name.upper().startswith("SYNTH"):
        return synthetic_from_name(name if survey == "tess" else f"{name}__{survey}")
    if survey == "kepler":
        return download_kepler(name, cadence=cadence)
    if survey == "tess":
        data = load_lightcurve_files(paths) if paths else None
        return data if data is not None else download_tess(name, cadence=cadence)
    if survey == "ztf":
        return download_ztf(name, band=band, cache=cache, hint=hint)
    return download_asassn(name, band=band, cache=cache)


def process_target(
    name: str,
    survey: str,
    paths: list[str],
    cadence: str,
    band: str,
    nshuffle: int,
    star_class: StarClass | str = DEFAULT_STAR_CLASS,
    period_hint: float | None = None,
    cache: lc_cache.CachePolicy | None = None,
    hint: dict | None = None,
) -> tuple[dict, pd.DataFrame | None]:
    """Top-level so ProcessPoolExecutor can pickle it. Returns (result, oc_df).

    `cache` travels with the work item rather than sitting on the module,
    because Windows spawns its workers and a module-level setting made in the
    parent is simply not there in the child -- the `--prewhiten` bug.
    """
    cls = resolve_cls(star_class)
    try:
        data = load_data(name, survey, paths, cadence, band, cache, hint)
    except Exception as ex:
        return {
            "target": name,
            "survey": survey,
            "star_class": cls.name,
            "error": f"download failed: {ex}",
        }, None
    if data is None:
        return {"target": name, "survey": survey, "star_class": cls.name, "error": "no data"}, None
    try:
        return run_target(name, data, nshuffle, cls, period_hint=period_hint)
    except Exception as ex:
        return {
            "target": name,
            "survey": survey,
            "star_class": cls.name,
            "error": f"analysis failed: {ex}",
        }, None


def print_result(res: dict) -> None:
    vlog(f"=== {res.get('target')} [{res.get('survey')}] ===")
    if res.get("error"):
        vlog(f"  {res['error']}")
        return
    keys = (
        "star_class",
        "mode",
        "period_days",
        "n_series",
        "oc_rms_seconds",
        "p_ljung_box",
        "jitter_sigma_t_seconds",
        "verdict",
        "note",
    )
    if res.get("prewhitened_modes"):
        vlog(
            "  prewhitened: "
            + ", ".join(
                f"{m['period_days']:.5f} d (A={m['amplitude']:.4f})"
                for m in res["prewhitened_modes"]
            )
        )
    vlog("  " + json.dumps({k: res[k] for k in keys if k in res}))


# ---- triage -----------------------------------------------------------------


def _p(val: object) -> float:
    try:
        return 1.0 if val is None else float(val)
    except (TypeError, ValueError):
        return 1.0


# Verdict fragments that are themselves a reason to follow a star up, so a
# rescore of a payload written before a statistic existed still routes on what
# that run actually concluded.
# The verdict the catalog route writes. The bare "decoded word matches a
# catalog sequence" of older payloads was written without the window test and
# is deliberately not in FOLLOWUP_VERDICTS: such a row is re-judged from its
# stored fields by `store.recompute_triage`, not promoted on the text.
CATALOG_VERDICT = "decoded word matches a catalog sequence, squarer than a sinusoid"

FOLLOWUP_VERDICTS = (
    "significant autocorrelation",
    "excess cycle-to-cycle",
    "two-level phase code",
    "catalog sequence, squarer",
    "sparse kick code",
    "favour prime lags",
)


def is_followup(res: dict) -> bool:
    """Does this star go to stage 2?

    Four routes, and they answer different questions. The first three ask
    whether the O-C is *unpredictable*: Ljung-Box and max|ACF| against a
    shuffle null, and in sparse mode the per-epoch jitter. The fourth asks
    whether it is *decodable*, which is not the same thing and was the gap --
    autocorrelation is powered by repetition, so a short redundant word lights
    it up and a long word carrying information does not. Measured on injected
    codes (research/injection-sensitivity.md): a 21-bit word at 3.3e-3 of a
    cycle reached this gate 5% of the time through the autocorrelation routes
    while its word decoded exactly, and matched a catalog sequence, in over
    half the injections.

    The decode route needs both halves. `p_decode_shuffle` says the two-level
    word is locked to the cycle count better than a reshuffle of the same
    values ever manages -- the draws search word length exactly as the
    observation does, so the choice of L is inside the p-value. `seq_square` is
    stage 2's own `is_square_code`, and it is what keeps smooth astrophysics
    out: a shuffle null is white, so Blazhko modulation and a light-travel-time
    sinusoid also beat it, and only the two-level fit separates them.

    The catalog route needs both halves too. `seq_p` is already a
    look-elsewhere correction over every packed stream, rotation, offset and
    length the matcher tries, and the gate is the one `score.sequence_bonus`
    pays out on; but a smooth wave can correlate with a catalog window by
    chance, so the window has to beat a sinusoid on the steps as well
    (`seq_window_square`, sequence.is_square_match). The test is asked of the
    catalog window rather than the fold's own word because a fixed pattern
    earns it at the two repeats the matcher needs, where `is_square_code`
    needs three: three of the four stars that cleared the bare gate sat at
    2.5 repeats with no squareness at all, two of them with a coherent
    secondary period. The cost: 56 of 1,165 injection recoveries reached
    triage on the catalog route alone, all with the word repeated two to
    three times, which the old `seq_square` guard would have refused outright;
    a clean two-repeat 21-bit counting staircase reads 0.998 on the window
    test against the 0.25 floor. 2 of 320 quiet control runs decoded to any
    catalog word at all and none reached the gate.
    """
    if res.get("error"):
        return False
    v = str(res.get("verdict") or "").lower()
    if any(frag in v for frag in FOLLOWUP_VERDICTS):
        return True
    draws = res.get("shuffle_draws")
    return (
        below_gate(res.get("p_ljung_box"), draws, ACF_P_FOLLOWUP, "Ljung-Box")
        or below_gate(res.get("p_max_acf_shuffle"), draws, ACF_P_FOLLOWUP, "ACF shuffle")
        or (
            below_gate(
                res.get("p_jitter_bootstrap"),
                res.get("jitter_nboot"),
                ACF_P_FOLLOWUP,
                "jitter bootstrap",
            )
            and (res.get("jitter_var_t_over_err") or 0) > JITTER_MIN_SNR
        )
        or (
            below_gate(res.get("p_decode_shuffle"), draws, ACF_P_FOLLOWUP, "decode MDL")
            and bool(res.get("seq_square"))
        )
        or (
            below_gate(res.get("p_quantised_shuffle"), draws, ACF_P_FOLLOWUP, "quantised steps")
            and bool(res.get("seq_square"))
        )
        # Analytic, so no draw count to check against; the Bonferroni over
        # the windows tried is inside the p-value already.
        or (catalog_p(res) < score.SEQ_P_GATE and bool(res.get("seq_window_square")))
        # The kicks: a one-sided set of rare large steps whose spacing reads
        # as a catalog sequence, or favours prime separations, against random
        # kick sets on the same cycles (kick_word).
        or bool(kick_flags(res))
        # The per-cycle fold of a sparse survey: the same two decode gates,
        # asked of the series the blocks averaged away (analyse_epoch_fold).
        or bool(epoch_flags(res))
    )


def catalog_p(res: dict) -> float:
    """The catalog-match p-value the triage gate reads.

    `seq_p` is the viewer's and the score bonus's copy of the matched-filter
    statistic; `p_catalog_matched` is stage 1's own, from its fold of the same
    steps. They normally agree, and a payload written before the sequence
    fields were stored beside it carries only the second.
    """
    for key in ("seq_p", "p_catalog_matched"):
        v = res.get(key)
        if v is not None:
            return float(v)
    return 1.0


def timings_from_frame(oc: pd.DataFrame) -> list[Timing]:
    """The stored per-cycle O-C rows back into Timing objects."""
    out = []
    for r in oc.itertuples(index=False):
        amp = getattr(r, "amp_ratio", 1.0)
        amp = 1.0 if amp is None else amp
        out.append(
            Timing(
                int(r.cycle_index),
                float(r.t_center),
                float(r.shift_days),
                float(r.shift_err_days),
                int(r.npts) if getattr(r, "npts", None) is not None else 0,
                float(r.chi2_red) if getattr(r, "chi2_red", None) is not None else 1.0,
                1.0 if amp is None or not np.isfinite(float(amp)) else float(amp),
            )
        )
    return out


def recompute_routes(
    rec: dict, oc: pd.DataFrame, nshuffle: int, epoch_oc: pd.DataFrame | None = None
) -> tuple[dict, np.ndarray] | None:
    """Re-run analyse_oc on a stored O-C series. No light curve, no network.

    Everything analyse_oc measures depends only on the O-C values and cycle
    indices, which the store has, so a change to a route -- the fold-first
    decode, the quantisation test -- reaches every star with a re-score
    rather than a re-download. The sparse jitter test is the exception (it
    needs the epochs), so its fields and its verdict fragment are kept from
    the stored record.

    The returned record carries `analyse_oc`'s `sequence_row`, which belongs to
    its own table rather than to the payload; `_recompute_one` takes it out and
    `store.payload_fields` is the backstop if a caller does not.
    """
    if oc is None or oc.empty or not rec.get("period_days"):
        return None
    timings = timings_from_frame(oc)
    period = float(rec["period_days"])
    K = int(rec.get("cycles_per_block") or 1)
    ground = not is_dense_survey(rec.get("survey"))
    stats_d, oc_d = analyse_oc(timings, K * period, period, ground, nshuffle=nshuffle)
    if stats_d.get("verdict") is None:
        return None
    old = str(rec.get("verdict") or "")
    jitter = [frag for frag in old.split("; ") if frag.startswith("excess cycle-to-cycle")]
    stats_d["verdict"] = _add_flags(stats_d["verdict"], jitter[:1])
    # The per-cycle series is stored too (stage1_epoch_oc), so the epoch fold
    # is re-read here as well; a row written before it exists keeps whatever
    # epoch fields it had, which is none.
    if epoch_oc is not None and not epoch_oc.empty:
        ep, _ep_d = analyse_epoch_fold(timings_from_frame(epoch_oc), nshuffle=nshuffle)
        stats_d.update(ep)
        stats_d["verdict"] = _add_flags(stats_d["verdict"], epoch_flags(ep))
    return {**rec, **stats_d}, oc_d


def _recompute_one(args: tuple) -> tuple:
    """Pool worker: (target, survey, rec, oc frame, nshuffle)
    -> (target, survey, (rec, [(cycle, detrended days)], sequence row)
        | None | ("error", msg)).

    Everything the store's write needs leaves here in the one result, because
    that write is the single serial thread of a `--recompute-routes` pass and
    whatever is left in it is paid once per star with every core waiting. So
    the sequence row comes out of `analyse_oc` rather than being decoded again
    in the parent -- a pure-Python fold of every word length against every bit
    period over every cycle, which was both the longest step of the write and a
    long enough hold on the GIL to stop the pool's queue feeder refilling the
    workers -- and the detrended O-C leaves as (cycle, days) pairs, so the
    parent does not have to read the frame back out of SQLite to recover the
    cycle indices it already sent us.
    """
    target, survey, rec, oc, nshuffle = args[:5]
    epoch_oc = args[5] if len(args) > 5 else None
    try:
        out = recompute_routes(rec, oc, nshuffle, epoch_oc)
        if out is None:
            return target, survey, None
        new, oc_d = out
        seq = new.pop("sequence_row")
        det = [(int(c), float(d)) for c, d in zip(oc["cycle_index"], oc_d, strict=True)]
        return target, survey, (new, det, seq)
    except Exception as ex:  # a broken row is a skipped row, not a dead run
        return target, survey, ("error", str(ex))


def triage_records(records: Sequence[dict]) -> list[dict]:
    hits = [r for r in records if is_followup(r)]
    hits.sort(
        key=lambda r: (
            min(
                _p(r.get("p_ljung_box")),
                _p(r.get("p_max_acf_shuffle")),
                _p(r.get("p_jitter_bootstrap")),
            ),
            str(r.get("target") or ""),
        )
    )
    return hits


def unpack_result(packed: Any) -> tuple[dict, pd.DataFrame | None]:
    if isinstance(packed, tuple) and len(packed) == 2:
        return packed
    return packed, None


def _period_hint(
    hints: dict[tuple[str, str], float] | None, name: str, survey: str
) -> float | None:
    if not hints:
        return None
    p = hints.get((name, survey))
    if p is None or not np.isfinite(p) or p <= 0:
        return None
    return float(p)


def _drop_skipped(
    names: Sequence[str], surveys: Sequence[str], skip_keys: set[tuple[str, str]]
) -> list[str]:
    """Keep names that still have at least one survey not in skip_keys."""
    if not skip_keys:
        return list(names)
    return [n for n in names if any((n, s) not in skip_keys for s in surveys)]


def run_targets(
    targets: Sequence[str],
    surveys: Sequence[str],
    cadence: str,
    band: str,
    nshuffle: int,
    jobs: int,
    file_index: dict[str, list[str]] | None = None,
    star_class: StarClass | str = DEFAULT_STAR_CLASS,
    store: Store | None = None,
    period_hints: dict[tuple[str, str], float] | None = None,
    fetch_hints: dict[str, dict] | None = None,
    skip_keys: set[tuple[str, str]] | None = None,
    cache: lc_cache.CachePolicy | None = None,
) -> list[dict]:
    jobs = max(1, jobs)
    cls = resolve_cls(star_class)
    file_index = file_index or {}
    skip_keys = skip_keys or set()
    # De-duplicate: a name repeated in a target file used to be downloaded and
    # analysed once per occurrence, and appear that many times in the returned
    # list, while the store kept only the last write.
    keys, seen = [], set()
    n_dup = n_skip = 0
    for n in targets:
        for s in surveys:
            key = (n, s)
            if key in seen:
                n_dup += 1
                continue
            seen.add(key)
            if key in skip_keys:
                n_skip += 1
                continue
            keys.append(key)
    if n_dup:
        log(f"dropped {n_dup} duplicate target/survey pair(s)")
    if n_skip:
        log(f"skip-existing: {n_skip} / {len(seen)} already processed")
    if not keys:
        log("Nothing left to process")
        return []
    by_key = {}
    n_hint = sum(1 for n, s in keys if _period_hint(period_hints, n, s) is not None)
    log(
        f"Processing {len(keys)} target/survey pairs with {jobs} worker(s) [star class: {cls.name}]"
        + (f"; Gaia P for {n_hint}" if n_hint else "")
    )
    pbar = tqdm(total=len(keys), desc="analyze", unit="lc")

    def _finish(key: tuple[str, str], packed: Any) -> None:
        res, oc_df = unpack_result(packed)
        by_key[key] = res
        print_result(res)
        if store is not None:
            store.upsert_stage1(res, oc_df, is_followup(res))
        pbar.update(1)

    try:
        if jobs == 1:
            for n, s in keys:
                _finish(
                    (n, s),
                    process_target(
                        n,
                        s,
                        file_index.get(n, []),
                        cadence,
                        band,
                        nshuffle,
                        cls,
                        _period_hint(period_hints, n, s),
                        cache,
                        (fetch_hints or {}).get(n),
                    ),
                )
        else:
            with worker_pool(jobs) as pool:
                futs = {
                    pool.submit(
                        process_target,
                        n,
                        s,
                        file_index.get(n, []),
                        cadence,
                        band,
                        nshuffle,
                        cls,
                        _period_hint(period_hints, n, s),
                        cache,
                        (fetch_hints or {}).get(n),
                    ): (n, s)
                    for (n, s) in keys
                }
                for fut in as_completed(futs):
                    key = futs[fut]
                    try:
                        packed = fut.result()
                    except Exception as ex:
                        packed = (
                            {"target": key[0], "survey": key[1], "error": f"worker crashed: {ex}"},
                            None,
                        )
                    _finish(key, packed)
    finally:
        pbar.close()
    return [by_key[k] for k in keys]


def _run_catalog_or_names(
    args: argparse.Namespace, cls: StarClass, surveys: Sequence[str], store: Store
) -> list[dict]:
    """Explicit names/file, else every matching stage0_targets row, else the demo list."""
    named = list(args.targets)
    if args.target_file:
        named += [l.strip() for l in open(args.target_file) if l.strip() and not l.startswith("#")]

    skip_keys = (
        store.stage1_skip_keys(retry_errors=args.retry_errors) if args.skip_existing else set()
    )

    def _go(
        names: Sequence[str],
        spec: StarClass,
        use_surveys: Sequence[str],
        file_index: dict[str, list[str]],
        period_hints: dict[tuple[str, str], float] | None = None,
    ) -> list[dict]:
        log(
            f"Star class: {spec.name} (P {spec.pmin}-{spec.pmax} d, <= {spec.max_nharm} harmonics, "
            f"prewhiten {spec.prewhiten})"
        )
        if period_hints is None:
            period_hints = store.period_hints(names, use_surveys)
        # Free here, and what lets a ZTF star skip the cone search entirely.
        fetch_hints = store.ztf_fetch_hints(list(names)) if "ztf" in use_surveys else {}
        return run_targets(
            names,
            use_surveys,
            args.cadence,
            args.band,
            args.nshuffle,
            args.jobs,
            file_index,
            star_class=spec,
            store=store,
            period_hints=period_hints,
            fetch_hints=fetch_hints,
            skip_keys=skip_keys,
            cache=cache_policy(args),
        )

    if named:
        named = _drop_skipped(named, surveys, skip_keys)
        if not named:
            log("skip-existing: every named target already has a stage-1 row")
            return []
        file_index = {}
        tess_named = _drop_skipped(named, ["tess"], skip_keys) if "tess" in surveys else []
        if tess_named and not args.no_bulk:
            file_index = prefetch_lightcurves(
                tess_named, args.cadence, Path(args.download_dir), args.download_jobs
            )
        return _go(named, cls, surveys, file_index)

    want_class = None if args.star_class == "auto" else cls.name
    catalog = store.list_stage0_targets(surveys, want_class)
    if catalog.empty:
        log(
            f"No stage0_targets for {surveys}"
            + (f" / {want_class}" if want_class else "")
            + f"; using the built-in {cls.name} list"
        )
        names = _drop_skipped(DEFAULT_TARGETS[cls.name], surveys, skip_keys)
        if not names:
            log("skip-existing: every demo target already has a stage-1 row")
            return []
        file_index = {}
        tess_named = _drop_skipped(names, ["tess"], skip_keys) if "tess" in surveys else []
        if tess_named and not args.no_bulk:
            file_index = prefetch_lightcurves(
                tess_named, args.cadence, Path(args.download_dir), args.download_jobs
            )
        return _go(names, cls, surveys, file_index)

    classes = [c for c in STAR_CLASSES if c in set(catalog["star_class"])]
    log(f"Store catalog: {len(catalog)} targets across {classes} / {surveys}")

    tess_names = catalog.loc[catalog["survey"] == "tess", "target"].astype(str).tolist()
    tess_names = _drop_skipped(tess_names, ["tess"], skip_keys)
    file_index = {}
    if tess_names and not args.no_bulk:
        file_index = prefetch_lightcurves(
            tess_names, args.cadence, Path(args.download_dir), args.download_jobs
        )

    recs = []
    for class_name in classes:
        sub = catalog[catalog["star_class"] == class_name]
        for survey in surveys:
            chunk = sub.loc[sub["survey"] == survey, ["target", "period_days"]]
            names = _drop_skipped(chunk["target"].astype(str).tolist(), [survey], skip_keys)
            if names:
                hints = {}
                by_name = dict(zip(chunk["target"].astype(str), chunk["period_days"].to_numpy()))
                for tgt in names:
                    p = by_name.get(tgt)
                    if p is not None and p == p and p > 0:
                        hints[(tgt, survey)] = float(p)
                spec = with_prewhiten(STAR_CLASSES[class_name], args.prewhiten)
                recs.extend(_go(names, spec, [survey], file_index, hints))
    return recs


def cache_policy(args: argparse.Namespace) -> lc_cache.CachePolicy:
    """The ASAS-SN / ZTF cache settings for this run, as one picklable value.

    `--lc-cache-age 0` means "never expire", which is the right spelling for a
    reprocessing campaign that wants the exact photometry an earlier pass saw.
    """
    age = getattr(args, "lc_cache_age", lc_cache.DEFAULT_MAX_AGE_DAYS)
    return lc_cache.CachePolicy(
        enabled=not getattr(args, "no_lc_cache", False),
        max_age_days=None if age is not None and age <= 0 else age,
        refresh=bool(getattr(args, "lc_cache_refresh", False)),
    )


def _human_bytes(n: float) -> str:
    """A size a reader can compare at a glance, rather than 0.00 GB."""
    for unit in ("B", "kB", "MB", "GB"):
        if abs(n) < 1000.0 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.2f} {unit}"
        n /= 1000.0
    return f"{n:.2f} GB"


def _print_cache_stats() -> int:
    rows = lc_cache.stats()
    if not rows:
        log(f"ASAS-SN / ZTF cache is empty ({layout.show(lc_cache.root())})")
        return 0
    log(f"ASAS-SN / ZTF raw-photometry cache in {layout.show(lc_cache.root())}")
    for name in sorted(rows):
        r = rows[name]
        oldest = time.strftime("%Y-%m-%d", time.localtime(r["oldest"])) if r["oldest"] else "-"
        newest = time.strftime("%Y-%m-%d", time.localtime(r["newest"])) if r["newest"] else "-"
        log(
            f"  {name:7s} {int(r['targets']):>7,d} targets  "
            f"{int(r['entries']):>7,d} entries  {int(r['rows']):>12,d} epochs  "
            f"{_human_bytes(r['bytes']):>9s}  fetched {oldest} .. {newest}"
        )
    return 0


def _prefetch_targets(
    args: argparse.Namespace, cls: StarClass, survey: str, store: Store
) -> list[str]:
    """Which names `--lc-cache-prefetch` should warm, chosen the way a run does.

    Explicit `--targets` / `--target-file` first, then the stage-0 catalog for
    this survey and class, then the built-in demo list -- the same order
    `_run_catalog_or_names` uses, so prefetching and running select the same
    stars rather than two overlapping sets.
    """
    named = list(args.targets)
    if args.target_file:
        with open(args.target_file, encoding="utf-8") as fh:
            named += [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    if named:
        return named
    want_class = None if args.star_class == "auto" else cls.name
    catalog = store.list_stage0_targets([survey], want_class)
    if not catalog.empty:
        return catalog["target"].astype(str).tolist()
    log(f"No stage0_targets for {survey}; using the built-in {cls.name} list")
    return list(DEFAULT_TARGETS[cls.name])


def _retry_download_errors(args: argparse.Namespace, store: Store) -> set[tuple[str, str]]:
    """Re-run stars whose last attempt died in download or analysis."""
    failed = store.stage1_retry_rows()
    if failed.empty:
        log("No download/analysis failures in the store")
        return []
    file_index = {}
    idx_path = Path(args.download_dir) / "index.json"
    if idx_path.exists():
        try:
            file_index = json.loads(idx_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            file_index = {}
    log(f"Retrying {len(failed)} download/analysis failure(s)")
    recs = []
    for (class_name, survey), chunk in failed.groupby(["star_class", "survey"], sort=False):
        names = chunk["target"].astype(str).tolist()
        spec = with_prewhiten(resolve_cls(class_name), args.prewhiten)
        log(
            f"Star class: {spec.name} (P {spec.pmin}-{spec.pmax} d, <= {spec.max_nharm} harmonics, "
            f"prewhiten {spec.prewhiten}); {len(names)} {survey} target(s)"
        )
        recs.extend(
            run_targets(
                names,
                [survey],
                args.cadence,
                args.band,
                args.nshuffle,
                args.jobs,
                file_index,
                star_class=spec,
                store=store,
                period_hints=store.period_hints(names, [survey]),
                fetch_hints=store.ztf_fetch_hints(list(names)) if survey == "ztf" else {},
                cache=cache_policy(args),
            )
        )
    return recs


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--survey",
        choices=["tess", "kepler", "asassn", "ztf", "both", "all"],
        default="all",
        help="'all' = tess+asassn+ztf (default); 'both' = tess+asassn (legacy); "
        "'kepler' is for named KIC targets only (no stage-0 catalog)",
    )
    ap.add_argument("--targets", nargs="*", default=[])
    ap.add_argument("--target-file")
    ap.add_argument(
        "--cadence",
        default="short",
        help="TESS: short (2-min) or fast (20-s; recommended for bcep where available)",
    )
    ap.add_argument(
        "--star-class",
        default="auto",
        choices=["auto", *STAR_CLASSES],
        help="analysis presets; 'auto' reads the '# star_class:' header of --target-file, "
        "else every class in the store (demo names if the catalog is empty)",
    )
    ap.add_argument(
        "--prewhiten",
        type=int,
        default=None,
        help="override the number of independent modes removed before timing (class default: "
        + ", ".join(f"{k}={v.prewhiten}" for k, v in STAR_CLASSES.items())
        + ")",
    )
    ap.add_argument(
        "--band",
        default=AUTO_BAND,
        choices=[AUTO_BAND, "g", "r", "V", "gr", "all"],
        help="'auto' (default) takes the survey's own default: ASAS-SN g "
        "(2018+; V is 2012-18), ZTF gr -- both filters as one series, with "
        "per-oid amplitude scales. Naming a band applies it to every survey.",
    )
    ap.add_argument("--nshuffle", type=int, default=2000)
    ap.add_argument(
        "--jobs",
        type=int,
        default=max(1, os.cpu_count() or 4),
        help="analysis workers (default: every core). Each one runs single-threaded "
        "numpy, so one per core is the right number rather than an oversubscription",
    )
    ap.add_argument("--download-dir", default=str(layout.LC_CACHE_DIR))
    ap.add_argument("--download-jobs", type=int, default=8)
    ap.add_argument("--no-bulk", action="store_true", help="skip bulk MAST prefetch")
    ap.add_argument(
        "--no-lc-cache",
        action="store_true",
        help="do not read or write the ASAS-SN / ZTF raw-photometry cache",
    )
    ap.add_argument(
        "--lc-cache-refresh",
        action="store_true",
        help="re-download ASAS-SN / ZTF even where a fresh entry exists, and rewrite it",
    )
    ap.add_argument(
        "--lc-cache-age",
        type=float,
        default=lc_cache.DEFAULT_MAX_AGE_DAYS,
        help=f"refuse a cached ASAS-SN / ZTF entry older than this many days "
        f"(default {lc_cache.DEFAULT_MAX_AGE_DAYS:g}; 0 = never expire). Both surveys "
        "keep appending epochs to the same object, so an entry is a snapshot",
    )
    ap.add_argument(
        "--lc-cache-prefetch",
        action="store_true",
        help="download ASAS-SN / ZTF photometry into the cache for the selected "
        "targets and exit, analysing nothing; what a sparse injection grid needs first",
    )
    ap.add_argument(
        "--lc-cache-stats",
        action="store_true",
        help="print what the ASAS-SN / ZTF cache holds and exit",
    )
    ap.add_argument(
        "--lc-cache-prune",
        action="store_true",
        help="delete aged-out and orphaned ASAS-SN / ZTF cache entries and exit; "
        "honours --lc-cache-age",
    )
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--inject", action="store_true")
    ap.add_argument("--triage", action="store_true", help="recompute triage flags in the store")
    ap.add_argument(
        "--recompute-routes",
        action="store_true",
        help="re-run every O-C statistic (ACF, Ljung-Box, fold-first decode, quantised "
        "steps, secondary periods) from the stored O-C series of every star -- no light "
        "curves, no network -- then rewrite the payload, the sequence row and the triage "
        "flag. Run after a change to a route; combine with --triage to rescore as well",
    )
    ap.add_argument(
        "--reseq",
        action="store_true",
        help="re-decode every stored O-C and rewrite its sequence row; run this "
        "after a change to the sequence scoring, before --triage",
    )
    ap.add_argument(
        "--retry-errors",
        action="store_true",
        help="re-run stars whose stage-1 row is 'download failed' or 'analysis failed'; "
        "with --skip-existing, re-run any error row",
    )
    ap.add_argument(
        "--skip-existing",
        action="store_true",
        help="skip (target, survey) pairs that already have a stage-1 row; "
        "combine with --retry-errors to still reprocess error rows",
    )
    ap.add_argument("--verbose", action="store_true", help="print a per-star result line")
    add_db_arg(ap)
    args = ap.parse_args()
    global VERBOSE
    VERBOSE = args.verbose
    # Resolved once, so the star-class header and the target lines are read
    # from the same file whether the name was a path or a bare list name.
    if args.target_file:
        args.target_file = layout.starlist(args.target_file)
    if args.nshuffle < MIN_SHUFFLE:
        log(
            f"--nshuffle {args.nshuffle} cannot resolve p < {ACF_P_FOLLOWUP:g}; "
            f"raising to {MIN_SHUFFLE}"
        )
        args.nshuffle = MIN_SHUFFLE

    surveys = {"both": ["tess", "asassn"], "all": ["tess", "asassn", "ztf"]}.get(
        args.survey, [args.survey]
    )

    # Cache maintenance touches no store and downloads nothing, so it runs
    # before anything opens the database.
    if args.lc_cache_stats:
        return _print_cache_stats()

    if args.lc_cache_prune:
        pol = cache_policy(args)
        removed = lc_cache.prune(max_age_days=pol.max_age_days)
        log(
            f"Pruned {len(removed)} file(s) from {layout.show(lc_cache.root())}"
            + (f" (older than {pol.max_age_days:g} d, plus orphans)" if pol.max_age_days else "")
        )
        return _print_cache_stats()

    if args.lc_cache_prefetch:
        sparse = [s for s in surveys if s in lc_cache.SURVEYS]
        if not sparse:
            log(f"--lc-cache-prefetch needs a sparse survey; got {surveys}")
            return 2
        cls = with_prewhiten(resolve_star_class(args.star_class, args.target_file), args.prewhiten)
        pol = cache_policy(args)
        with Store(args.db) as store:
            for survey in sparse:
                names = _prefetch_targets(args, cls, survey, store)
                if not names:
                    log(f"no {survey} targets selected")
                    continue
                hints = store.ztf_fetch_hints(list(names)) if survey == "ztf" else {}
                prefetch_sparse(names, survey, args.download_jobs, pol, hints)
        return _print_cache_stats()

    if args.recompute_routes:
        with Store(args.db) as store:
            total = store.count_recompute_routes()
            with tqdm(total=total, desc="routes", unit="star") as bar:
                n = store.recompute_routes(
                    _recompute_one,
                    is_followup,
                    args.nshuffle,
                    jobs=args.jobs,
                    progress=lambda _i: bar.update(1),
                )
            print(f"Recomputed the routes for {n} / {total} O-C series in {store.path}")
        if not args.triage:
            return 0

    if args.reseq:
        with Store(args.db) as store:
            total = store.count_stage1()
            with tqdm(total=total, desc="re-decode", unit="star") as bar:
                n = store.recompute_sequence(progress=lambda _i: bar.update(1))
            print(f"Re-decoded {n} / {total} O-C series in the store ({store.path})")
        if not args.triage:
            return 0

    if args.triage:
        with Store(args.db) as store:
            n = store.recompute_triage(is_followup)
            total = store.count_stage1()
            print(f"Marked {n} / {total} follow-ups in the store ({store.path})")
        return 0

    if args.retry_errors and not args.skip_existing:
        with Store(args.db) as store:
            recs = _retry_download_errors(args, store)
            n_triage = sum(1 for r in recs if is_followup(r))
            print(
                f"\nRetried {len(recs)} download/analysis failures "
                f"({n_triage} triage) in {store.path}"
            )
        return 0

    cls = with_prewhiten(resolve_star_class(args.star_class, args.target_file), args.prewhiten)

    with Store(args.db) as store:
        if args.synthetic:
            log(
                f"Star class: {cls.name} (P {cls.pmin}-{cls.pmax} d, <= {cls.max_nharm} harmonics, "
                f"prewhiten {cls.prewhiten})"
            )
            name = f"SYNTH {cls.name}" + ("" if args.inject else " null")
            skip_keys = (
                store.stage1_skip_keys(retry_errors=args.retry_errors)
                if args.skip_existing
                else set()
            )
            recs = run_targets(
                [name],
                surveys,
                args.cadence,
                args.band,
                args.nshuffle,
                1,
                star_class=cls,
                store=store,
                skip_keys=skip_keys,
            )
            if args.inject and recs:
                missed = [
                    f"{r.get('target')} [{r.get('survey')}]" for r in recs if not is_followup(r)
                ]
                if missed:
                    print(
                        "\nFAILED: injected control did not reach triage: " + ", ".join(missed),
                        file=sys.stderr,
                    )
                    print(
                        "  This is a positive control. If it does not fire, the "
                        "detection path is not being exercised at all.",
                        file=sys.stderr,
                    )
                    return 1
                log(f"injected control reached triage in all {len(recs)} run(s)")
        else:
            recs = _run_catalog_or_names(args, cls, surveys, store)
        n_triage = sum(1 for r in recs if is_followup(r))
        print(f"\nStored {len(recs)} summaries ({n_triage} triage) in {store.path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
