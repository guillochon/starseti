#!/usr/bin/env python3
"""Figure: the three built-in synthetic controls, timed by the production
stage-1 code, showing what a phase-coded pulsator would look like.

Each row is one star class (classical Cepheid, RR Lyrae, beta Cephei) with the
pipeline's own injected word (primes, counting, Fibonacci) at dt/P = 1e-3 per set
bit, sampled at the TESS 2-min cadence with sector gaps.  Left: a stretch of the
light curve.  Middle: the per-cycle O-C from the template time-shift fit, with the
injected staircase (a set bit resets the clock and every later cycle inherits the
advance) drawn through it.  Right: the first difference of the O-C -- the steps --
which is the series the decoder reads; the injected bit of each cycle colours the
point and the word the decoder recovered is printed.

    uv run python paper/figures/make_synthetic.py      # writes synthetic_examples.pdf

Runs the real `stage1_pulse.synthetic` -> `run_target` path; the Cepheid row takes
a few minutes because 2-min sampling puts thousands of points in each 5.4 d cycle.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import sequence  # noqa: E402
import stage1_pulse as s1  # noqa: E402

OUT = Path(__file__).with_name("synthetic_examples.pdf")

ROWS = [
    ("cepheid", "Classical Cepheid, $P=5.37$ d", "primes 2, 3, 5, ..., 19"),
    ("rrlyrae", "RR Lyrae (RRab), $P=0.567$ d", "counting 1, 2, ..., 8"),
    ("bcep", r"$\beta$ Cephei, $P=0.201$ d", "Fibonacci 1, 1, 2, ..., 21"),
]
BIT0, BIT1 = "#3b6ea5", "#d94f3d"


def detrend(x: np.ndarray, y: np.ndarray, deg: int = 2) -> np.ndarray:
    ok = np.isfinite(y)
    if ok.sum() < deg + 2:
        return y
    c = np.polyfit(x[ok], y[ok], deg)
    return y - np.polyval(c, x)


def one_row(star_class: str, nshuffle: int = 200):
    lcd = s1.synthetic("tess", inject=True, star_class=star_class)
    cls = s1.STAR_CLASSES[star_class]
    res, df = s1.run_target(f"SYNTH {star_class}", lcd, nshuffle, cls)
    if df is None:
        raise RuntimeError(f"{star_class}: {res.get('note')}")
    period = res["period_days"]
    # the injected code, as the light curve carried it
    cyc_all = np.floor((lcd.t - lcd.t.min()) / period).astype(int)
    code = s1.synth_inject_code(cyc_all, star_class)
    cyc = df["cycle_index"].to_numpy()
    # the stage-1 timer numbers cycles from its own t0; align by nearest centre
    t_c = df["t_center"].to_numpy()
    inj_cyc = np.floor((t_c - lcd.t.min()) / period).astype(int).clip(0, len(code) - 1)
    inj_shift_days = s1.code_shift(code, inj_cyc, s1.SYNTH_INJECT_FRAC, period)
    meas = df["shift_days"].to_numpy() * 86400.0
    inj = inj_shift_days * 86400.0
    steps_meas = np.full_like(meas, np.nan)
    steps_bit = np.full(meas.shape, -1)
    for i in range(1, len(cyc)):
        if cyc[i] - cyc[i - 1] == 1:
            steps_meas[i] = meas[i] - meas[i - 1]
            steps_bit[i] = code[inj_cyc[i]]
    dec = sequence.decode_bit_code(list(df["oc_detrended_days"] * 86400.0), list(cyc))
    return lcd, res, df, period, meas, inj, steps_meas, steps_bit, dec


def main() -> None:
    plt.rcParams.update({"font.size": 8, "font.family": "serif"})
    fig, axes = plt.subplots(3, 3, figsize=(7.1, 7.2), gridspec_kw={"width_ratios": [1.0, 1.35, 1.35], "hspace": 0.62, "wspace": 0.45})
    handles_steps = ([], [])
    for (star_class, title, word), axrow in zip(ROWS, axes):
        lcd, res, df, period, meas, inj, steps, bits, dec = one_row(star_class)
        cyc = df["cycle_index"].to_numpy()
        err = df["shift_err_days"].to_numpy() * 86400.0
        # --- light curve excerpt: ~3 cycles from the first sector
        ax = axrow[0]
        t0 = lcd.t.min()
        m = (lcd.t - t0) < 3.0 * period
        ax.plot(lcd.t[m] - t0, lcd.f[m], ".", ms=1.2, color="0.25", rasterized=True)
        ax.set_xlabel("time (d)")
        ax.set_ylabel("relative flux")
        ax.set_title(title, loc="left", fontsize=8)
        # --- O-C with the injected staircase
        ax = axrow[1]
        x = cyc.astype(float)
        m_d = detrend(x, meas)
        i_d = detrend(x, inj)
        ax.errorbar(cyc, m_d, yerr=err, fmt="o", ms=2.0, color="0.2", ecolor="0.6", elinewidth=0.6, capsize=0, label="measured per-cycle O$-$C", zorder=3)
        ax.step(cyc, i_d, where="mid", color=BIT1, lw=0.9, alpha=0.8, label="injected staircase (detrended)", zorder=2)
        ax.set_xlabel("pulsation cycle")
        ax.set_ylabel("O$-$C (s)")
        ax.set_title(f"O$-$C, {len(cyc)} timed cycles", loc="left", fontsize=8)
        # --- steps, coloured by the injected bit
        ax = axrow[2]
        ok = np.isfinite(steps)
        for b, col, lab in ((0, BIT0, "bit 0"), (1, BIT1, "bit 1")):
            sel = ok & (bits == b)
            ax.plot(cyc[sel], steps[sel], "o", ms=2.6, color=col, label=lab, zorder=3)
        mid = sequence.two_level_mid([float(v) for v in steps[ok]])
        if mid is not None:
            ax.axhline(mid, color="0.5", lw=0.7, ls="--")
        ax.axhline(0, color="0.8", lw=0.5)
        ax.set_xlabel("pulsation cycle")
        ax.set_ylabel("O$-$C step (s)")
        label = "decoder: no word"
        if dec:
            name = dec.get("seq_label") or dec.get("label") or dec.get("seq_name") or ""
            label = f"decoded: {dec.get('n', '?')}-bit word" + (f"\nnamed: {name}" if name else "")
        ax.set_title(f"steps; injected {word}", loc="left", fontsize=8)
        ax.margins(y=0.30)
        ax.text(0.02, 0.97, label, transform=ax.transAxes, fontsize=7, va="top", ha="left",
                bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "0.7", "lw": 0.5})
        if star_class == "cepheid":
            handles_steps = ax.get_legend_handles_labels()
    h1, l1 = axes[0][1].get_legend_handles_labels()
    h2, l2 = handles_steps
    fig.legend(h1 + h2, l1 + [f"injected {x}" for x in l2], loc="upper center", ncol=4, fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.985))
    fig.subplots_adjust(top=0.93)
    fig.savefig(OUT, bbox_inches="tight", dpi=200)
    fig.savefig(OUT.with_suffix(".png"), bbox_inches="tight", dpi=170)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
