#!/usr/bin/env python3
"""Figure: the best-scoring star beside the star with the most significant
named word outside the counting family, drawn from the stored pipeline
products.

Left column: TIC 152825245, first in the viewer's default score order, with no
named word.  Right column: TIC 147576267, the first TESS star in the viewer's
p-value order once the counting windows are set aside (triangular 561, 595,
minimal binary, inverted).  Rows:

  (a) the phase-folded TESS light curve with the Fourier template;
  (b) the per-cycle O-C from the template fit, with the stage-2 steepest-branch
      crossing timing over it (the independent method) where stage 2 ran, per
      TESS segment;
  (c) Lomb-Scargle periodograms, in cycles^-1, of the per-cycle O-C and of the
      per-cycle amplitude -- a second pulsation mode modulates both at the same
      frequency, a phase-only transmitter would modulate one;
  (d) the steps (first differences of the O-C, winsorized as the decoder reads
      them) averaged per residue of the decoded word length L, with the
      two-level word, the best single sinusoid at any frequency up to the
      once-per-cycle Nyquist limit folded the same way, and, on the right, the
      catalog window the matched filter named, at its best alignment.

    uv run python paper/figures/make_candidates.py     # writes candidates.pdf

Reads data/cepheid.sqlite read-only; needs no light curves or network.
"""

from __future__ import annotations

import json
import sqlite3
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.timeseries import LombScargle
from matplotlib.gridspec import GridSpecFromSubplotSpec

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import sequence  # noqa: E402
from store import fold_from_blobs  # noqa: E402

DB = ROOT / "data" / "cepheid.sqlite"
OUT = Path(__file__).with_name("candidates.pdf")

STARS = [
    ("TIC 152825245", "TIC 152825245, RRc", "top score, no named word"),
    ("TIC 147576267", "TIC 147576267, RRc", "top non-counting named word"),
]
C_OC, C_AMP, C_IND, C_WORD, C_SIN, C_WIN = (
    "#222222",
    "#d94f3d",
    "#3b6ea5",
    "#3b6ea5",
    "#e08a1e",
    "#7a7a7a",
)


def load(con: sqlite3.Connection, target: str) -> dict:
    s1 = json.loads(
        con.execute("SELECT payload FROM stage1_summary WHERE target = ?", (target,)).fetchone()[0]
    )
    row = con.execute("SELECT payload FROM stage2_vetting WHERE target = ?", (target,)).fetchone()
    s2 = json.loads(row[0]) if row else None  # QUIET stars never reach stage 2
    code = json.loads(
        con.execute(
            "SELECT seq_code_json FROM stage1_sequence WHERE target = ?", (target,)
        ).fetchone()[0]
    )
    fold = fold_from_blobs(
        *con.execute(
            "SELECT phase, flux, template_flux FROM stage1_fold WHERE target = ?", (target,)
        ).fetchone()
    )
    oc = np.array(
        con.execute(
            "SELECT cycle_index, t_center, oc_detrended_days, shift_err_days, amp_ratio "
            "FROM stage1_oc WHERE target = ? AND survey = 'tess' ORDER BY cycle_index",
            (target,),
        ).fetchall(),
        float,
    )
    return {"s1": s1, "s2": s2, "code": code, "fold": fold, "oc": oc}


def segments(t: np.ndarray, gap: float = 10.0) -> list[np.ndarray]:
    cut = np.where(np.diff(t) > gap)[0]
    return np.split(np.arange(t.size), cut + 1)


def panel_lc(ax, d: dict, title: str, sub: str) -> None:
    f = d["fold"]
    ph, fl = np.asarray(f["phase"]), np.asarray(f["flux"])
    ax.plot(ph, fl, ".", ms=1.0, color="0.55", rasterized=True)
    ax.plot(f["template_phase"], f["template_flux"], "-", color=C_AMP, lw=1.2)
    ax.set_xlim(0, 1)
    ax.set_xlabel("phase")
    ax.set_ylabel("relative flux")
    P = d["s1"]["period_days"]
    ax.set_title(f"{title}, $P={P:.4f}$ d\n{sub}", fontsize=8.5)


def panel_oc(cell, fig, d: dict) -> None:
    oc = d["oc"]
    cyc, t, y, e = oc[:, 0], oc[:, 1], oc[:, 2] * 86400.0, oc[:, 3] * 86400.0
    cr = ((d["s2"] or {}).get("robustness") or {}).get("crossing") or {"t": [], "oc_s": []}
    tc, yc = np.asarray(cr["t"], float), np.asarray(cr["oc_s"], float)
    segs = segments(t)
    widths = [max(float(np.ptp(t[s])), 1.0) for s in segs]
    sub = GridSpecFromSubplotSpec(1, len(segs), subplot_spec=cell, width_ratios=widths, wspace=0.06)
    lim = np.nanpercentile(np.abs(np.r_[y, yc]), 99.5) * 1.15
    for i, s in enumerate(segs):
        ax = fig.add_subplot(sub[0, i])
        lo, hi = t[s].min() - 0.5, t[s].max() + 0.5
        m = (tc >= lo) & (tc <= hi)
        if tc.size:
            ax.plot(
                tc[m],
                yc[m],
                "o",
                ms=2.6,
                mfc="none",
                mec=C_IND,
                mew=0.6,
                label="steepest-branch crossing",
            )
        ax.errorbar(
            t[s], y[s], e[s], fmt=".", ms=2.5, color=C_OC, elinewidth=0.5, label="template fit"
        )
        ax.set_xlim(lo, hi)
        ax.set_ylim(-lim, lim)
        ax.axhline(0, color="0.8", lw=0.5, zorder=0)
        if i == 0:
            ax.set_ylabel("O$-$C (s)")
        else:
            ax.tick_params(labelleft=False)
            ax.spines["left"].set_linestyle((0, (2, 2)))
        if i < len(segs) - 1:
            ax.spines["right"].set_linestyle((0, (2, 2)))
        if len(segs) > 2:  # narrow sectors: one tick each, at the middle
            ax.set_xticks([10 * round((lo + hi) / 20)])
        else:
            ax.locator_params(axis="x", nbins=3 if len(segs) > 1 else 5)
        if i == len(segs) // 2:
            ax.set_xlabel(f"BTJD (d); {len(cyc)} timed cycles")
        if i == 0:
            ax.legend(
                loc="lower left",
                bbox_to_anchor=(0.0, 1.0),
                ncol=2,
                fontsize=6,
                frameon=False,
                handletextpad=0.2,
                columnspacing=0.8,
                borderaxespad=0.1,
            )
            ax.text(
                -0.02,
                1.02,
                "(b)",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=8,
                weight="bold",
            )


def panel_pgram(ax, d: dict, marks: list[tuple[float, str]]) -> None:
    oc = d["oc"]
    cyc, y, amp = oc[:, 0], oc[:, 2], oc[:, 4]
    f = np.linspace(0.004, 0.5, 6000)
    for series, col, lab in ((y, C_OC, "O$-$C"), (amp, C_AMP, "amplitude")):
        pw = LombScargle(cyc, series).power(f)
        ax.plot(f, pw / pw.max(), "-", color=col, lw=0.9, label=lab)
    for fx, lab in marks:
        ax.axvline(fx, color="0.5", ls=":", lw=0.8)
        right = fx > 0.35
        ax.text(
            fx - 0.006 if right else fx + 0.006,
            1.02,
            lab,
            fontsize=6,
            va="bottom",
            color="0.3",
            ha="right" if right else "left",
        )
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, 1.18)
    ax.set_xlabel("frequency (cycle$^{-1}$)")
    ax.set_ylabel("normalized power")
    ax.legend(loc="upper right", fontsize=6.5, frameon=False, bbox_to_anchor=(1.0, 0.9))


def fold_steps(d: dict):
    """The decoder's own step series, folded at its word length."""
    oc = d["oc"]
    pts = sequence.code_points(list(oc[:, 2] * 86400.0), [int(c) for c in oc[:, 0]], "staircase")
    c = np.array([p["c"] for p in pts], int)
    ys = np.array([p["y"] for p in pts], float)
    L, k = int(d["code"]["n"]), int(d["code"]["k"])
    res = ((c // k) % L + L) % L
    means = np.array([ys[res == r].mean() for r in range(L)])
    sem = np.array(
        [
            ys[res == r].std(ddof=1) / np.sqrt((res == r).sum()) if (res == r).sum() > 1 else np.nan
            for r in range(L)
        ]
    )
    counts = np.array([(res == r).sum() for r in range(L)])
    return c, ys, res, means, sem, counts, L


def two_level(means: np.ndarray, counts: np.ndarray, bits: np.ndarray) -> np.ndarray:
    lv = {b: np.average(means[bits == b], weights=counts[bits == b]) for b in (0, 1)}
    return np.where(bits == 1, lv[1], lv[0])


def r2(model: np.ndarray, means: np.ndarray, counts: np.ndarray) -> float:
    mu = np.average(means, weights=counts)
    return 1.0 - np.sum(counts * (means - model) ** 2) / np.sum(counts * (means - mu) ** 2)


def panel_word(ax, d: dict, window: str | None, window_label: str = "") -> dict:
    c, ys, res, means, sem, counts, L = fold_steps(d)
    # the fold's own word: the best two-cluster split of the residue means
    order = np.sort(means)
    cost = [np.var(order[:i]) * i + np.var(order[i:]) * (L - i) for i in range(1, L)]
    thr = 0.5 * (order[int(np.argmin(cost))] + order[int(np.argmin(cost)) + 1])
    bits = (means > thr).astype(int)
    word_model = two_level(means, counts, bits)
    # the best single sinusoid over the unfolded steps, folded the same way
    rss_f = sequence._best_sinusoid_rss(list(c), list(ys))
    nu = rss_f[1]
    ang = 2 * np.pi * nu * c
    basis = np.column_stack([np.ones(c.size), np.cos(ang), np.sin(ang)])
    coef, *_ = np.linalg.lstsq(basis, ys, rcond=None)
    pred = basis @ coef
    sin_model = np.array([pred[res == r].mean() for r in range(L)])
    x = np.arange(L)
    ax.errorbar(
        x,
        means,
        sem,
        fmt="o",
        ms=3,
        color=C_OC,
        elinewidth=0.6,
        zorder=3,
        label="folded step means",
    )
    ax.step(
        np.r_[x, L] - 0.5,
        np.r_[word_model, word_model[-1]],
        where="post",
        color=C_WORD,
        lw=1.1,
        label=f"two-level word ($R^2={r2(word_model, means, counts):.2f}$)",
    )
    ax.plot(
        x,
        sin_model,
        "-",
        color=C_SIN,
        lw=1.1,
        label=f"one sinusoid, {1 / nu:.2f} cycles ($R^2={r2(sin_model, means, counts):.2f}$)",
    )
    out = {
        "nu": nu,
        "r2_word": r2(word_model, means, counts),
        "r2_sin": r2(sin_model, means, counts),
        "counts": (counts.min(), counts.max()),
        "L": L,
    }
    if window:
        w = np.array([int(b) for b in window])
        best = None
        for sgn in (1, -1):
            for r in range(L):
                wr = np.roll(w, r)
                if wr.min() == wr.max():
                    continue
                rho = sgn * np.corrcoef(wr, means)[0, 1]
                if best is None or rho > best[0]:
                    best = (rho, r, sgn)
        rho, r, sgn = best
        wr = np.roll(w, r) if sgn > 0 else 1 - np.roll(w, r)
        win_model = two_level(means, counts, wr)
        miss = np.where(wr != bits)[0]
        ax.step(
            np.r_[x, L] - 0.5,
            np.r_[win_model, win_model[-1]],
            where="post",
            color=C_WIN,
            lw=0.9,
            ls="--",
            label=f"named window, {window_label}",
        )
        for m in miss:
            ax.axvspan(m - 0.5, m + 0.5, color=C_WIN, alpha=0.18, lw=0)
        out.update({"r2_window": r2(win_model, means, counts), "mismatch": len(miss), "rho": rho})
    ax.axhline(0, color="0.8", lw=0.5, zorder=0)
    ax.set_xlim(-0.8, L - 0.2)
    ax.set_xlabel(f"residue (cycle mod {L})")
    ax.set_ylabel("mean O$-$C step (s)")
    lo, hi = np.nanmin(means - np.nan_to_num(sem)), np.nanmax(means + np.nan_to_num(sem))
    span = hi - lo
    ax.set_ylim(lo - 0.05 * span, hi + 0.62 * span)
    ax.legend(loc="upper left", fontsize=6, frameon=False, handlelength=1.6)
    ax.text(
        0.99,
        0.02,
        f"{counts.min()}$-${counts.max()} steps per residue",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6,
    )
    return out


def tag(ax, lab: str) -> None:
    ax.text(
        -0.02,
        1.02,
        f"({lab})",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        weight="bold",
    )


def main() -> None:
    con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    plt.rcParams.update({"font.size": 7.5, "axes.titlesize": 8.5, "font.family": "serif"})
    fig = plt.figure(figsize=(7.1, 7.4))
    gs = fig.add_gridspec(4, 2, hspace=0.72, wspace=0.28, height_ratios=[0.85, 1, 0.85, 1.25])
    marks = {"TIC 152825245": [(0.187, "half the 0.61-mode alias")], "TIC 147576267": []}
    for col, (target, title, sub) in enumerate(STARS):
        d = load(con, target)
        ax = fig.add_subplot(gs[0, col])
        panel_lc(ax, d, title, sub)
        tag(ax, "a")
        panel_oc(gs[1, col], fig, d)
        ax = fig.add_subplot(gs[2, col])
        panel_pgram(ax, d, marks[target])
        tag(ax, "c")
        recog = d["code"]["recog"]
        window = recog.get("window")
        family = recog.get("label", "").split(" [")[0]
        nums = ", ".join(str(n) for n in recog.get("nums", []))
        ax = fig.add_subplot(gs[3, col])
        info = panel_word(ax, d, window, f"{family} {nums}")
        tag(ax, "d")
        print(target, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in info.items()})
    fig.savefig(OUT, bbox_inches="tight", dpi=300)
    fig.savefig(OUT.with_suffix(".png"), bbox_inches="tight", dpi=150)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
