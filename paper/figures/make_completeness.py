#!/usr/bin/env python3
"""Figure: injection-recovery completeness of the stage-1 search on real hosts,
from data/injection.sqlite.

Top row: TESS, the `W-<class>-fold` runs (staircase waveform, one bit per cycle,
40 quiet hosts per class, an 8-bit random word and the 28-bit primes word,
amplitudes scanned in log dt/P). Two gates per class: `flagged` (the star would
have been sent to vetting) and `decoded` (>= 90 per cent of the injected bits
read back correctly).

Bottom row: ASAS-SN and ZTF RR Lyrae, the `R1-asassn-rrlyrae` and
`R2-ztf-rrlyrae` runs (40 quiet hosts each, the same two words, the code held
for 1, 8 and 32 pulsation cycles per bit). Only the flagged gate is drawn: no
ground-based injection was decoded (research/response1.md, M2).

Vertical lines mark the Learned et al. (2008) trigger scaling and the kick
ceiling of a self-powered transmitter without energy storage, per class
(src/physics.py; research/modulation-physics.md; research/response1.md, M1).

    uv run python paper/figures/make_completeness.py     # writes completeness.pdf
"""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import physics  # noqa: E402

DB = ROOT / "data" / "injection.sqlite"
OUT = Path(__file__).with_name("completeness.pdf")
CLASSES = [("rrlyrae", "RR Lyrae"), ("cepheid", "classical Cepheids"), ("bcep", r"$\beta$ Cephei")]
WORDS = [("random:8", "8-bit random word", "#3b6ea5"), ("primes", "28-bit primes word", "#d94f3d")]
SPARSE = [("R1-asassn-rrlyrae", "ASAS-SN RR Lyrae"), ("R2-ztf-rrlyrae", "ZTF RR Lyrae")]
BIT_CYCLES = [(1, ":"), (8, "--"), (32, "-")]


def binned(frac: np.ndarray, hit: np.ndarray, edges: np.ndarray):
    x, y, lo, hi = [], [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (frac >= a) & (frac < b)
        n = int(m.sum())
        if n == 0:
            continue
        k = int(hit[m].sum())
        p = k / n
        # Wilson interval
        z = 1.0
        den = 1 + z**2 / n
        cen = (p + z**2 / (2 * n)) / den
        half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / den
        x.append(np.sqrt(a * b))
        y.append(p)
        lo.append(max(0.0, cen - half))
        hi.append(min(1.0, cen + half))
    return np.array(x), np.array(y), np.array(lo), np.array(hi)


def style(ax, kick: float) -> None:
    ax.axvline(1e-3, color="0.4", lw=0.8, ls=":")
    ax.axvline(kick, color="0.4", lw=0.8, ls="-.")
    ax.set_xscale("log")
    ax.set_xlim(8e-6, 0.06)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(r"injected $\delta t / P$ per set bit")
    ax.grid(alpha=0.2, lw=0.5)


def main() -> None:
    con = sqlite3.connect(f"file:{DB.as_posix()}?mode=ro", uri=True)
    cur = con.cursor()
    plt.rcParams.update({"font.size": 8, "font.family": "serif"})
    fig = plt.figure(figsize=(7.1, 5.2))
    gs = fig.add_gridspec(2, 6, hspace=0.55, wspace=0.55)
    top = [fig.add_subplot(gs[0, 0:2]), None, None]
    top[1] = fig.add_subplot(gs[0, 2:4], sharey=top[0])
    top[2] = fig.add_subplot(gs[0, 4:6], sharey=top[0])
    bottom = [fig.add_subplot(gs[1, 0:3], sharey=top[0]), fig.add_subplot(gs[1, 3:6], sharey=top[0])]
    edges = np.logspace(-5, np.log10(0.05) + 1e-9, 9)
    # Per class: the annular coupling differs by class once the cross-section
    # turns over (review 1, M1), so each panel carries its own kick ceiling.
    kicks = {cls: physics.kick_ceiling(physics.coupling(cls, physics.DEFAULT_GEOMETRY)) for cls, _ in CLASSES}

    for (cls, clabel), ax in zip(CLASSES, top):
        n_hosts = cur.execute(
            "select count(distinct target) from injections where run_id=? and injected=1", (f"W-{cls}-fold",)
        ).fetchone()[0]
        for word, wlabel, col in WORDS:
            rows = cur.execute(
                "select frac, triage, bit_recall, stats_ran, error from injections "
                "where run_id=? and word=? and injected=1",
                (f"W-{cls}-fold", word),
            ).fetchall()
            if not rows:
                continue
            frac = np.array([r[0] for r in rows], float)
            flagged = np.array([1 if r[1] else 0 for r in rows], float)
            decoded = np.array([1 if (r[2] is not None and r[2] >= 0.9) else 0 for r in rows], float)
            for hit, ls, mk, glabel in ((flagged, "-", "o", "flagged"), (decoded, "--", "s", "decoded")):
                x, y, lo, hi = binned(frac, hit, edges)
                ax.plot(x, y, ls=ls, marker=mk, ms=3, color=col, lw=1.1, label=f"{wlabel}, {glabel}")
                ax.fill_between(x, lo, hi, color=col, alpha=0.10, lw=0)
        style(ax, kicks[cls])
        ax.set_title(f"{clabel} ({n_hosts} TESS hosts)", fontsize=8)
    for ax in top[1:]:
        plt.setp(ax.get_yticklabels(), visible=False)
    top[0].set_ylabel("fraction recovered")
    top[0].text(1e-3 / 1.12, 0.93, "Learned et al.\n$10^{-3}$", fontsize=6.5, color="0.3", ha="right", va="top")
    top[0].text(
        kicks[CLASSES[0][0]] * 1.08,
        0.98,
        "$G=1$\nkick ceiling",
        fontsize=6.5,
        color="0.3",
        ha="left",
        va="top",
        linespacing=0.95,
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", pad=1.15),
    )
    h, l = top[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=4, fontsize=6.5, frameon=False, bbox_to_anchor=(0.5, 0.97))

    for (run, slabel), ax in zip(SPARSE, bottom):
        n_hosts = cur.execute(
            "select count(distinct target) from injections where run_id=? and injected=1", (run,)
        ).fetchone()[0]
        for word, wlabel, col in WORDS:
            for k, ls in BIT_CYCLES:
                rows = cur.execute(
                    "select frac, triage from injections where run_id=? and word=? and cycles_per_bit=? "
                    "and injected=1 and stats_ran=1 and error is null",
                    (run, word, k),
                ).fetchall()
                if not rows:
                    continue
                frac = np.array([r[0] for r in rows], float)
                flagged = np.array([1 if r[1] else 0 for r in rows], float)
                x, y, lo, hi = binned(frac, flagged, edges)
                ax.plot(x, y, ls=ls, marker="o", ms=2.5, color=col, lw=1.0)
                if k == 32:
                    ax.fill_between(x, lo, hi, color=col, alpha=0.10, lw=0)
        style(ax, kicks["rrlyrae"])
        ax.set_title(f"{slabel} ({n_hosts} hosts), flagged", fontsize=8)
    plt.setp(bottom[1].get_yticklabels(), visible=False)
    bottom[0].set_ylabel("fraction flagged")
    handles = [Line2D([], [], color="0.25", ls=ls, lw=1.0, label=f"{k} cycle{'s' if k > 1 else ''} per bit") for k, ls in BIT_CYCLES]
    bottom[1].legend(handles=handles, loc="upper left", fontsize=6.5, frameon=False)

    fig.savefig(OUT, bbox_inches="tight")
    fig.savefig(OUT.with_suffix(".png"), bbox_inches="tight", dpi=170)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
