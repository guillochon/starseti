#!/usr/bin/env python3
"""Figure: the transmitter geometry -- host pulsator, driving zone, collector
swarm and the hollow-cone neutrino beam that grazes the limb.

The drawing follows research/modulation-physics.md section 5.  It is a schematic,
not to scale: the swarm sits at tens of stellar radii and the driving zone is a
few per cent of the radius deep, so both are exaggerated to be legible.  Labels
only name the components; what each one does is in the caption.  Every label
sits clear of the drawing with a thin leader to the thing it names, except the
two that sit on the coloured areas they name (pulsator, photosphere), and
`check_labels` refuses to write the figure if any other label touches a drawn
line or shape.  The collector tiles sit exactly on their orbit, and each beam
is fired by one of those tiles.  The chords in panel (b) follow the formulae in
the note (radial chord H per crossing, grazing chord 2 sqrt(2 r H)), and the
grazing chord is drawn whole.

    python paper/figures/make_geometry.py         # writes geometry.pdf next to it
"""

from __future__ import annotations

from itertools import pairwise
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Wedge
from matplotlib.text import Text
from matplotlib.transforms import Bbox

OUT = Path(__file__).with_name("geometry.pdf")

STAR_FILL = plt.cm.YlOrBr(0.25)  # type: ignore[attr-defined]
STAR_EDGE = "#c98a12"
STAR_INK = "#5a3a00"
ZONE = "#d94f3d"
SWARM = "#3b6ea5"
BEAM = "#7b3fa0"
INK = "#222222"
LEADER = {"arrowstyle": "-", "lw": 0.7, "shrinkA": 1, "shrinkB": 2}


def unit(deg: float) -> np.ndarray:
    a = np.radians(deg)
    return np.array([np.cos(a), np.sin(a)])


def panel_a(ax) -> list:
    """The star, its driving zone, the swarm and two firing tiles."""
    R = 1.0
    r_zone = 0.86  # He II partial-ionisation zone, exaggerated inward
    R_swarm = 2.45  # collector orbit, ~43 R_star in reality (broken scale)
    ax.set_aspect("equal")
    ax.set_xlim(-4.1, 2.9)
    ax.set_ylim(-3.4, 3.3)
    ax.axis("off")

    ax.add_patch(Circle((0, 0), R, facecolor=STAR_FILL, edgecolor=STAR_EDGE, lw=1.0, zorder=2))
    ax.add_patch(
        Wedge((0, 0), r_zone + 0.045, 0, 360, width=0.09, facecolor=ZONE, alpha=0.85, zorder=4)
    )
    on_colour = [
        ax.text(0, 0, "pulsator", ha="center", va="center", fontsize=8, color=STAR_INK, zorder=5)
    ]

    # collector swarm: tiles exactly on a dash-dotted orbit
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(
        R_swarm * np.cos(th),
        R_swarm * np.sin(th),
        ls=(0, (4, 2, 1, 2)),
        lw=0.7,
        color=SWARM,
        zorder=1,
    )
    n_tiles = 24
    tile_deg = [360.0 * i / n_tiles for i in range(n_tiles)]
    for deg in tile_deg:
        c = R_swarm * unit(deg)
        tvec = 0.11 * unit(deg + 90)
        ax.plot(
            [c[0] - tvec[0], c[0] + tvec[0]],
            [c[1] - tvec[1], c[1] + tvec[1]],
            color=SWARM,
            lw=2.4,
            zorder=6,
            solid_capstyle="round",
        )

    # two of the tiles fire hollow-cone beams whose edges graze the driving zone
    emitters = [tile_deg[14], tile_deg[19]]  # 210 deg and 285 deg
    for deg in emitters:
        E = R_swarm * unit(deg)
        d = np.linalg.norm(E)
        half = np.arcsin(r_zone / d)  # cone half-angle that grazes the shell
        toward = np.arctan2(-E[1], -E[0])
        L_t = np.sqrt(d**2 - r_zone**2)  # emitter to tangent point
        tangents = []
        for s in (+1, -1):
            u = np.array([np.cos(toward + s * half), np.sin(toward + s * half)])
            tangents.append(E + L_t * u)
            far = E + 2 * L_t * u  # the grazing ray leaves the shell on the far side
            ax.plot([E[0], far[0]], [E[1], far[1]], color=BEAM, lw=1.0, zorder=5)
        ax.add_patch(
            Polygon(
                np.array([E, tangents[0], tangents[1]]),
                closed=True,
                facecolor=BEAM,
                alpha=0.14,
                edgecolor="none",
                zorder=1.5,
            )
        )
        axis_end = E + (d - R - 0.08) * np.array([np.cos(toward), np.sin(toward)])
        ax.plot([E[0], axis_end[0]], [E[1], axis_end[1]], color=BEAM, lw=0.6, ls=":", zorder=5)
        ax.plot(*E, "o", ms=5, mfc="white", mec=BEAM, mew=1.2, zorder=7)  # the firing tile

    # ---- labels: names only, placed in the empty corners outside the orbit
    E0 = R_swarm * unit(emitters[0])
    labels = [
        ax.annotate(
            "collector swarm",
            xy=R_swarm * unit(60) + 0.12 * unit(60),
            xytext=(0.15, 3.05),
            fontsize=8,
            color=SWARM,
            ha="left",
            va="center",
            arrowprops={**LEADER, "color": SWARM},
        ),
        ax.annotate(
            "He II\ndriving zone",
            xy=(r_zone + 0.05) * unit(135),
            xytext=(-4.05, 2.25),
            fontsize=8,
            color=ZONE,
            ha="left",
            va="center",
            arrowprops={**LEADER, "color": ZONE},
        ),
        ax.annotate(
            "hollow-cone\n$\\nu$ beam",
            xy=E0 + 0.7 * unit(emitters[0] + 180),
            xytext=(-4.05, -2.55),
            fontsize=8,
            color=BEAM,
            ha="left",
            va="center",
            arrowprops={**LEADER, "color": BEAM},
        ),
        ax.annotate(
            "emitting\ntile",
            xy=E0 + 0.1 * unit(emitters[0] + 250),
            xytext=(-4.05, -1.35),
            fontsize=8,
            color=BEAM,
            ha="left",
            va="center",
            arrowprops={**LEADER, "color": BEAM},
        ),
    ]
    tag = ax.text(-4.05, 3.25, "(a)", fontsize=10, fontweight="bold", va="top")
    return [*labels, tag], on_colour


def panel_b(ax) -> tuple[list, list]:
    """A shell of thickness H crossed radially and along the grazing chord."""
    ax.set_aspect("equal")
    ax.axis("off")
    r0, H = 1.0, 0.11  # r0 is the outer edge of the driving zone, H its thickness
    r_ph = r0 + 0.22  # the photosphere; the layer above the zone is widened to hold its label
    b = r0 - H  # the grazing ray is tangent to the inner edge of the zone
    x_half = np.sqrt(r0**2 - b**2)
    lo, hi = 38.0, 142.0  # arc window, symmetric about the tangent point at 90 deg
    tt = np.radians(np.linspace(lo, hi, 400))
    r_in = 0.55  # where the drawing stops inward
    # the star below its photosphere, in the colour of panel (a)
    ring = (
        np.r_[r_ph * np.cos(tt), r_in * np.cos(tt[::-1])],
        np.r_[r_ph * np.sin(tt), r_in * np.sin(tt[::-1])],
    )
    ax.add_patch(
        Polygon(np.column_stack(ring), closed=True, facecolor=STAR_FILL, edgecolor="none", zorder=0)
    )
    zone = (
        np.r_[r0 * np.cos(tt), (r0 - H) * np.cos(tt[::-1])],
        np.r_[r0 * np.sin(tt), (r0 - H) * np.sin(tt[::-1])],
    )
    ax.add_patch(
        Polygon(
            np.column_stack(zone),
            closed=True,
            facecolor=ZONE,
            alpha=0.22,
            edgecolor="none",
            zorder=1,
        )
    )
    ax.plot(r_ph * np.cos(tt), r_ph * np.sin(tt), color=STAR_EDGE, lw=1.0, zorder=2)
    for rr in (r0, r0 - H):
        ax.plot(rr * np.cos(tt), rr * np.sin(tt), color=ZONE, lw=1.2, zorder=2)

    # radial ray at 60 degrees: crosses H of shell once
    p_in, p_out = (r0 - H) * unit(60), r0 * unit(60)
    ext = (r_ph + 0.14) * unit(60)
    deep = (r_in + 0.02) * unit(60)
    ax.plot([ext[0], deep[0]], [ext[1], deep[1]], color="0.35", lw=0.8, ls="--", zorder=3)
    ax.plot(
        [p_out[0], p_in[0]], [p_out[1], p_in[1]], color=INK, lw=2.6, solid_capstyle="butt", zorder=4
    )

    # grazing ray tangent to the inner edge: the whole chord through the zone
    ax.plot([-x_half - 0.3, x_half + 0.3], [b, b], color="0.35", lw=0.8, ls="--", zorder=3)
    ax.plot([-x_half, x_half], [b, b], color=BEAM, lw=2.6, solid_capstyle="butt", zorder=4)

    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(0.17, 1.52)
    labels = [
        ax.annotate(
            "radial ray",
            xy=ext,
            xytext=(0.5, 1.42),
            fontsize=8,
            color=INK,
            ha="left",
            va="center",
            arrowprops={**LEADER, "color": INK},
        ),
        ax.annotate(
            "grazing ray",
            xy=(-x_half - 0.25, b),
            xytext=(-1.22, 1.25),
            fontsize=8,
            color=BEAM,
            ha="left",
            va="center",
            arrowprops={**LEADER, "color": BEAM},
        ),
        ax.annotate(
            "driving zone",
            xy=(r0 - 0.5 * H) * unit(136),
            xytext=(-1.22, 0.25),
            fontsize=8,
            color=ZONE,
            ha="left",
            va="center",
            arrowprops={**LEADER, "color": ZONE},
        ),
    ]
    # the two lengths, set right against their segments: H below and right of the
    # radial crossing, the grazing chord's length just below the chord
    on_colour = [
        ax.text(
            *(p_in + np.array([0.02, -0.125])),
            r"$\mathbfit{H}$",
            fontsize=9,
            color=INK,
            ha="left",
            va="top",
            zorder=5,
        ),
        ax.text(
            0.0,
            b - 0.07,
            r"$\mathbfit{2\sqrt{2\,r H}}$",
            fontsize=9,
            color=BEAM,
            ha="center",
            va="top",
            zorder=5,
        ),
        ax.text(
            0.0,
            0.5 * (r0 + r_ph) + 0.005,
            "photosphere",
            fontsize=7.5,
            color=STAR_INK,
            ha="center",
            va="center",
            zorder=5,
        ),
    ]
    tag = ax.text(-1.22, 1.52, "(b)", fontsize=10, fontweight="bold", va="top")
    return [*labels, tag], on_colour


def drawn_points(ax, exclude: set, fills: bool = True) -> np.ndarray:
    """Display-space points densely sampled along every line and patch edge,
    and, with `fills`, over the interior of every filled shape."""
    pts = []
    for art in list(ax.lines) + list(ax.patches):
        if art in exclude:
            continue
        if isinstance(art, FancyArrowPatch):
            # an arrow's path is only its rendered outline; its shaft is posA -> posB
            verts = ax.transData.transform(np.asarray(art._posA_posB, float))
        elif hasattr(art, "get_xydata"):
            verts = ax.transData.transform(np.asarray(art.get_xydata(), float))
        elif art.get_linewidth() == 0 or art.get_edgecolor()[3] == 0:
            verts = np.empty((0, 2))  # an unstroked fill has no outline to hit
        else:
            verts = None
        # a patch outline is sampled polygon by polygon, so a CLOSEPOLY or a
        # jump between sub-paths is never read as a stroke
        polys = (
            [verts]
            if verts is not None
            else art.get_transform().transform_path(art.get_path()).to_polygons(closed_only=False)
        )
        for poly in polys:
            for a, b in pairwise(poly):
                n = max(2, int(np.hypot(*(b - a)) / 1.0))
                pts.append(np.linspace(a, b, n))
        if fills and isinstance(art, Polygon | Circle | Wedge):
            # filled shapes: also test the interior with a coarse grid
            ext = art.get_window_extent()
            gx, gy = np.meshgrid(np.linspace(ext.x0, ext.x1, 60), np.linspace(ext.y0, ext.y1, 60))
            grid = np.column_stack([gx.ravel(), gy.ravel()])
            inside = art.get_path().contains_points(grid, transform=art.get_transform())
            pts.append(grid[inside])
    return np.vstack(pts) if pts else np.empty((0, 2))


def check_labels(fig, axes_labels: list[tuple]) -> None:
    """Refuse to write the figure if a label's text box touches the drawing.

    `axes_labels` holds (axes, labels, on_colour) triples. A label must clear
    every line, edge and filled shape; an on-colour label sits on the fill it
    names, so it only has to clear the lines and edges.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bad = []
    for ax, labels, on_colour in axes_labels:
        leaders = {lab.arrow_patch for lab in labels if getattr(lab, "arrow_patch", None)}
        everything = drawn_points(ax, leaders)
        lines_only = drawn_points(ax, leaders, fills=False)
        for lab in [*labels, *on_colour]:
            pts = lines_only if lab in on_colour else everything
            # Text's own extent: an Annotation's would include its leader
            box = Text.get_window_extent(lab, renderer)
            bb = Bbox.from_extents(box.x0 - 1, box.y0 - 1, box.x1 + 1, box.y1 + 1)
            hit = pts[
                (pts[:, 0] >= bb.x0)
                & (pts[:, 0] <= bb.x1)
                & (pts[:, 1] >= bb.y0)
                & (pts[:, 1] <= bb.y1)
            ]
            if len(hit):
                bad.append(lab.get_text().replace("\n", " "))
            if not ax.bbox.contains(box.x0, box.y0) or not ax.bbox.contains(box.x1, box.y1):
                bad.append(lab.get_text().replace("\n", " ") + " (outside its panel)")
    if bad:
        raise SystemExit("labels overlap the drawing: " + "; ".join(bad))


def main() -> None:
    plt.rcParams.update({"font.size": 8.5, "font.family": "serif"})
    fig, (ax, inset) = plt.subplots(
        1, 2, figsize=(7.0, 3.3), gridspec_kw={"width_ratios": [1.05, 1.0], "wspace": 0.06}
    )
    labels_a, colour_a = panel_a(ax)
    labels_b, colour_b = panel_b(inset)
    check_labels(fig, [(ax, labels_a, colour_a), (inset, labels_b, colour_b)])
    fig.savefig(OUT, bbox_inches="tight")
    fig.savefig(OUT.with_suffix(".png"), bbox_inches="tight", dpi=170)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
