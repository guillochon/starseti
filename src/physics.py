"""What driving a pulsator costs, and what a transmitter can afford.

The numbers behind the viewer's energy-per-bit box and the injection report's
power section, derived in research/modulation-physics.md and corrected for the
neutrino cross-section in research/response1.md (M1). Everything here is
stdlib-only so `serve_candidates.py` can import it next to `sequence.py`.

The model, in one line: a set bit is a triggered early pulse, and triggering it
costs the beam

    E_trigger = L * dt / (G * f_abs * f_zone)

where L*dt is the stellar energy the phase advance is worth, G is the gain
Learned et al. invoke and never quantify (taken as 1 here), f_abs is the
fraction of the beam the star absorbs at all and f_zone the fraction of *that*
which lands in the He II driving zone. The product 1/(G f_abs f_zone) is the
coupling constant the viewer used to carry as an unexplained 0.1; it is a
*divisor* and is >= 1 unless the star amplifies.

Three beam geometries, because f_zone and f_abs are where the physics lives:

* **central** -- a beam through the disk. Neutrinos deposit in proportion to
  mass, so nearly everything lands in the core and f_zone is the mass fraction
  above the driving zone, ~6e-4 for a Cepheid. The coupling is in the
  thousands and a transmitter at unit gain has to beam more than the star's
  own luminosity, which a transmitter powered by that star cannot collect.
* **annular** -- a hollow cone at the limb, impact parameter at the *base*
  of the driving zone, the layer whose overlying thermal time equals one
  period (the Cox transition-region criterion). A ray at impact parameter b
  samples only radii >= b, so the beam crosses the zone and the layers above
  it and nothing else. The column there is set by the star's envelope
  structure: a MESA horizontal-branch model (0.63 Msun, 46 Lsun, 6725 K;
  research/response2.md) puts the base of an RRab's He II zone at ~55,000 K
  with 750 g/cm^2 above it, a density of 6e-8 g/cm^3 and H/R = 0.03, i.e. a
  mass fraction f_above ~ 8e-7 rather than the 1e-3 the hydrostatic guess
  before 2026-09-26 gave; a grey radiative envelope reproduces the model to a
  factor of two and supplies the Cepheid (1e-4) and beta Cep (1.4e-6) values.
  Every grazing limb is therefore optically thin to any neutrino: the tunable
  beam at 1 EeV absorbs 1e-4 to 2e-3 of its energy, and the design that
  applies to every class is the resonant one below. Because tau << 1 the
  deposition follows the local density along the chord, 84 per cent of it
  within one scale height above the tangent point, i.e. inside the zone whose
  base the ray grazes: f_zone = F_ZONE_ANNULAR = 1. (Review 2's f_zone = 0.5
  was for a tau ~ 1 path, which no longer exists.)
* **glashow** -- the annular geometry with a monochromatic anti-electron-
  neutrino beam at the Glashow resonance (6.3 PeV, sigma ~ 5e-31 cm^2 per
  electron; IceCube saw one such event, Nature 591, 220), the largest
  cross-section a neutrino has. On the grazing columns at the base of the
  zone it gives tau = f_abs of 0.023 (Cepheid), 3e-3 (RRab) and 0.064 (beta
  Cep), thirty times the 1 EeV beam, and it is the default geometry.

The transmitter is in the system -- at P_beam/L ~ 1e-3 the only plausible
power source is the star -- so it is a collector swarm intercepting a fraction
f_coll of the star's sky and converting it at efficiency eta:

    P_beam / L = duty * (dt/P) / (G f_abs f_zone) = duty * (dt/P) * coupling
    f_coll     = P_beam / (L * eta)

That puts a ceiling on the phase code itself: a self-powered transmitter cannot
drive dt/P above eta * f_coll / (duty * coupling), and with f_coll <= 1 (a
complete Dyson sphere) that is the largest amplitude the hypothesis allows.
The injection search measures the amplitude the *pipeline* needs; this module
says whether a transmitter could have supplied it.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence

LSUN_W = 3.826e26  # W
DAY_S = 86400.0
YEAR_S = 3.156e7

# The gain the star provides beyond delivering the energy. Learned et al. call
# the star "a huge signal amplifier" and decline to quantify it; the whole of
# this module is at G = 1, and every figure scales as 1/G.
GAIN = 1.0
# Fraction of the energy an annular (or Glashow) beam deposits that lands in
# the driving zone. The limb is optically thin, so the deposition follows the
# density along the chord and sits within a scale height of the tangent point
# at the base of the zone (module docstring): effectively all of it.
F_ZONE_ANNULAR = 1.0

# The fraction of a beam absorbed at the energy where the path is one
# interaction length deep -- 1 - 1/e, rounded the way the note rounds it. On a
# path that never reaches tau = 1 the absorbed fraction is tau itself.
F_ABS_TAU1 = 0.6

# Fraction of the star's time a message keeps the bit set. A random word at
# fair-coin density; the injection grids use the same value.
DUTY = 0.5

# Efficiencies of turning collected starlight into a beamed neutrino flux at
# which the collector fraction is quoted. Nobody knows this number; unity is
# the physical limit and the paper allows it, so that an RRab transmitter is
# judged against the whole of its star's output rather than a tenth of it.
ETA_GRID = (1.0, 0.1, 0.01)

GEOMETRIES = ("annular", "central", "glashow")

# What a transmitter with no energy store can put into one kick: one cycle's
# intake, eta x f_coll x L x P. The paper's ceiling is the physical maximum, a
# complete sphere at unit efficiency, so that dt/P <= G f_abs f_zone whatever
# the duty cycle; it is the number a measured threshold is compared with.
STORAGE_ETA = 1.0
STORAGE_F_COLL = 1.0
DEFAULT_GEOMETRY = "glashow"

# --- the neutrino-nucleon cross-section ------------------------------------
#
# Linear in energy at and below a few TeV, sigma ~ 7e-36 cm^2 (E / TeV); above
# ~10 PeV the W propagator turns the deep-inelastic cross-section over to
# sigma_CC+NC ~ 7.8e-36 cm^2 (E / GeV)^0.363 (Gandhi, Quigg, Reno & Sarcevic
# 1996, Astropart. Phys. 5, 81; 1998, PRD 58, 093009). The two fits cross near
# 60 TeV; taking the smaller of the two is within the softening between them
# that Bustamante & Connolly (2019) and IceCube (2021) measure.
SIGMA_1TEV_CM2 = 7e-36
GANDHI_NORM_CM2 = 7.8e-36
GANDHI_INDEX = 0.363
M_P_G = 1.6726e-24

# The highest beam energy the annular design is allowed, by convention: 1 EeV.
# The highest *observed* neutrino energy is the KM3NeT event of 2025
# (~220 PeV; Nature 638, 376), which the paper quotes beside it.
BEAM_ENERGY_MAX_GEV = 1e9
BEAM_ENERGY_OBSERVED_GEV = 2.2e8

# The Glashow resonance: anti-nu_e on an atomic electron forming a W-, peaking
# at 6.3 PeV with ~5e-31 cm^2 per electron. An H/He envelope with X = 0.7
# carries (1 + X) / 2 ~ 0.85 electrons per proton mass.
GLASHOW_ENERGY_GEV = 6.3e6
GLASHOW_SIGMA_CM2 = 5e-31
ELECTRONS_PER_M_P = 0.85


def sigma_nu_n(energy_gev: float) -> float:
    """Neutrino-nucleon cross-section in cm^2 at `energy_gev`."""
    e = float(energy_gev)
    return min(SIGMA_1TEV_CM2 * e / 1e3, GANDHI_NORM_CM2 * e**GANDHI_INDEX)


def energy_for_sigma(sigma_cm2: float) -> float:
    """The beam energy in GeV at which the cross-section reaches `sigma_cm2`.

    `sigma_nu_n` is the smaller of two rising functions, so its inverse is the
    larger of the two inverses. Past ~1e12 GeV the answer is a statement that
    no beam reaches it, not an energy anyone would build.
    """
    s = float(sigma_cm2)
    return max(s / SIGMA_1TEV_CM2 * 1e3, (s / GANDHI_NORM_CM2) ** (1.0 / GANDHI_INDEX))


def format_energy(energy_gev: float) -> str:
    """'9 TeV', '100 PeV', '1.0 EeV', '1.6e+06 EeV' -- for tables and tooltips."""
    e = float(energy_gev)
    for unit, scale in (("EeV", 1e9), ("PeV", 1e6), ("TeV", 1e3), ("GeV", 1.0)):
        if e >= scale:
            v = e / scale
            if v >= 1e4:
                return f"{v:.1e} {unit}"
            return f"{v:.0f} {unit}" if v >= 10 else f"{v:.1f} {unit}"
    return f"{e:.2g} GeV"


# Reference star per class (M/Msun, R/Rsun, L/Lsun, P days) and the two
# quantities of the envelope at the base of the driving zone, the layer whose
# overlying thermal time is one period: the mass fraction above it (`f_above`,
# which is f_zone for a central beam) and its scale height in units of the
# radius (`h_over_r`), which sets the chord enhancement sqrt(2R/H) of a
# grazing ray. The RRab values are read off a MESA horizontal-branch model
# (0.63 Msun, 46 Lsun, 6725 K: 750 g/cm^2 above 55,000 K); the other two come
# from a grey radiative envelope that reproduces that model to a factor of two
# (research/modulation-physics.md, 2026-09-26 correction). The tau = 1
# energies are computed from the columns, which is what made the Table 1
# inconsistency the first referee found impossible to repeat.
REFERENCE = {
    "cepheid": {
        "label": "classical Cepheid",
        "mass_msun": 5.0,
        "radius_rsun": 50.0,
        "lum_lsun": 3000.0,
        "period_days": 10.0,
        "f_above": 1.0e-4,
        "h_over_r": 0.042,
    },
    "rrlyrae": {
        "label": "RR Lyrae (RRab)",
        "mass_msun": 0.7,
        "radius_rsun": 5.0,
        "lum_lsun": 50.0,
        "period_days": 0.57,
        "f_above": 8e-7,
        "h_over_r": 0.030,
    },
    "bcep": {
        "label": "beta Cephei",
        "mass_msun": 10.0,
        "radius_rsun": 7.0,
        "lum_lsun": 15000.0,
        "period_days": 0.20,
        "f_above": 1.4e-6,
        "h_over_r": 0.011,
    },
}
DEFAULT_CLASS = "cepheid"

_G_CGS = 6.674e-8
_MSUN_G = 1.989e33
_RSUN_CM = 6.957e10
_LSUN_CGS = 3.828e33


def reference(star_class: str | None) -> dict:
    """The reference star for a class, defaulting to the Cepheid."""
    return REFERENCE.get(str(star_class or "").lower(), REFERENCE[DEFAULT_CLASS])


def log_g(star_class: str | None) -> float:
    """Surface gravity of the reference star, cgs. f_abs at fixed beam energy
    scales with it, which is why the instability strip is transparent and the
    compact pulsators are not."""
    r = reference(star_class)
    g = _G_CGS * r["mass_msun"] * _MSUN_G / (r["radius_rsun"] * _RSUN_CM) ** 2
    return math.log10(g)


def kelvin_helmholtz_yr(star_class: str | None) -> float:
    """GM^2 / (R L) for the reference star, in years: ~5e3 for the Cepheid and
    ~6e4 for the RRab. The heating of section 2.3 cancels the code's period
    offset after a third of this."""
    r = reference(star_class)
    m = r["mass_msun"] * _MSUN_G
    return _G_CGS * m * m / (r["radius_rsun"] * _RSUN_CM * r["lum_lsun"] * _LSUN_CGS) / YEAR_S


def disk_column(star_class: str | None) -> float:
    """M / (pi R^2), the disk-averaged column a central ray sees, in g cm^-2."""
    r = reference(star_class)
    return r["mass_msun"] * _MSUN_G / (math.pi * (r["radius_rsun"] * _RSUN_CM) ** 2)


def chord_enhancement(star_class: str | None) -> float:
    """sqrt(2R/H): how much longer a grazing ray's chord through the driving
    zone is than a radial ray's (~8 Cepheid, ~10 RRab, ~30 beta Cep)."""
    return math.sqrt(2.0 / float(reference(star_class)["h_over_r"]))


def grazing_column(star_class: str | None) -> float:
    """The column along a chord grazing the driving zone: the mass above the
    zone, f_above of the disk column, stretched along the chord, in g cm^-2.
    Disk over grazing is 2 / (sqrt(2R/H) f_above): ~400, ~200 and ~60."""
    r = reference(star_class)
    return disk_column(star_class) * float(r["f_above"]) * chord_enhancement(star_class) / 2.0


def _check_geometry(geometry: str) -> None:
    if geometry not in GEOMETRIES:
        raise ValueError(f"unknown beam geometry {geometry!r}; one of {GEOMETRIES}")


def column(star_class: str | None, geometry: str = DEFAULT_GEOMETRY) -> float:
    """The column the beam has to be absorbed in, per geometry."""
    _check_geometry(geometry)
    return disk_column(star_class) if geometry == "central" else grazing_column(star_class)


def tau1_energy_gev(star_class: str | None, geometry: str = DEFAULT_GEOMETRY) -> float:
    """The beam energy at which the path is one interaction length deep."""
    _check_geometry(geometry)
    if geometry == "glashow":
        return GLASHOW_ENERGY_GEV
    return energy_for_sigma(M_P_G / column(star_class, geometry))


def beam_energy_gev(
    star_class: str | None, geometry: str = DEFAULT_GEOMETRY, e_max: float = BEAM_ENERGY_MAX_GEV
) -> float:
    """The energy the design runs at: tau = 1 where that is reachable, else the
    highest energy allowed. The resonance has no energy to choose."""
    _check_geometry(geometry)
    if geometry == "glashow":
        return GLASHOW_ENERGY_GEV
    return min(tau1_energy_gev(star_class, geometry), float(e_max))


def optical_depth(
    star_class: str | None,
    geometry: str = DEFAULT_GEOMETRY,
    energy_gev: float | None = None,
    e_max: float = BEAM_ENERGY_MAX_GEV,
) -> float:
    """tau = Sigma sigma / m_p along the geometry's path at the beam energy."""
    _check_geometry(geometry)
    if geometry == "glashow":
        return column(star_class, geometry) * ELECTRONS_PER_M_P * GLASHOW_SIGMA_CM2 / M_P_G
    e = beam_energy_gev(star_class, geometry, e_max) if energy_gev is None else float(energy_gev)
    return column(star_class, geometry) * sigma_nu_n(e) / M_P_G


def f_abs(
    star_class: str | None, geometry: str = DEFAULT_GEOMETRY, e_max: float = BEAM_ENERGY_MAX_GEV
) -> float:
    """Fraction of the beam the star absorbs at all: F_ABS_TAU1 where the path
    reaches tau = 1, and tau itself on a path that stays thin."""
    tau = optical_depth(star_class, geometry, e_max=e_max)
    return F_ABS_TAU1 if tau >= 1.0 - 1e-9 else float(tau)


def f_zone(star_class: str | None, geometry: str = DEFAULT_GEOMETRY) -> float:
    """Fraction of the absorbed beam that lands in the driving zone."""
    _check_geometry(geometry)
    if geometry == "central":
        return float(reference(star_class)["f_above"])
    return F_ZONE_ANNULAR


def coupling(
    star_class: str | None = None,
    geometry: str = DEFAULT_GEOMETRY,
    gain: float = GAIN,
    e_max: float = BEAM_ENERGY_MAX_GEV,
) -> float:
    """1 / (G f_abs f_zone): what a bit costs in units of L*dt.

    Per class for the resonant beam at the base of the zone (~16 beta Cep,
    ~330 RRab, ~44 Cepheid); thirty times more for the tunable beam at 1 EeV;
    in the thousands to millions for a central beam, where it is the
    deposition penalty in person.
    """
    return 1.0 / (float(gain) * f_abs(star_class, geometry, e_max) * f_zone(star_class, geometry))


def beam_power_ratio(frac: float, c: float, duty: float = DUTY) -> float:
    """P_beam / L for a code of amplitude dt/P = frac at coupling c."""
    return float(duty) * float(frac) * float(c)


def collector_fraction(frac: float, c: float, eta: float, duty: float = DUTY) -> float:
    """Fraction of the star's sky a self-powered transmitter has to intercept.

    Above 1 the transmitter would have to collect more than the star emits,
    which no efficiency fixes: the hypothesis is impossible at that amplitude.
    """
    return beam_power_ratio(frac, c, duty) / float(eta)


def frac_ceiling(c: float, eta: float, f_coll: float = 1.0, duty: float = DUTY) -> float:
    """The largest dt/P a self-powered transmitter can drive.

    A complete Dyson sphere (f_coll = 1) at efficiency eta is the most it can
    have; a code the pipeline needs more amplitude than this to notice is one
    the hypothesis cannot supply.
    """
    return float(eta) * float(f_coll) / (float(duty) * float(c))


def kick_ceiling(c: float, eta: float = STORAGE_ETA, f_coll: float = STORAGE_F_COLL) -> float:
    """The largest dt/P a single kick can carry without energy storage.

    Banking cycles of intake for one large kick needs a reservoir of order
    1e29 J per cycle at an RR Lyrae, so a self-powered transmitter kicks
    every cycle with what it collects in one: dt/P <= eta f_coll / coupling.
    The constraint is per kick, so no duty cycle enters (the duty-1 case of
    `frac_ceiling`). With the defaults, a complete sphere at unit efficiency,
    it is the physical maximum: dt/P <= f_abs f_zone, 3e-3 for an RRab.
    """
    return frac_ceiling(c, eta, f_coll, duty=1.0)


def energy_per_bit_lsun_s(frac: float, c: float, lum_lsun: float, period_days: float) -> float:
    """E_bit in solar-luminosity-seconds: (dt/P) * coupling * L * P."""
    return float(frac) * float(c) * float(lum_lsun) * float(period_days) * DAY_S


def mean_power_w(frac: float, c: float, lum_lsun: float, period_days: float) -> float:
    """Average transmitter power for one bit per cycle, in watts."""
    e = energy_per_bit_lsun_s(frac, c, lum_lsun, period_days)
    return e / (float(period_days) * DAY_S) * LSUN_W


def sagan_k(watts: float) -> float | None:
    """Sagan's continuous Kardashev index, K = (log10 P[W] - 6) / 10."""
    if not (watts > 0):
        return None
    return (math.log10(watts) - 6.0) / 10.0


def _energy_label(star_class: str, geometry: str) -> str:
    """The tau = 1 energy as a string, saying so when it is past the ceiling."""
    e1 = tau1_energy_gev(star_class, geometry)
    if e1 <= BEAM_ENERGY_MAX_GEV:
        return format_energy(e1)
    return (
        f"unreachable (tau = 1 needs {format_energy(e1)}; "
        f"run at {format_energy(BEAM_ENERGY_MAX_GEV)})"
    )


def class_table(gain: float = GAIN) -> list[dict]:
    """One row per class for /api/meta: the couplings, the ceilings, the beam.

    The viewer reads `coupling` (annular, the default) for the energy box and
    the rest for the tooltip, so the assumptions travel with the number.
    """
    out = []
    for cls in REFERENCE:
        ref = reference(cls)
        c_ann = coupling(cls, "annular", gain)
        c_cen = coupling(cls, "central", gain)
        c_gla = coupling(cls, "glashow", gain)
        out.append(
            {
                "star_class": cls,
                "label": ref["label"],
                "log_g": round(log_g(cls), 2),
                "f_above": ref["f_above"],
                "h_over_r": ref["h_over_r"],
                "coupling": c_gla,
                "coupling_annular": c_ann,
                "coupling_central": c_cen,
                "coupling_glashow": c_gla,
                "f_abs_annular": f_abs(cls, "annular"),
                "f_abs_central": f_abs(cls, "central"),
                "f_abs_glashow": f_abs(cls, "glashow"),
                "tau_glashow": optical_depth(cls, "glashow"),
                "beam_energy_annular_gev": beam_energy_gev(cls, "annular"),
                "beam_energy_central_gev": beam_energy_gev(cls, "central"),
                "tau1_annular_gev": tau1_energy_gev(cls, "annular"),
                "tau1_central_gev": tau1_energy_gev(cls, "central"),
                "beam_tau1_annular": _energy_label(cls, "annular"),
                "beam_tau1_central": _energy_label(cls, "central"),
                "column_disk": disk_column(cls),
                "column_grazing": grazing_column(cls),
                "chord_enhancement": round(chord_enhancement(cls), 1),
                "kelvin_helmholtz_yr": kelvin_helmholtz_yr(cls),
                "frac_ceiling": {f"{eta:g}": frac_ceiling(c_gla, eta) for eta in ETA_GRID},
                "kick_ceiling": kick_ceiling(c_gla),
                "kick_ceiling_eta": STORAGE_ETA,
                "kick_ceiling_f_coll": STORAGE_F_COLL,
                "reference": {
                    "lum_lsun": ref["lum_lsun"],
                    "period_days": ref["period_days"],
                    "mass_msun": ref["mass_msun"],
                    "radius_rsun": ref["radius_rsun"],
                },
            }
        )
    return out


def viewer_block(gain: float = GAIN) -> dict:
    """Everything the page needs, with the constants named."""
    return {
        "gain": float(gain),
        "f_abs_tau1": F_ABS_TAU1,
        "duty": DUTY,
        "eta_grid": list(ETA_GRID),
        "geometry": DEFAULT_GEOMETRY,
        "beam_energy_max_gev": BEAM_ENERGY_MAX_GEV,
        "beam_energy_observed_gev": BEAM_ENERGY_OBSERVED_GEV,
        "cross_section": (
            f"sigma = {SIGMA_1TEV_CM2:g} cm^2 (E/TeV) at low energy, "
            f"{GANDHI_NORM_CM2:g} cm^2 (E/GeV)^{GANDHI_INDEX:g} above ~10 PeV "
            "(Gandhi et al. 1996, 1998)"
        ),
        "note": "research/modulation-physics.md; research/response1.md (M1)",
        "classes": class_table(gain),
    }


def _pct(x: float) -> str:
    if x >= 10:
        return f"{x:.0f}x"
    if x >= 1:
        return f"{x:.1f}x"
    return f"{100 * x:.1f}%" if x >= 1e-3 else f"{100 * x:.1e}%"


def _sci(x: float) -> str:
    return f"{x:.1e}"


def _times(x: float) -> str:
    """A ratio to the kick ceiling: '20x', '1.3x', '0.05x'."""
    if x >= 10:
        return f"{x:.0f}x"
    if x >= 1:
        return f"{x:.1f}x"
    return f"{x:.2f}x"


def power_lines(thresholds: Iterable[dict], duty: float = DUTY) -> list[str]:
    """The injection report's power section: what each threshold costs.

    For every measured threshold, the beam power in units of the star's own
    luminosity and the share of a Dyson sphere the transmitter would need at
    each eta, for the annular geometry -- and the central-beam figure beside
    it, which is what the paper's own picture implies. A threshold above the
    self-powered ceiling is one the hypothesis cannot supply, and the line
    says so.
    """
    rows = list(thresholds)
    if not rows:
        return []
    out = [
        "## Power the threshold implies",
        "",
        "| class | gate | dt/P | P_beam/L (annular) | "
        + " | ".join(f"f_coll @ eta={eta:g}" for eta in ETA_GRID)
        + " | P_beam/L (central) | E_bit (L_sun s) | Sagan K | supplied? | vs kick ceiling |",
        "|---|---|---:|---:|" + "---:|" * len(ETA_GRID) + "---:|---:|---:|---|---:|",
    ]
    for t in sorted(rows, key=lambda r: (str(r["star_class"]), str(r["gate"]))):
        cls = str(t["star_class"])
        frac = float(t["frac"])
        bound = bool(t.get("is_bound"))
        ref = reference(cls)
        c_ann = coupling(cls, DEFAULT_GEOMETRY)
        c_cen = coupling(cls, "central")
        ratio = beam_power_ratio(frac, c_ann, duty)
        colls = [collector_fraction(frac, c_ann, eta, duty) for eta in ETA_GRID]
        e_bit = energy_per_bit_lsun_s(frac, c_ann, ref["lum_lsun"], ref["period_days"])
        k = sagan_k(mean_power_w(frac, c_ann, ref["lum_lsun"], ref["period_days"]))
        ge = "&ge; " if bound else ""
        best_eta = max(ETA_GRID)
        if collector_fraction(frac, c_ann, best_eta, duty) > 1.0:
            supplied = f"no: needs > 1 sky even at eta = {best_eta:g}"
        else:
            ok = [eta for eta, f in zip(ETA_GRID, colls) if f <= 1.0]
            supplied = f"yes, down to eta = {min(ok):g}"
        out.append(
            f"| {cls} | {t['gate']} | {ge}{_sci(frac)} | {ge}{_sci(ratio)} | "
            + " | ".join(ge + _pct(f) for f in colls)
            + f" | {ge}{_sci(beam_power_ratio(frac, c_cen, duty))} | {ge}{e_bit:.1e} "
            f"| {ge}{'--' if k is None else f'{k:.1f}'} | {supplied} "
            f"| {ge}{_times(frac / kick_ceiling(c_ann))} |"
        )
    ann = ", ".join(
        f"{coupling(c, DEFAULT_GEOMETRY):.3g} for {c} (tau = f_abs = "
        f"{f_abs(c, DEFAULT_GEOMETRY):.2g})"
        for c in REFERENCE
    )
    kicks = ", ".join(
        f"{kick_ceiling(coupling(c, DEFAULT_GEOMETRY)):.1e} for {c}" for c in REFERENCE
    )
    out += [
        "",
        f"At unit gain and duty {duty:g}, with a resonant anti-electron-neutrino beam at "
        f"{format_energy(GLASHOW_ENERGY_GEV)} grazing the base of the driving zone, whose "
        f"column is thin for every class (coupling 1/(f_abs f_zone) = {ann}) and a central beam "
        "(coupling "
        + ", ".join(f"{coupling(c, 'central'):.0f} for {c}" for c in REFERENCE)
        + "). `P_beam/L` is the beam power in units of the star's own luminosity; "
        "`f_coll` is the fraction of the star's sky a self-powered collector would have "
        "to intercept at each conversion efficiency, so a value above 1 is a transmitter "
        "that has to gather more than the star emits. `E_bit` and the Sagan K are for the "
        "class's reference star (see research/modulation-physics.md); the viewer computes "
        "them per star. A `&ge;` row inherits the bound from the threshold. `vs kick "
        f"ceiling` is the threshold over the largest dt/P a transmitter with no energy "
        f"store can kick in one cycle (eta = {STORAGE_ETA:g}, collecting {STORAGE_F_COLL:.0%} "
        f"of the star's sky: dt/P = {kicks}); "
        "a beacon designed under that constraint sits that many times below where the "
        "search listens.",
        "",
    ]
    return out


def ceiling_lines(star_classes: Sequence[str] = tuple(REFERENCE), duty: float = DUTY) -> list[str]:
    """The amplitude a self-powered transmitter can reach, per class and eta."""
    out = [
        "| geometry | coupling | "
        + " | ".join(f"dt/P ceiling @ eta={eta:g}" for eta in ETA_GRID)
        + " |",
        "|---|---:|" + "---:|" * len(ETA_GRID),
    ]
    for cls in star_classes:
        c_ann = coupling(cls, DEFAULT_GEOMETRY)
        out.append(
            f"| annular, {cls} | {c_ann:.3g} | "
            + " | ".join(_sci(frac_ceiling(c_ann, eta, 1.0, duty)) for eta in ETA_GRID)
            + " |"
        )
    for cls in star_classes:
        c_cen = coupling(cls, "central")
        out.append(
            f"| central, {cls} | {c_cen:.0f} | "
            + " | ".join(_sci(frac_ceiling(c_cen, eta, 1.0, duty)) for eta in ETA_GRID)
            + " |"
        )
    return out
