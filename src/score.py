"""SETI modulation score. Matches candidates.html scoreStar."""

from __future__ import annotations

import json
import math
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

CLASS_BASE = {
    "SURVIVOR": 0.70,
    "WEAK": 0.45,
    "TRIAGE": 0.32,
    "QUIET": 0.06,
    "CATALOG": 0,
    "ASTRO": 0.12,
    "RED NOISE": 0.10,
    "P-MULTIPLE": 0.08,
    "SYSTEMATICS": 0.07,
    "SHAPE": 0.06,
    "ARTIFACT": 0.04,
    "ERROR": 0.0,
}
WHY = {
    "SURVIVOR": "survived every stage-2 veto",
    "WEAK": "non-white O–C, not independently confirmed",
    "TRIAGE": "stage-1 flag; not yet vetted",
    "QUIET": "no stage-1 follow-up flag",
    "CATALOG": "in the stage-0 catalog; not yet timed",
    "ASTRO": "explained by known astrophysics",
    "RED NOISE": "short-range correlated noise",
    "P-MULTIPLE": "2P/4P sampling or period doubling",
    "SYSTEMATICS": "TESS orbit/sector period",
    "SHAPE": "light-curve shape change, not timing",
    "ARTIFACT": "structure is a pipeline artifact",
    "ERROR": "download or analysis failed",
}
_VETO = {"ASTRO", "RED NOISE", "P-MULTIPLE", "SYSTEMATICS", "SHAPE", "ARTIFACT", "ERROR"}
# Sparse timing is intrinsically weaker evidence than per-cycle timing: the
# O-C is block-averaged, there are far fewer points per cycle, and ground-based
# sampling aliases. This penalty is what is left of that after stage 2 grew a
# real sparse battery (re-blocking, richer template, block bootstrap, season
# consistency, LTTE, labelled secondary periods) -- it used to be 0.40/0.10/0.08
# and was standing in for tests that were never run at all.
# NOTE: these want re-calibrating against the sparse/dense SURVIVOR rates from a
# full re-vet; they are a defensible starting point, not a measured one.
_SPARSE_PENALTY = {"SURVIVOR": 0.15, "WEAK": 0.05, "TRIAGE": 0.04}
SEQ_P_GATE = 0.01
SEQ_BONUS_SCALE = 0.20


def _num(v: object) -> float | None:
    if v is None or v == "":
        return None
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def _pval(v: object) -> float:
    x = _num(v)
    return 1.0 if x is None else x


def _truthy(v: object) -> bool:
    return v is True or v == 1 or v == "1" or v == "true"


def _js_fixed(x: float, n: int = 4) -> float:
    return float(f"{x:.{n}f}")


def _json_obj(raw: str | bytes | dict | None) -> dict:
    if not raw:
        return {}
    if isinstance(raw, dict):
        return raw
    try:
        return json.loads(raw)
    except (TypeError, json.JSONDecodeError, ValueError):
        return {}


EVIDENCE_GATE = 0.05  # above this, a p-value is no evidence at all
EVIDENCE_FULL_LOG10 = -6.0  # log10 p at which the curve saturates
_LOG10_GATE = math.log10(EVIDENCE_GATE)
# Only reached when a caller passes a p of exactly 0 with no log companion;
# well past saturation, so the precise value never matters.
_LOG10_MIN = -300.0


def _evidence_log10(log10_p: object) -> float:
    """Log-evidence curve, taking log10(p) directly.

    Working in log space is what lets an underflowed analytic tail (scipy's
    .sf() returns 0.0 for an extreme statistic) stay distinguishable from a
    p-value that really sat at 1e-12.
    """
    if log10_p == -math.inf:
        return 1.0
    x = _num(log10_p)
    if x is None or x >= _LOG10_GATE:
        return 0.0
    return min(1.0, max(0.0, x / EVIDENCE_FULL_LOG10))


def _evidence(p: object, log10_p: object = None) -> float:
    if log10_p is not None and _num(log10_p) is not None:
        return _evidence_log10(log10_p)
    x = _pval(p)
    if x >= EVIDENCE_GATE:
        return 0.0
    return _evidence_log10(math.log10(x) if x > 0 else _LOG10_MIN)


def _mean_evidence(terms: Iterable[object]) -> float:
    """Mean log-evidence over (p, log10_p) terms; a bare p is also accepted."""
    xs = []
    for term in terms:
        p, log10_p = term if isinstance(term, tuple) else (term, None)
        if _num(log10_p) is not None:
            xs.append(_evidence_log10(log10_p))
            continue
        if p is None or p == "":
            continue
        if _num(p) is None:
            continue
        xs.append(_evidence(p))
    if not xs:
        return 0.0
    return sum(xs) / len(xs)


def _snr_factor(rms: object, err: object) -> float:
    r, e = _num(rms), _num(err)
    if r is None or e is None or e <= 0:
        return 1.0
    return min(1.0, max(0.0, r / (3.0 * e)))


# Secondary-period kinds that must not cost the score a second time. Each
# either drove a hard veto (so the _VETO branch has already capped the score),
# says the period is not modulation at all, is a harmonic whose fundamental is
# charged instead, or -- for a code candidate -- is the signal rather than a
# contaminant.
_SEC_ACCOUNTED_KINDS = frozenset(
    {"systematics", "p-multiple", "beat", "astro-range", "secular", "harmonic", "code-candidate"}
)
# Fallback for payloads written before `kind` existed.
_SEC_ACCOUNTED_LABELS = ("systematics", "P-multiple", "(ASTRO)", "Blazhko", "longer than baseline")


def _real_secondary_periods(secs: Sequence[dict] | None) -> list[dict]:
    """Secondary periods that should still cost the score.

    alias_suspect is only ever set for ground-based data -- is_alias() returns
    False outright when the survey is not ground-based -- so for dense TESS
    candidates it is always False and discriminates nothing. Stage 2 attaches a
    `kind` instead; only the kinds that are genuinely unexplained coherent
    modulation are charged. Stage-1-only rows carry neither and keep the
    original behaviour.
    """
    out = []
    for s in secs or []:
        if not s or s.get("alias_suspect"):
            continue
        kind = s.get("kind")
        if kind is not None:
            if kind not in _SEC_ACCOUNTED_KINDS:
                out.append(s)
            continue
        label = str(s.get("label") or "")
        if any(key in label for key in _SEC_ACCOUNTED_LABELS):
            continue
        out.append(s)
    return out


def _fmt_p(p: float) -> str:
    if p == 0:
        return "0"
    return f"{p:.2g}"


def sequence_bonus(p_seq: object) -> float:
    """Score add-on for a catalog bit-string match with p_seq < 0.01.

    Same log-evidence curve as Ljung-Box/ACF (full scale at 1e-6), times 0.20.
    """
    x = _pval(p_seq)
    if x >= SEQ_P_GATE:
        return 0.0
    e = min(1.0, max(0.0, math.log10(max(x, 1e-12)) / -6.0))
    return SEQ_BONUS_SCALE * e


def classify(
    s1: Mapping[str, Any], s2: Mapping[str, Any] | None, p1: Mapping[str, Any] | None
) -> str:
    if s1.get("error") or (p1 or {}).get("error"):
        return "ERROR"
    if s2 and s2.get("classification"):
        return str(s2["classification"])
    if _truthy(s1.get("in_triage")):
        return "TRIAGE"
    return "QUIET"


def score_star(
    s1: dict, s2: dict | None = None, p1: dict | None = None, p2: dict | None = None
) -> dict:
    s1 = s1 or {}
    p1 = p1 or {}
    p2 = p2 or {}
    s2 = s2 or None
    label = classify(s1, s2, p1)
    base = CLASS_BASE.get(label, 0.05)
    p_lb = p1.get("p_ljung_box", s1.get("p_ljung_box"))
    p_acf = p1.get("p_max_acf_shuffle")
    p_jit = p1.get("p_jitter_bootstrap")
    p_boot = p2.get("block_bootstrap_p")
    # The O-C read as a clock-locked bit word, calibrated against the same
    # shuffle null as the ACF. Without it a star promoted by the decode route
    # alone would carry no evidence and rank below the noise it was picked out
    # of. Absent terms are skipped, so a row written before the statistic
    # existed scores exactly as it did.
    p_dec = p1.get("p_decode_shuffle")
    # The two routes built for a designed beacon: quantised steps (shuffle
    # null) and the folded steps matched against the catalogs (analytic,
    # Bonferroni over every trial). Absent on rows written before they existed.
    p_quant = p1.get("p_quantised_shuffle")
    p_cat = p1.get("p_catalog_matched")
    # Monte-Carlo p-values are floored at 1/(draws+1) so they have no log
    # companion; an analytic tail does, and passing it keeps an underflow from
    # collapsing to p = 0.
    ev = _mean_evidence(
        [
            (p_lb, p1.get("log10_p_ljung_box")),
            (p_acf, None),
            (p_boot, None),
            (p_dec, None),
            (p_quant, None),
            (p_cat, None if p_cat is None or p_cat <= 0 else math.log10(p_cat)),
        ]
    )
    snr = _snr_factor(
        p1.get("oc_rms_seconds", s1.get("oc_rms_seconds")), p1.get("median_timing_err_seconds")
    )
    ev_snr = ev * snr
    robust = 0.0
    notes = []
    bshift = p2.get("boundary_shift") or {}
    rich = p2.get("rich_template") or {}
    indep = p2.get("independent_timing") or {}
    segs = p2.get("segments") or {}
    if bshift.get("structure_persists"):
        robust += 0.06
        notes.append("O–C survives cycle-boundary shift")
    if rich.get("structure_persists"):
        robust += 0.06
        notes.append("O–C survives richer template")
    if indep.get("agrees") is True:
        robust += 0.08
        notes.append("independent timing agrees")
    elif indep.get("agrees") is False:
        robust -= 0.08
        notes.append("independent timing disagrees")
    if segs.get("independent_confirmation"):
        robust += 0.06
        notes.append("structure in ≥2 sectors")
    secs = p2.get("secondary_periods") or p1.get("secondary_periods") or []
    if _real_secondary_periods(secs):
        robust -= 0.12
        notes.append("coherent secondary period (likely multi-mode / Blazhko)")
    if snr < 1:
        notes.append("O–C rms is " + f"{snr * 3:.1f}" + "× the timing error")
    mode = str((p1 or {}).get("mode") or s1.get("mode") or "")
    sparse_pen = 0.0
    if mode == "sparse":
        sparse_pen = _SPARSE_PENALTY.get(label, 0.0)
        if sparse_pen:
            notes.append("sparse timing (weaker per-cycle / vet checks)")
    p_seq = s1.get("seq_p")
    if p_seq is None:
        p_seq = p1.get("seq_p")
    veto = label in _VETO
    seq_bonus = 0.0 if (veto or label == "ERROR") else sequence_bonus(p_seq)
    if seq_bonus:
        notes.append("sequence p=" + _fmt_p(_pval(p_seq)))
    if label == "ERROR":
        score = 0.0
    elif veto:
        score = min(1.0, max(0.0, base + 0.10 * ev_snr + min(robust, 0.04)))
    else:
        score = min(1.0, max(0.0, base + 0.30 * ev_snr + robust + seq_bonus))
    if sparse_pen:
        score = max(0.0, score - sparse_pen)
    parts = {
        "class_base": _js_fixed(base),
        "evidence": _js_fixed(ev),
        "snr": _js_fixed(snr),
        "robust": _js_fixed(robust),
        "sparse_pen": _js_fixed(sparse_pen),
        "seq_bonus": _js_fixed(seq_bonus),
        "seq_p": _pval(p_seq),
        "why": WHY.get(label, ""),
        "notes": notes,
        "p_min": min(_pval(p_lb), _pval(p_acf), _pval(p_jit), _pval(p_boot)),
    }
    return {
        "score": score,
        "score_pct": round(score * 1000) / 10,
        "classification": label,
        "p_min": parts["p_min"],
        "parts": parts,
        "stage2": bool(s2 and s2.get("classification")),
    }


def score_from_payloads(
    s1: dict, s1_payload: str | dict | None, s2: dict | None, s2_payload: str | dict | None
) -> dict:
    p1 = _json_obj(s1_payload)
    p2 = _json_obj(s2_payload)
    s2_row = None
    if s2 or p2:
        s2_row = dict(s2 or {})
        if p2.get("classification") and not s2_row.get("classification"):
            s2_row["classification"] = p2["classification"]
        if not s2_row.get("classification") and not p2:
            s2_row = None
    return score_star(s1, s2_row, p1, p2)


# ---- presentation of the stage-2 test battery -----------------------------
# The pass/fail thresholds mirror stage2_vet's own verdict logic; they live
# here, next to the score, rather than being restated in the viewer.
VET_BOOTSTRAP_MAX_P = 0.01
VET_JITTER_MAX_P = 1e-3
VET_JITTER_MIN_SNR = 4.0
VET_LTTE_MAX_P = 1e-4
VET_LTTE_MIN_VAR = 0.5


def catalog_parts() -> dict:
    """`parts` block for a stage-0 row that has not been timed yet."""
    return {
        "class_base": 0.0,
        "evidence": 0.0,
        "snr": 1.0,
        "robust": 0.0,
        "sparse_pen": 0.0,
        "seq_bonus": 0.0,
        "seq_p": 1.0,
        "why": WHY["CATALOG"],
        "notes": [],
        "p_min": 1.0,
    }


def _fmt(v: object, nd: int | None = None) -> str:
    x = _num(v)
    if x is None:
        return "\u2014"
    return f"{x:.{nd}f}" if nd is not None else str(x)


def _fmt_pval(v: object) -> str:
    x = _num(v)
    if x is None:
        return "\u2014"
    if x == 0:
        return "0"
    if x < 1e-3:
        return f"{x:.2e}"
    return f"{x:.3g}"


def vetting_test_rows(p2: dict | None) -> list[dict]:
    """One row per stage-2 test: name, value, pass flag, detail."""
    p2 = p2 or {}
    if not p2:
        return []
    rows = []
    b = p2.get("boundary_shift")
    if b:
        rows.append(
            {
                "name": "Cycle boundary +P/2",
                "value": "persists" if b.get("structure_persists") else "vanishes",
                "ok": bool(b.get("structure_persists")),
                "detail": "Ljung-Box p=" + _fmt_pval(b.get("p_ljung_box")),
            }
        )
    r = p2.get("rich_template")
    if r:
        rows.append(
            {
                "name": f"Richer template ({r.get('nharm', '?')} harm.)",
                "value": "persists" if r.get("structure_persists") else "vanishes",
                "ok": bool(r.get("structure_persists")),
                "detail": (
                    "rms="
                    + _fmt(r.get("oc_rms_seconds"), 1)
                    + " s  p="
                    + _fmt_pval(r.get("p_ljung_box"))
                ),
            }
        )
    ind = p2.get("independent_timing")
    if ind:
        agrees = ind.get("agrees")
        rows.append(
            {
                "name": "Independent steepest-branch timing",
                "value": (
                    "agrees"
                    if agrees is True
                    else "disagrees"
                    if agrees is False
                    else f"n={ind.get('n', 0)}"
                ),
                "ok": agrees is True,
                "detail": (
                    "\u03c1="
                    + _fmt(ind.get("spearman_vs_template_oc"), 2)
                    + "  p="
                    + _fmt_pval(ind.get("p_spearman"))
                ),
            }
        )
    base = p2.get("baseline") or {}
    if base.get("p_quantised_shuffle") is not None:
        pq = float(base["p_quantised_shuffle"])
        rows.append(
            {
                "name": "Quantised steps",
                "value": "two step sizes" if pq < 1e-3 else "continuous",
                "ok": pq < 1e-3,
                "detail": "Otsu share="
                + _fmt(base.get("quantised_stat"), 2)
                + "  p="
                + _fmt_pval(pq),
            }
        )
    if base.get("p_catalog_matched") is not None:
        pc = float(base["p_catalog_matched"])
        rows.append(
            {
                "name": "Folded steps vs catalog",
                "value": str(base.get("catalog_matched_label") or "no match")
                if pc < 1e-3
                else "no match",
                "ok": pc < 1e-3,
                "detail": "p=" + _fmt_pval(pc),
            }
        )
    fs = p2.get("fold_stability")
    if fs and fs.get("tested"):
        rows.append(
            {
                "name": "Word stable across halves",
                "value": "same word" if fs.get("stable") else "differs",
                "ok": bool(fs.get("stable")),
                "detail": "agreement="
                + _fmt(fs.get("agreement"), 2)
                + (
                    f"  ({fs['first']['n']} vs {fs['second']['n']} bits)"
                    if fs.get("first") and fs.get("second")
                    else ""
                ),
            }
        )
    sh = p2.get("shape_change")
    if sh:
        rows.append(
            {
                "name": "Shape vs timing",
                "value": "shape-driven" if sh.get("shape_driven") else "timing-like",
                "ok": not sh.get("shape_driven"),
                "detail": (
                    "|O\u2013C|\u2013\u03c7\u00b2 \u03c1="
                    + _fmt(sh.get("spearman_absOC_vs_chi2"), 2)
                    + "  O\u2013C\u2013amp \u03c1="
                    + _fmt(sh.get("spearman_OC_vs_amplitude"), 2)
                ),
            }
        )
    if p2.get("block_bootstrap_p") is not None:
        rows.append(
            {
                "name": "Block-bootstrap ACF",
                "value": _fmt_pval(p2.get("block_bootstrap_p")),
                "ok": _pval(p2.get("block_bootstrap_p")) <= VET_BOOTSTRAP_MAX_P,
                "detail": "long-range structure beyond short-range red noise",
            }
        )
    segs = p2.get("segments")
    if segs:
        rows.append(
            {
                "name": "Per-sector confirmation",
                "value": (
                    f"{segs.get('n_significant', 0)}/{segs.get('n_segments_tested', 0)} sectors"
                ),
                "ok": bool(segs.get("independent_confirmation")),
                "detail": "structure in at least two independent segments",
            }
        )
    ltte = p2.get("ltte")
    if ltte:
        astro = (
            _pval(ltte.get("p_ftest")) < VET_LTTE_MAX_P
            and (_num(ltte.get("variance_explained")) or 0) > VET_LTTE_MIN_VAR
        )
        rows.append(
            {
                "name": "LTTE / long sinusoid",
                "value": (
                    "P="
                    + _fmt(ltte.get("period_days"), 2)
                    + " d  A="
                    + _fmt(ltte.get("semi_amplitude_seconds"), 1)
                    + " s"
                ),
                "ok": not astro,
                "detail": (
                    "F-test p="
                    + _fmt_pval(ltte.get("p_ftest"))
                    + "  var="
                    + _fmt(ltte.get("variance_explained"), 2)
                ),
            }
        )
    jit = p2.get("baseline_jitter")
    if jit:
        rows.append(
            {
                "name": "Sparse-mode jitter",
                "value": _fmt(jit.get("jitter_sigma_t_seconds"), 1) + " s",
                "ok": (
                    _pval(jit.get("p_jitter_bootstrap")) < VET_JITTER_MAX_P
                    and (_num(jit.get("jitter_var_t_over_err")) or 0) > VET_JITTER_MIN_SNR
                ),
                "detail": (
                    "p="
                    + _fmt_pval(jit.get("p_jitter_bootstrap"))
                    + "  var/err="
                    + _fmt(jit.get("jitter_var_t_over_err"), 2)
                ),
            }
        )
    return rows
