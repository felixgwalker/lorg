from dataclasses import dataclass
import numpy as np


@dataclass
class ViabilityResult:
    delta_froh: float
    froh_reduction_pct: float
    viability_gain: float   # logistic function output in [0, 1]
    recommended: bool
    risk_category: str      # "low" | "moderate" | "high"
    froh_before: float
    froh_after: float
    ne_before: float
    ne_after: float


def score_rescue(
    froh_before: float,
    froh_after: float,
    ne_before: float,
    ne_after: float,
) -> ViabilityResult:
    """Score the benefit of a genetic rescue intervention.

    Parameters
    ----------
    froh_before : FROH of recipient population before rescue.
    froh_after  : projected FROH after rescue / admixture.
    ne_before   : effective population size before rescue.
    ne_after    : projected effective population size after rescue.

    Returns
    -------
    ViabilityResult dataclass with delta_froh, froh_reduction_pct,
    viability_gain (logistic, range [0,1]), recommendation flag, and
    risk_category for the post-rescue population.
    """
    delta_froh = float(froh_before - froh_after)

    if froh_before > 0:
        froh_reduction_pct = delta_froh / froh_before * 100.0
    else:
        froh_reduction_pct = 0.0

    # Logistic viability gain: 1 / (1 + exp(-10 * (delta_froh - 0.1)))
    # Centred at delta_froh = 0.1 with steep slope
    exponent = -10.0 * (delta_froh - 0.1)
    viability_gain = float(1.0 / (1.0 + np.exp(exponent)))

    recommended = viability_gain > 0.15

    if froh_after < 0.125:
        risk_category = "low"
    elif froh_after < 0.25:
        risk_category = "moderate"
    else:
        risk_category = "high"

    return ViabilityResult(
        delta_froh=round(delta_froh, 6),
        froh_reduction_pct=round(froh_reduction_pct, 4),
        viability_gain=round(viability_gain, 6),
        recommended=recommended,
        risk_category=risk_category,
        froh_before=round(float(froh_before), 6),
        froh_after=round(float(froh_after), 6),
        ne_before=float(ne_before),
        ne_after=float(ne_after),
    )


# ---------------------------------------------------------------------------
# Legacy scoring used by the pipeline
# ---------------------------------------------------------------------------

GRADE_THRESHOLDS = [
    (80, "A", "Excellent viability improvement; rescue highly effective."),
    (60, "B", "Good viability improvement; rescue likely beneficial."),
    (40, "C", "Moderate viability improvement; rescue may help."),
    (20, "D", "Low viability improvement; rescue has limited impact."),
    (0,  "F", "Minimal or no viability improvement detected."),
]


def compute_viability_score(sim_result, initial_F):
    mean_F_final = float(sim_result["mean_F"][-1])
    if initial_F <= 0:
        score = 0.0
    else:
        improvement = (initial_F - mean_F_final) / initial_F
        score = float(np.clip(improvement * 100, 0, 100))

    grade = "F"
    interpretation = GRADE_THRESHOLDS[-1][2]
    for threshold, g, interp in GRADE_THRESHOLDS:
        if score >= threshold:
            grade = g
            interpretation = interp
            break

    mean_H_initial = float(sim_result["mean_H"][0])
    mean_H_final = float(sim_result["mean_H"][-1])

    return {
        "score": round(score, 2),
        "grade": grade,
        "interpretation": interpretation,
        "initial_F": round(float(initial_F), 4),
        "final_mean_F": round(mean_F_final, 4),
        "initial_mean_H": round(mean_H_initial, 4),
        "final_mean_H": round(mean_H_final, 4),
        "F_reduction_pct": round((initial_F - mean_F_final) / max(initial_F, 1e-9) * 100, 2),
    }
