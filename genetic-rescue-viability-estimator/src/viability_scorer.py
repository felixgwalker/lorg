import numpy as np


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
