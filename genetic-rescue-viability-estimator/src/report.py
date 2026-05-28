import json
import os
from datetime import date
import pandas as pd
import numpy as np


# ---------------------------------------------------------------------------
# Text summary report
# ---------------------------------------------------------------------------

def write_summary_report(
    population_data: list,
    sim_result: dict,
    score_result: dict,
    output_dir: str,
    params: dict = None,
) -> str:
    """Write a human-readable plain-text summary report.

    Parameters
    ----------
    population_data : list of per-individual dicts (from vcf_parser / froh pipeline).
    sim_result      : dict returned by run_wf_simulation.
    score_result    : dict returned by compute_viability_score.
    output_dir      : directory to write the report into.
    params          : optional dict of run parameters to include.

    Returns path to the written file.
    """
    path = os.path.join(output_dir, "rescue_summary_report.txt")

    n_indiv = len(population_data)
    froh_values = [d.get("FROH", d.get("F_initial", float("nan"))) for d in population_data]
    mean_froh = float(np.nanmean(froh_values))
    min_froh = float(np.nanmin(froh_values))
    max_froh = float(np.nanmax(froh_values))

    mean_H_initial = float(sim_result["mean_H"][0])
    mean_H_final = float(sim_result["mean_H"][-1])
    mean_F_initial = float(sim_result["mean_F"][0])
    mean_F_final = float(sim_result["mean_F"][-1])
    ci_lo_final = float(np.percentile(sim_result["H_trajectories"][:, -1], 2.5))
    ci_hi_final = float(np.percentile(sim_result["H_trajectories"][:, -1], 97.5))

    n_gens = len(sim_result["generations"]) - 1
    n_reps = sim_result["H_trajectories"].shape[0]

    lines = [
        "=" * 70,
        "  GENETIC RESCUE VIABILITY ESTIMATOR — SUMMARY REPORT",
        "=" * 70,
        f"  Report generated : {date.today().isoformat()}",
        "",
        "POPULATION OVERVIEW",
        "-" * 70,
        f"  Individuals sampled     : {n_indiv}",
        f"  Mean FROH               : {mean_froh:.4f}",
        f"  Min FROH                : {min_froh:.4f}",
        f"  Max FROH                : {max_froh:.4f}",
        "",
    ]

    if params:
        lines += [
            "SIMULATION PARAMETERS",
            "-" * 70,
        ]
        for k, v in params.items():
            lines.append(f"  {k:<28}: {v}")
        lines.append("")

    lines += [
        "SIMULATION RESULTS",
        "-" * 70,
        f"  Generations simulated   : {n_gens}",
        f"  Replicates              : {n_reps}",
        f"  Initial mean H          : {mean_H_initial:.4f}",
        f"  Final mean H (rescue)   : {mean_H_final:.4f}",
        f"  Initial mean F          : {mean_F_initial:.4f}",
        f"  Final mean F (rescue)   : {mean_F_final:.4f}",
        f"  Final H 95% CI          : [{ci_lo_final:.4f}, {ci_hi_final:.4f}]",
        "",
        "VIABILITY SCORE",
        "-" * 70,
        f"  Score                   : {score_result['score']:.1f} / 100",
        f"  Grade                   : {score_result['grade']}",
        f"  F reduction             : {score_result['F_reduction_pct']:.1f}%",
        f"  Interpretation          : {score_result['interpretation']}",
        "",
        "INBREEDING CLASSIFICATION (FROH thresholds)",
        "-" * 70,
        "  < 0.05   : Low inbreeding",
        "  0.05-0.125: Moderate inbreeding",
        "  0.125-0.25: High inbreeding",
        "  > 0.25   : Very high inbreeding",
        "",
    ]

    # Recommendation
    if score_result["grade"] in ("A", "B"):
        rec = "RECOMMENDED — genetic rescue is expected to substantially improve viability."
    elif score_result["grade"] == "C":
        rec = "CONSIDER — genetic rescue may offer moderate benefit; further modelling advised."
    else:
        rec = "NOT RECOMMENDED — expected benefit is limited given current parameters."

    lines += [
        "RECOMMENDATION",
        "-" * 70,
        f"  {rec}",
        "",
        "=" * 70,
    ]

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    return path


# ---------------------------------------------------------------------------
# CSV / JSON outputs (unchanged from original)
# ---------------------------------------------------------------------------

def write_froh_csv(population_data, output_dir):
    path = os.path.join(output_dir, "froh_per_individual.csv")
    rows = []
    for indiv in population_data:
        rows.append({
            "individual": indiv["individual"],
            "F_initial": indiv["F_initial"],
            "H_initial": indiv["H_initial"],
        })
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)
    return path


def write_trajectories_csv(sim_result, output_dir):
    path = os.path.join(output_dir, "rescue_trajectories.csv")
    gens = sim_result["generations"]
    df = pd.DataFrame({
        "generation": gens,
        "mean_H": np.round(sim_result["mean_H"], 6),
        "ci_lower": np.round(sim_result["ci_lower"], 6),
        "ci_upper": np.round(sim_result["ci_upper"], 6),
        "mean_F": np.round(sim_result["mean_F"], 6),
    })
    df.to_csv(path, index=False)
    return path


def write_viability_json(score_result, output_dir):
    path = os.path.join(output_dir, "viability_score.json")
    with open(path, "w") as fh:
        json.dump(score_result, fh, indent=2)
    return path
