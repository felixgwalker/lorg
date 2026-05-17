import json
import os
import pandas as pd
import numpy as np


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
