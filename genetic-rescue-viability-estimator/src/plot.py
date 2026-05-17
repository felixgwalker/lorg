import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os


def plot_trajectories(sim_result, baseline_result, output_dir, fmt="png"):
    gens = sim_result["generations"]
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.fill_between(gens, baseline_result["ci_lower"], baseline_result["ci_upper"],
                    alpha=0.2, color="tab:red", label="No-rescue 95% CI")
    ax.plot(gens, baseline_result["mean_H"], color="tab:red", linewidth=2,
            linestyle="--", label="No-rescue mean H")

    ax.fill_between(gens, sim_result["ci_lower"], sim_result["ci_upper"],
                    alpha=0.2, color="tab:blue", label="Rescue 95% CI")
    ax.plot(gens, sim_result["mean_H"], color="tab:blue", linewidth=2,
            label="Rescue mean H")

    ax.axvline(x=1, color="green", linestyle=":", linewidth=1.5, label="Rescue event")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Mean Heterozygosity (H)")
    ax.set_title("Wright-Fisher Simulation: Heterozygosity Trajectories")
    ax.legend(loc="best")
    ax.set_ylim(0, 1)
    fig.tight_layout()

    path = os.path.join(output_dir, f"trajectory_plot.{fmt}")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_distribution(sim_result, baseline_result, output_dir, fmt="png"):
    H_rescue_final = sim_result["H_trajectories"][:, -1]
    H_baseline_final = baseline_result["H_trajectories"][:, -1]

    fig, ax = plt.subplots(figsize=(8, 6))

    positions = [1, 2]
    data = [H_baseline_final, H_rescue_final]
    parts = ax.violinplot(data, positions=positions, showmedians=True, showextrema=True)

    for pc in parts["bodies"]:
        pc.set_alpha(0.7)
    parts["bodies"][0].set_facecolor("tab:red")
    parts["bodies"][1].set_facecolor("tab:blue")

    ax.set_xticks(positions)
    ax.set_xticklabels(["No Rescue", "With Rescue"])
    ax.set_ylabel("Final Heterozygosity (H)")
    ax.set_title("Distribution of Final Heterozygosity Across Replicates")
    ax.set_ylim(0, 1)
    fig.tight_layout()

    path = os.path.join(output_dir, f"distribution_plot.{fmt}")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
