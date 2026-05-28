import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

def _save(fig, output_dir: str, stem: str, fmt: str) -> list:
    """Save figure as PNG and/or SVG.  Always saves PNG; also SVG when requested.

    Returns list of saved file paths.
    """
    paths = []
    for ext in (["png", "svg"] if fmt == "svg" else ["png"]):
        path = os.path.join(output_dir, f"{stem}.{ext}")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        paths.append(path)
    return paths


# ---------------------------------------------------------------------------
# FROH comparison bar chart
# ---------------------------------------------------------------------------

def plot_froh_comparison(
    froh_before: list,
    froh_after: list,
    individual_ids: list = None,
    output_dir: str = ".",
    fmt: str = "png",
) -> list:
    """Bar chart comparing per-individual FROH before and after rescue.

    Parameters
    ----------
    froh_before     : list of FROH values before rescue (length n_individuals).
    froh_after      : list of FROH values after rescue  (length n_individuals).
    individual_ids  : optional labels; defaults to IND_001 … IND_n.
    output_dir      : directory for output file(s).
    fmt             : "png" (PNG only) or "svg" (PNG + SVG).

    Returns list of saved file paths.
    """
    n = len(froh_before)
    if individual_ids is None:
        individual_ids = [f"IND_{i+1:03d}" for i in range(n)]

    x = np.arange(n)
    bar_width = 0.4

    fig, ax = plt.subplots(figsize=(max(10, n * 0.6), 6))

    bars_before = ax.bar(x - bar_width / 2, froh_before, bar_width,
                         label="FROH before rescue", color="tab:red", alpha=0.8)
    bars_after = ax.bar(x + bar_width / 2, froh_after, bar_width,
                        label="FROH after rescue", color="tab:blue", alpha=0.8)

    # Threshold reference lines
    ax.axhline(y=0.25, color="darkred",   linestyle="--", linewidth=1, label="Very high (0.25)")
    ax.axhline(y=0.125, color="orange",   linestyle="--", linewidth=1, label="High (0.125)")
    ax.axhline(y=0.05,  color="green",    linestyle="--", linewidth=1, label="Low (0.05)")

    if n <= 30:
        ax.set_xticks(x)
        ax.set_xticklabels(individual_ids, rotation=45, ha="right", fontsize=8)
    else:
        ax.set_xticks([])
        ax.set_xlabel("Individuals")

    ax.set_ylabel("FROH")
    ax.set_title("Per-individual FROH Before vs After Genetic Rescue")
    ax.set_ylim(0, min(1.05, max(max(froh_before), max(froh_after)) * 1.15 + 0.05))
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()

    paths = _save(fig, output_dir, "froh_comparison", fmt)
    plt.close(fig)
    return paths


# ---------------------------------------------------------------------------
# WF simulation trajectory plot
# ---------------------------------------------------------------------------

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

    paths = _save(fig, output_dir, "trajectory_plot", fmt)
    plt.close(fig)
    # Return single path string for backwards compatibility with pipeline
    return paths[0]


# ---------------------------------------------------------------------------
# Final heterozygosity distribution (violin) plot
# ---------------------------------------------------------------------------

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

    paths = _save(fig, output_dir, "distribution_plot", fmt)
    plt.close(fig)
    return paths[0]


# ---------------------------------------------------------------------------
# WF allele-frequency trajectory plot (uses simulate_wright_fisher output)
# ---------------------------------------------------------------------------

def plot_wf_trajectories(
    trajectories: np.ndarray,
    output_dir: str = ".",
    fmt: str = "png",
    title: str = "Wright-Fisher Allele Frequency Trajectories",
    max_shown: int = 50,
) -> list:
    """Plot allele-frequency trajectories from simulate_wright_fisher.

    Parameters
    ----------
    trajectories : np.ndarray shape (n_reps, n_gen).
    output_dir   : directory for output file(s).
    fmt          : "png" or "svg".
    title        : plot title.
    max_shown    : maximum number of individual trajectories to draw (for clarity).

    Returns list of saved file paths.
    """
    n_reps, n_gen = trajectories.shape
    gens = np.arange(1, n_gen + 1)

    fig, ax = plt.subplots(figsize=(10, 6))

    shown = min(n_reps, max_shown)
    for i in range(shown):
        ax.plot(gens, trajectories[i], color="steelblue", alpha=0.15, linewidth=0.8)

    mean_freq = trajectories.mean(axis=0)
    ci_lo = np.percentile(trajectories, 2.5, axis=0)
    ci_hi = np.percentile(trajectories, 97.5, axis=0)

    ax.fill_between(gens, ci_lo, ci_hi, alpha=0.3, color="steelblue", label="95% CI")
    ax.plot(gens, mean_freq, color="navy", linewidth=2, label="Mean frequency")

    ax.axhline(y=0.0, color="black", linewidth=0.5, linestyle="--")
    ax.axhline(y=1.0, color="black", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Allele Frequency")
    ax.set_title(title)
    ax.set_ylim(-0.02, 1.02)
    ax.legend(loc="best")
    fig.tight_layout()

    paths = _save(fig, output_dir, "wf_allele_freq", fmt)
    plt.close(fig)
    return paths
