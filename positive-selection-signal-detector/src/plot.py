import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os


def plot_dnds(tested_genes, output_dir, fmt="png"):
    omegas = [g["omega"] for g in tested_genes]
    significant = [g.get("significant", False) for g in tested_genes]
    dN_vals = [g["mean_dN"] for g in tested_genes]
    dS_vals = [g["mean_dS"] for g in tested_genes]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    n_bins = max(5, len(omegas) // 3)
    ax1.hist(omegas, bins=n_bins, color="steelblue", edgecolor="white", alpha=0.8)
    ax1.axvline(x=1.0, color="red", linestyle="--", linewidth=1.5, label="omega=1 (neutral)")
    ax1.axvline(x=1.5, color="orange", linestyle=":", linewidth=1.5, label="omega=1.5 threshold")
    ax1.set_xlabel("omega (dN/dS)")
    ax1.set_ylabel("Number of genes")
    ax1.set_title("Distribution of dN/dS Ratios")
    ax1.legend()

    colors = ["tab:red" if sig else "tab:blue" for sig in significant]
    ax2.scatter(dS_vals, dN_vals, c=colors, alpha=0.7, s=60, edgecolors="none")
    max_val = max(max(dN_vals + [0.001]), max(dS_vals + [0.001]))
    ax2.plot([0, max_val], [0, max_val], "k--", linewidth=1, label="dN=dS (neutral)")
    ax2.set_xlabel("dS (synonymous substitutions/site)")
    ax2.set_ylabel("dN (nonsynonymous substitutions/site)")
    ax2.set_title("dN vs dS Scatter Plot")
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="tab:red", markersize=8, label="Significant"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="tab:blue", markersize=8, label="Not significant"),
        Line2D([0], [0], color="k", linestyle="--", label="dN=dS"),
    ]
    ax2.legend(handles=legend_elements)

    fig.tight_layout()
    path = os.path.join(output_dir, f"dnds_plot.{fmt}")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
