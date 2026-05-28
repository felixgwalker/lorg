import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os


def plot_omega_scan(windows, output_dir, alpha=0.05, fmt="png"):
    """Generate a sliding-window omega scan plot.

    Parameters
    ----------
    windows : list of dict
        Each dict must contain:
            'window_pos'   : int — alignment position of the window centre/start
            'omega'        : float — dN/dS estimate for the window
        Optionally:
            'p_adj'        : float — FDR-adjusted p-value (used to highlight
                             significant windows).  If absent, 'p_raw' is used.
            'p_raw'        : float — raw LRT p-value
    output_dir : str
        Directory where the PNG/SVG will be written.
    alpha : float
        Significance threshold for highlighting windows (default 0.05).
    fmt : str
        'png' or 'svg' (default 'png').

    Returns
    -------
    str
        Path to the saved figure.
    """
    os.makedirs(output_dir, exist_ok=True)

    positions = [w["window_pos"] for w in windows]
    omegas = [min(float(w["omega"]) if not np.isinf(w["omega"]) else 10.0, 10.0)
              for w in windows]

    # Determine significance for each window
    sig_flags = []
    for w in windows:
        pval = w.get("p_adj", w.get("p_raw", 1.0))
        sig_flags.append(float(pval) < alpha and w.get("omega", 1.0) > 1.0)

    fig, ax = plt.subplots(figsize=(12, 5))

    # Line plot of omega
    ax.plot(positions, omegas, color="steelblue", linewidth=1.5, zorder=2, label="omega (dN/dS)")

    # Highlight significant windows
    sig_positions = [p for p, s in zip(positions, sig_flags) if s]
    sig_omegas = [o for o, s in zip(omegas, sig_flags) if s]
    if sig_positions:
        ax.scatter(sig_positions, sig_omegas, color="crimson", zorder=3,
                   s=50, label=f"Significant (FDR < {alpha})", edgecolors="none")

    # Neutral reference line at omega = 1
    ax.axhline(y=1.0, color="black", linestyle="--", linewidth=1.2, zorder=1,
               label="omega = 1 (neutral)")

    ax.set_xlabel("Alignment position (nt)")
    ax.set_ylabel("omega (dN/dS)")
    ax.set_title("Sliding-window dN/dS scan")
    ax.legend(loc="upper right")
    ax.set_ylim(bottom=0)

    fig.tight_layout()
    path = os.path.join(output_dir, f"omega_scan.{fmt}")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


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
