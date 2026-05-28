"""Plots for off-target sites: Manhattan plot and mismatch-position heatmap."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from pathlib import Path


MISMATCH_COLORS = {0: "#2196F3", 1: "#4CAF50", 2: "#FF9800", 3: "#F44336"}
KMER_LEN = 20


def plot_manhattan(hits: list[dict], output_dir: Path, fmt: str = "png") -> Path:
    """Manhattan-style scatter plot: genomic position vs CFD score.

    Points are coloured by number of mismatches.
    Returns the path to the saved figure.
    """
    out_path = output_dir / f"off_target_manhattan.{fmt}"

    fig, ax = plt.subplots(figsize=(12, 5))

    if not hits:
        ax.text(0.5, 0.5, "No off-target sites found", ha="center", va="center",
                transform=ax.transAxes, fontsize=14)
        ax.set_xlabel("Genome position")
        ax.set_ylabel("CFD score")
        ax.set_title("Off-target site Manhattan plot")
        fig.tight_layout()
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        # Also generate the heatmap (empty)
        plot_mismatch_heatmap(hits, output_dir, fmt=fmt)
        return out_path

    df = pd.DataFrame(hits)
    df["mismatches"] = df["mismatches"].astype(int)

    for mm in sorted(df["mismatches"].unique()):
        sub = df[df["mismatches"] == mm]
        color = MISMATCH_COLORS.get(mm, "#9E9E9E")
        size = max(20, 80 - mm * 20)
        ax.scatter(sub["pos"], sub["cfd_score"], c=color, s=size,
                   alpha=0.75, linewidths=0.5, edgecolors="white",
                   label=f"{mm} mismatch{'es' if mm != 1 else ''}", zorder=3)

    ax.set_xlabel("Genome position (bp)", fontsize=11)
    ax.set_ylabel("CFD score", fontsize=11)
    ax.set_ylim(-0.05, 1.1)
    ax.set_title("CRISPR Off-target Sites — Manhattan Plot", fontsize=13, fontweight="bold")
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8, alpha=0.6, label="CFD=0.5 threshold")
    ax.grid(axis="y", alpha=0.3)

    legend_patches = [
        mpatches.Patch(color=MISMATCH_COLORS[mm], label=f"{mm} mismatch{'es' if mm != 1 else ''}")
        for mm in sorted(MISMATCH_COLORS)
        if mm in df["mismatches"].values
    ]
    legend_patches.append(
        plt.Line2D([0], [0], color="gray", linestyle="--", linewidth=0.8, label="CFD=0.5 threshold")
    )
    ax.legend(handles=legend_patches, loc="upper right", fontsize=9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    # Always generate the heatmap alongside the Manhattan plot
    plot_mismatch_heatmap(hits, output_dir, fmt=fmt)

    return out_path


def plot_mismatch_heatmap(hits: list[dict], output_dir: Path, fmt: str = "png") -> Path:
    """Heatmap showing which guide positions are most frequently mismatched.

    Rows = individual guides; columns = guide positions 1-20 (1 = PAM-distal,
    20 = PAM-proximal).  Cell value = number of off-target sites with a
    mismatch at that position.

    Returns the path to the saved figure.
    """
    out_path = output_dir / f"mismatch_position_heatmap.{fmt}"

    # Filter to only off-target hits (mismatches > 0)
    off_hits = [h for h in hits if h.get("mismatches", 0) > 0]

    if not off_hits:
        fig, ax = plt.subplots(figsize=(10, 3))
        ax.text(0.5, 0.5, "No off-target sites to display", ha="center", va="center",
                transform=ax.transAxes, fontsize=13)
        ax.set_title("Mismatch Position Heatmap")
        fig.tight_layout()
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        return out_path

    # Collect unique guide names (preserving order)
    seen: list[str] = []
    for h in off_hits:
        gn = h["guide_name"]
        if gn not in seen:
            seen.append(gn)
    guide_names = seen

    # Build count matrix: shape (n_guides, KMER_LEN)
    matrix = np.zeros((len(guide_names), KMER_LEN), dtype=int)
    guide_idx = {gn: i for i, gn in enumerate(guide_names)}

    for h in off_hits:
        gn = h["guide_name"]
        mm_pos = h.get("mismatch_positions", [])
        # mismatch_positions may be a list of ints, or a semicolon-separated string
        if isinstance(mm_pos, str):
            mm_pos = [int(x) for x in mm_pos.split(";") if x.strip()]
        for pos in mm_pos:
            col = int(pos) - 1  # convert 1-based to 0-based
            if 0 <= col < KMER_LEN:
                matrix[guide_idx[gn], col] += 1

    n_guides = len(guide_names)
    fig_height = max(3, 1.0 * n_guides + 1.5)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    im = ax.imshow(matrix, aspect="auto", cmap="YlOrRd", interpolation="nearest")
    plt.colorbar(im, ax=ax, label="Off-target count with mismatch at position")

    ax.set_xticks(range(KMER_LEN))
    ax.set_xticklabels([str(i + 1) for i in range(KMER_LEN)], fontsize=8)
    ax.set_yticks(range(n_guides))
    ax.set_yticklabels(guide_names, fontsize=9)

    ax.set_xlabel("Guide position (1 = PAM-distal, 20 = PAM-proximal)", fontsize=10)
    ax.set_ylabel("Guide", fontsize=10)
    ax.set_title("Mismatch Position Frequency Heatmap", fontsize=13, fontweight="bold")

    # Annotate cells with counts > 0
    for r in range(n_guides):
        for c in range(KMER_LEN):
            val = matrix[r, c]
            if val > 0:
                ax.text(c, r, str(val), ha="center", va="center",
                        fontsize=7, color="black" if val < matrix.max() * 0.7 else "white")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path
