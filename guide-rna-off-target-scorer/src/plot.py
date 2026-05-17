"""Manhattan-style plot for off-target sites."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from pathlib import Path


MISMATCH_COLORS = {0: "#2196F3", 1: "#4CAF50", 2: "#FF9800", 3: "#F44336"}


def plot_manhattan(hits: list[dict], output_dir: Path, fmt: str = "png") -> Path:
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
    return out_path
