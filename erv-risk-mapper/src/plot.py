"""Chromosome density bar chart for ERV risk tiers."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from pathlib import Path
from collections import defaultdict


TIER_COLORS = {"low": "#4CAF50", "moderate": "#FF9800", "high": "#F44336"}
TIER_ORDER = ["low", "moderate", "high"]


def plot_chromosome_density(hits: list[dict], output_dir: Path, fmt: str = "png") -> Path:
    out_path = output_dir / f"erv_chromosome_density.{fmt}"

    fig, ax = plt.subplots(figsize=(10, 6))

    if not hits:
        ax.text(0.5, 0.5, "No ERV elements detected", ha="center", va="center",
                transform=ax.transAxes, fontsize=14)
        ax.set_title("ERV Chromosome Density")
        fig.tight_layout()
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        return out_path

    df = pd.DataFrame(hits)
    chroms = sorted(df["chrom"].unique(), key=lambda c: (
        int(c.replace("chr", "")) if c.replace("chr", "").isdigit() else 999
    ))

    counts: dict[str, dict[str, int]] = {c: defaultdict(int) for c in chroms}
    for _, row in df.iterrows():
        counts[row["chrom"]][row["risk_tier"]] += 1

    x = range(len(chroms))
    width = 0.25
    offsets = [-width, 0, width]

    for tier, offset in zip(TIER_ORDER, offsets):
        values = [counts[c].get(tier, 0) for c in chroms]
        bars = ax.bar(
            [xi + offset for xi in x], values, width=width,
            label=f"{tier.capitalize()} risk",
            color=TIER_COLORS[tier], alpha=0.85, edgecolor="white", linewidth=0.8,
        )
        for bar, val in zip(bars, values):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.05,
                        str(val), ha="center", va="bottom", fontsize=8, fontweight="bold")

    ax.set_xticks(list(x))
    ax.set_xticklabels(chroms, fontsize=10)
    ax.set_xlabel("Chromosome", fontsize=11)
    ax.set_ylabel("Number of ERV elements", fontsize=11)
    ax.set_title("ERV Element Density by Chromosome and Risk Tier", fontsize=13, fontweight="bold")
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.grid(axis="y", alpha=0.3)

    legend_patches = [
        mpatches.Patch(color=TIER_COLORS[t], label=f"{t.capitalize()} risk") for t in TIER_ORDER
    ]
    ax.legend(handles=legend_patches, loc="upper right", fontsize=10)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path
