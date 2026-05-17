"""Mechanism score bar charts for synonymous variants."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from src.aggregator import ScoredVariant


def plot_mechanism_scores(scored: list[ScoredVariant], output_dir: Path) -> tuple[Path, Path]:
    """Render grouped bar chart of per-mechanism scores for each variant."""
    if not scored:
        raise ValueError("No scored variants to plot.")

    mechanisms = ["splicing_score", "codon_usage_score", "mrna_stability_score", "folding_score"]
    labels = ["Splicing", "Codon Usage", "mRNA Stability", "Folding"]
    ids = [s.variant_id for s in scored]
    n_vars = len(scored)
    n_mech = len(mechanisms)

    data = np.array([
        [getattr(s, m) for m in mechanisms]
        for s in scored
    ])

    fig, axes = plt.subplots(1, 2, figsize=(max(10, n_vars * 0.8 + 4), 6))

    x = np.arange(n_vars)
    width = 0.18
    colors = ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]

    ax = axes[0]
    for i, (mech, color) in enumerate(zip(labels, colors)):
        offset = (i - n_mech / 2 + 0.5) * width
        ax.bar(x + offset, data[:, i], width, label=mech, color=color, alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(ids, rotation=45, ha="right", fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Score")
    ax.set_title("Per-Mechanism Scores by Variant")
    ax.legend(loc="upper right", fontsize=7)
    ax.axhline(0.3, color="orange", linestyle="--", linewidth=0.8, label="MODERATE")
    ax.axhline(0.6, color="red", linestyle="--", linewidth=0.8, label="HIGH")

    ax2 = axes[1]
    composites = [s.composite_score for s in scored]
    bar_colors = [_tier_color(s.impact_tier) for s in scored]
    ax2.bar(x, composites, color=bar_colors, alpha=0.85)
    ax2.set_xticks(x)
    ax2.set_xticklabels(ids, rotation=45, ha="right", fontsize=8)
    ax2.set_ylim(0, 1.05)
    ax2.set_ylabel("Composite Impact Index")
    ax2.set_title("Composite Impact Index by Variant")
    ax2.axhline(0.3, color="orange", linestyle="--", linewidth=0.8)
    ax2.axhline(0.6, color="red", linestyle="--", linewidth=0.8)
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#55A868", label="LOW"),
        Patch(facecolor="#DD8452", label="MODERATE"),
        Patch(facecolor="#C44E52", label="HIGH"),
    ]
    ax2.legend(handles=legend_elements, fontsize=7)

    plt.tight_layout()
    png_path = output_dir / "mechanism_scores.png"
    svg_path = output_dir / "mechanism_scores.svg"
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    return png_path, svg_path


def _tier_color(tier: str) -> str:
    return {"HIGH": "#C44E52", "MODERATE": "#DD8452", "LOW": "#55A868"}.get(tier, "#888888")
