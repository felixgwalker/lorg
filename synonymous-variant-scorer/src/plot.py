"""Mechanism score bar charts for synonymous variants."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from src.aggregator import ScoredVariant

# Tier colour map (aligned with new tier names: benign/uncertain/likely_functional)
_TIER_COLORS: dict[str, str] = {
    "likely_functional": "#C44E52",  # red
    "uncertain":         "#DD8452",  # orange
    "benign":            "#55A868",  # green
}

# Composite score thresholds matching _assign_tier
_THRESH_UNCERTAIN   = 0.2
_THRESH_FUNCTIONAL  = 0.5


def plot_mechanism_scores(scored: list[ScoredVariant], output_dir: Path) -> tuple[Path, Path]:
    """Render grouped bar chart of per-mechanism scores for each variant.

    Two side-by-side axes are produced:
      - Left: grouped bars for the four mechanism scores per variant.
      - Right: single bar per variant coloured by impact tier, showing the
        composite functional impact index.

    Parameters
    ----------
    scored:
        List of ``ScoredVariant`` instances to plot.
    output_dir:
        Directory where PNG and SVG outputs will be written.

    Returns
    -------
    tuple[Path, Path]
        (png_path, svg_path) of the written figure files.

    Raises
    ------
    ValueError
        If *scored* is empty.
    """
    if not scored:
        raise ValueError("No scored variants to plot.")

    mechanisms = [
        "splicing_score",
        "codon_usage_score",
        "mrna_stability_score",
        "folding_score",
    ]
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

    # ---- Left axis: per-mechanism grouped bars ----
    ax = axes[0]
    for i, (mech_label, color) in enumerate(zip(labels, colors)):
        offset = (i - n_mech / 2 + 0.5) * width
        ax.bar(x + offset, data[:, i], width, label=mech_label, color=color, alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(ids, rotation=45, ha="right", fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Score")
    ax.set_title("Per-Mechanism Scores by Variant")
    ax.legend(loc="upper right", fontsize=7)
    ax.axhline(_THRESH_UNCERTAIN,  color="orange", linestyle="--", linewidth=0.8)
    ax.axhline(_THRESH_FUNCTIONAL, color="red",    linestyle="--", linewidth=0.8)

    # ---- Right axis: composite impact index ----
    ax2 = axes[1]
    composites = [s.composite_score for s in scored]
    bar_colors = [_tier_color(s.impact_tier) for s in scored]
    ax2.bar(x, composites, color=bar_colors, alpha=0.85)
    ax2.set_xticks(x)
    ax2.set_xticklabels(ids, rotation=45, ha="right", fontsize=8)
    ax2.set_ylim(0, 1.05)
    ax2.set_ylabel("Composite Impact Index")
    ax2.set_title("Composite Impact Index by Variant")
    ax2.axhline(_THRESH_UNCERTAIN,  color="orange", linestyle="--", linewidth=0.8)
    ax2.axhline(_THRESH_FUNCTIONAL, color="red",    linestyle="--", linewidth=0.8)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=_TIER_COLORS["benign"],           label="benign"),
        Patch(facecolor=_TIER_COLORS["uncertain"],        label="uncertain"),
        Patch(facecolor=_TIER_COLORS["likely_functional"],label="likely_functional"),
    ]
    ax2.legend(handles=legend_elements, fontsize=7)

    plt.tight_layout()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / "mechanism_scores.png"
    svg_path = output_dir / "mechanism_scores.svg"
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    return png_path, svg_path


def _tier_color(tier: str) -> str:
    """Return hex colour for an impact tier string."""
    return _TIER_COLORS.get(tier, "#888888")
