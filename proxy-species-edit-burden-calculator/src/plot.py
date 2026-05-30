"""Plot generation for proxy-species edit-burden calculator.

Produces:
  - edit_density_plot.<fmt>  — bar chart (counts by variant class) +
                               positional density plot along chromosomes
  - impact_pie_plot.<fmt>    — pie chart of variants by impact category
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# Colour palette for impact categories
# ---------------------------------------------------------------------------
_IMPACT_COLOURS = {
    "coding_nonsense": "#d62728",
    "coding_nonsynonymous": "#ff7f0e",
    "exonic_nonsynonymous": "#ff7f0e",   # legacy alias
    "splice_site": "#9467bd",
    "regulatory": "#17becf",
    "coding_synonymous": "#2ca02c",
    "exonic_synonymous": "#2ca02c",      # legacy alias
    "intronic": "#8c564b",
    "intergenic": "#7f7f7f",
}


def _get_impact(v) -> str:
    """Return impact_category from either a dict or a Variant dataclass."""
    if hasattr(v, "impact_category"):
        return v.impact_category or "intergenic"
    return v.get("impact_category", "intergenic")


def _get_pos(v):
    if hasattr(v, "position"):
        return v.position
    return v.get("pos", 0)


def _get_chrom(v) -> str:
    if hasattr(v, "chrom"):
        return v.chrom or "chrUn"
    return v.get("chrom", "chrUn")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def plot_edit_burden(variants, burden, chrom_lengths, output_dir, fmt="png"):
    """Create the main edit-burden figure (variant class bar + density line).

    Parameters
    ----------
    variants : list
        Annotated variant list.
    burden : dict
        Output of compute_burden().
    chrom_lengths : dict
        {chrom: length_bp}
    output_dir : str
    fmt : str
        "png" or "svg".

    Returns
    -------
    str
        Path to the saved figure.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    class_counts = burden["class_counts"]
    classes = list(class_counts.keys())
    counts = [class_counts[c] for c in classes]
    weights = {
        "SNV": 1, "SMALL_INS": 3, "SMALL_DEL": 3,
        "LARGE_INS": 10, "LARGE_DEL": 10, "SV_INS": 50, "SV_DEL": 50,
    }
    weighted_counts = [class_counts[c] * weights.get(c, 1) for c in classes]

    x = np.arange(len(classes))
    width = 0.35
    bars1 = ax1.bar(x - width / 2, counts, width, label="Raw count",
                    color="steelblue", alpha=0.8)
    ax1.bar(x + width / 2, weighted_counts, width, label="Weighted burden",
            color="coral", alpha=0.8)
    ax1.set_xticks(x)
    ax1.set_xticklabels(classes, rotation=30, ha="right")
    ax1.set_ylabel("Count / Burden")
    ax1.set_title("Edit Counts by Variant Class")
    ax1.legend()

    for bar in bars1:
        h = bar.get_height()
        if h > 0:
            ax1.text(
                bar.get_x() + bar.get_width() / 2, h, str(int(h)),
                ha="center", va="bottom", fontsize=7,
            )

    chroms = list(chrom_lengths.keys())
    if not chroms:
        ax2.text(0.5, 0.5, "No chromosome data", ha="center", va="center",
                 transform=ax2.transAxes)
    else:
        max_density = 0.0
        for ci, chrom in enumerate(chroms):
            clen = chrom_lengths[chrom]
            chrom_vars = [v for v in variants if _get_chrom(v) == chrom]
            if not chrom_vars:
                continue
            positions = [_get_pos(v) for v in chrom_vars]
            n_bins = max(10, clen // 500)
            counts_hist, bin_edges = np.histogram(
                positions, bins=n_bins, range=(0, clen)
            )
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            density = counts_hist / max(clen / 1000, 1)
            local_max = float(density.max()) if density.max() > 0 else 1.0
            ax2.plot(
                bin_centers,
                density + ci * local_max * 1.5,
                label=chrom,
                linewidth=1.5,
            )
            max_density = max(max_density, local_max)

        ax2.set_xlabel("Genomic position (bp)")
        ax2.set_ylabel("Variant density (per kb, offset by chrom)")
        ax2.set_title("Edit Density Along Chromosomes")
        ax2.legend(loc="upper right")

    fig.suptitle(
        f"Edit Burden: {burden['total_edits']} variants, "
        f"weighted={burden['weighted_burden']}, "
        f"{burden['normalized_burden_per_mb']:.1f}/Mb",
        fontsize=11,
    )
    fig.tight_layout()

    path = os.path.join(output_dir, f"edit_density_plot.{fmt}")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_impact_pie(variants, output_dir, fmt="png"):
    """Create a pie chart of variant impact categories.

    Parameters
    ----------
    variants : list
        Annotated variant list.
    output_dir : str
    fmt : str
        "png" or "svg".

    Returns
    -------
    str
        Path to the saved figure.
    """
    # Tally by impact category.
    counts: dict[str, int] = {}
    for v in variants:
        cat = _get_impact(v)
        counts[cat] = counts.get(cat, 0) + 1

    if not counts:
        fig, ax = plt.subplots(figsize=(7, 7))
        ax.text(0.5, 0.5, "No variants", ha="center", va="center",
                transform=ax.transAxes, fontsize=14)
        ax.axis("off")
    else:
        labels = list(counts.keys())
        sizes = [counts[k] for k in labels]
        colours = [_IMPACT_COLOURS.get(k, "#aec7e8") for k in labels]

        fig, ax = plt.subplots(figsize=(8, 7))
        wedges, texts, autotexts = ax.pie(
            sizes,
            labels=None,
            autopct=lambda pct: f"{pct:.1f}%\n({int(round(pct / 100 * sum(sizes)))})",
            colors=colours,
            startangle=140,
            pctdistance=0.75,
        )
        for at in autotexts:
            at.set_fontsize(8)

        ax.legend(
            wedges,
            [f"{lbl} ({counts[lbl]})" for lbl in labels],
            loc="lower left",
            fontsize=9,
            bbox_to_anchor=(-0.05, -0.05),
        )
        ax.set_title("Variant Impact Category Distribution", fontsize=13, pad=16)

    fig.tight_layout()
    path = os.path.join(output_dir, f"impact_pie_plot.{fmt}")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
