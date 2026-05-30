"""ROH distribution histogram, FROH barplot, Ne trajectory, karyotype, and
FROH distribution plots."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

from pathlib import Path
from typing import Any

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from src.froh_calculator import FROHResult
from src.roh_detector import ROHSegment

# Length-class colour palette (shared across all plotting functions)
_CLASS_COLORS: dict[str, str] = {
    "short":  "#4C72B0",   # blue  — ancient
    "medium": "#DD8452",   # orange — moderate
    "long":   "#C44E52",   # red   — recent
}
_CLASS_LABELS: dict[str, str] = {
    "short":  "Short (<100 kb)",
    "medium": "Medium (100 kb–1 Mb)",
    "long":   "Long (>1 Mb)",
}


def plot_roh(
    segments: list[ROHSegment],
    froh_results: list[FROHResult],
    ne_estimates: list[Any],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Render ROH histogram, FROH barplot, and Ne trajectory; save PNG and SVG."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    _plot_roh_histogram(axes[0], segments)
    _plot_froh_barplot(axes[1], froh_results)
    _plot_ne_trajectory(axes[2], ne_estimates)

    plt.tight_layout()
    png_path = output_dir / "roh_analysis.png"
    svg_path = output_dir / "roh_analysis.svg"
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    return png_path, svg_path


# ---------------------------------------------------------------------------
# Public karyotype and distribution plots (per spec)
# ---------------------------------------------------------------------------

def plot_roh_karyotype(
    roh_segments: list[ROHSegment],
    genome_size: dict[str, int] | int,
    output_dir: Path,
) -> tuple[Path, Path]:
    """Horizontal bar chart showing ROH as filled rectangles per chromosome.

    Each chromosome is drawn as a thin grey background bar spanning
    [0, chrom_length].  ROH segments are overlaid as coloured rectangles
    (blue = short, orange = medium, red = long).

    Parameters
    ----------
    roh_segments:
        ROH segments for one or more individuals.  If multiple individuals
        are present the ROH from all individuals are overlaid (useful for a
        population-level view).
    genome_size:
        Either a dict mapping chromosome names to their lengths (bp), or a
        single integer used as the length of every chromosome present in the
        segments.
    output_dir:
        Directory for output files.  Created if absent.

    Returns
    -------
    (png_path, svg_path) — absolute Path objects for the saved files.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not roh_segments:
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.text(0.5, 0.5, "No ROH segments to display",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_title("ROH Karyotype")
        png_path = output_dir / "roh_karyotype.png"
        svg_path = output_dir / "roh_karyotype.svg"
        fig.savefig(png_path, dpi=150, bbox_inches="tight")
        fig.savefig(svg_path, bbox_inches="tight")
        plt.close(fig)
        return png_path, svg_path

    # Determine chromosomes and their lengths
    chroms: list[str] = sorted(
        {s.chrom for s in roh_segments},
        key=_chrom_sort_key,
    )

    def chrom_length(c: str) -> int:
        if isinstance(genome_size, dict):
            return genome_size.get(c, 0)
        if isinstance(genome_size, int):
            return genome_size
        return max(
            (s.end_pos for s in roh_segments if s.chrom == c),
            default=1,
        )

    n_chroms = len(chroms)
    fig_height = max(4, 0.5 * n_chroms)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    bar_height = 0.6
    for y_idx, chrom in enumerate(chroms):
        clen = chrom_length(chrom)
        # Grey background bar representing the full chromosome
        ax.barh(
            y_idx, clen, left=0, height=bar_height,
            color="#CCCCCC", edgecolor="#999999", linewidth=0.4, zorder=1,
        )
        # ROH rectangles
        chrom_segs = [s for s in roh_segments if s.chrom == chrom]
        for seg in chrom_segs:
            color = _CLASS_COLORS.get(seg.length_class, "#888888")
            ax.barh(
                y_idx,
                seg.length_bp,
                left=seg.start_pos,
                height=bar_height * 0.8,
                color=color,
                alpha=0.85,
                zorder=2,
            )

    ax.set_yticks(range(n_chroms))
    ax.set_yticklabels(chroms, fontsize=9)
    ax.set_xlabel("Genomic position (bp)", fontsize=10)
    ax.set_title("ROH Karyotype", fontsize=12)

    # Legend
    legend_patches = [
        mpatches.Patch(color=_CLASS_COLORS[cls], label=_CLASS_LABELS[cls])
        for cls in ("short", "medium", "long")
    ]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=8)

    plt.tight_layout()
    png_path = output_dir / "roh_karyotype.png"
    svg_path = output_dir / "roh_karyotype.svg"
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    return png_path, svg_path


def plot_froh_distribution(
    froh_values: list[FROHResult],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Plot FROH distribution across multiple individuals.

    Produces a stacked bar chart showing the short / medium / long
    contributions to total FROH for each individual, plus a histogram of
    overall FROH values.

    Parameters
    ----------
    froh_values:
        List of :class:`~src.froh_calculator.FROHResult` objects, one per
        individual.
    output_dir:
        Directory for output files.

    Returns
    -------
    (png_path, svg_path)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not froh_values:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No FROH data", ha="center", va="center",
                transform=ax.transAxes)
        ax.set_title("FROH Distribution")
        png_path = output_dir / "froh_distribution.png"
        svg_path = output_dir / "froh_distribution.svg"
        fig.savefig(png_path, dpi=150, bbox_inches="tight")
        fig.savefig(svg_path, bbox_inches="tight")
        plt.close(fig)
        return png_path, svg_path

    ids = [fr.individual_id for fr in froh_values]
    froh_short = np.array([fr.froh_short for fr in froh_values])
    froh_medium = np.array([fr.froh_medium for fr in froh_values])
    froh_long = np.array([fr.froh_long for fr in froh_values])
    froh_total = np.array([fr.froh for fr in froh_values])

    x = np.arange(len(ids))
    width = 0.55

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # --- Left panel: stacked bar chart per individual ---
    ax = axes[0]
    p1 = ax.bar(x, froh_short, width, label=_CLASS_LABELS["short"],
                color=_CLASS_COLORS["short"], alpha=0.85)
    p2 = ax.bar(x, froh_medium, width, bottom=froh_short,
                label=_CLASS_LABELS["medium"],
                color=_CLASS_COLORS["medium"], alpha=0.85)
    p3 = ax.bar(x, froh_long, width, bottom=froh_short + froh_medium,
                label=_CLASS_LABELS["long"],
                color=_CLASS_COLORS["long"], alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(ids, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("F_ROH")
    ax.set_title("FROH per Individual (by length class)")
    ax.legend(fontsize=8)
    ax.set_ylim(0, min(1.0, max(froh_total) * 1.35 + 0.005))

    # Annotate total FROH above each bar
    for xi, v in zip(x, froh_total):
        ax.text(xi, v + 0.001, f"{v:.3f}", ha="center", va="bottom", fontsize=7)

    # --- Right panel: FROH histogram ---
    ax2 = axes[1]
    if len(froh_total) > 1:
        bins = min(10, len(froh_total))
        ax2.hist(froh_total, bins=bins, color="#4C72B0", edgecolor="white",
                 alpha=0.85)
        ax2.axvline(float(np.median(froh_total)), color="red", linestyle="--",
                    linewidth=1.2, label=f"Median={np.median(froh_total):.3f}")
        ax2.legend(fontsize=8)
    else:
        ax2.bar([ids[0]], froh_total, color="#4C72B0", alpha=0.85)
    ax2.set_xlabel("F_ROH (total)")
    ax2.set_ylabel("Count")
    ax2.set_title("FROH Distribution")

    plt.tight_layout()
    png_path = output_dir / "froh_distribution.png"
    svg_path = output_dir / "froh_distribution.svg"
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    return png_path, svg_path


# ---------------------------------------------------------------------------
# Internal helpers used by plot_roh()
# ---------------------------------------------------------------------------

def _plot_roh_histogram(ax: plt.Axes, segments: list[ROHSegment]) -> None:
    if not segments:
        ax.text(0.5, 0.5, "No ROH detected", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("ROH Length Distribution")
        return
    lengths_mb = [s.length_bp / 1e6 for s in segments]
    min_l = max(0.01, min(lengths_mb))
    max_l = max(lengths_mb) + 0.1
    bins = np.logspace(np.log10(min_l), np.log10(max_l), 25)
    for cls in ("short", "medium", "long"):
        sub = [s.length_bp / 1e6 for s in segments if s.length_class == cls]
        if sub:
            ax.hist(sub, bins=bins, alpha=0.7,
                    label=_CLASS_LABELS[cls],
                    color=_CLASS_COLORS[cls])
    ax.set_xscale("log")
    ax.set_xlabel("ROH Length (Mb)")
    ax.set_ylabel("Count")
    ax.set_title("ROH Length Distribution")
    ax.legend(fontsize=7)
    ax.axvline(0.1, color="gray", linestyle=":", linewidth=0.8)
    ax.axvline(1.0, color="orange", linestyle=":", linewidth=0.8)


def _plot_froh_barplot(ax: plt.Axes, froh_results: list[FROHResult]) -> None:
    if not froh_results:
        ax.text(0.5, 0.5, "No FROH data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("FROH per Individual")
        return
    ids = [fr.individual_id for fr in froh_results]
    frohs = [fr.froh for fr in froh_results]
    x = np.arange(len(ids))
    bars = ax.bar(x, frohs, color="#4C72B0", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ids, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("F_ROH")
    ax.set_title("FROH per Individual")
    ax.set_ylim(0, min(1.0, max(frohs) * 1.3 + 0.01))
    for bar, v in zip(bars, frohs):
        ax.text(bar.get_x() + bar.get_width() / 2, v + 0.001, f"{v:.4f}", ha="center", fontsize=7)


def _plot_ne_trajectory(ax: plt.Axes, ne_estimates: list[Any]) -> None:
    if not ne_estimates:
        ax.text(0.5, 0.5, "No Ne data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("Ne Trajectory")
        return
    individuals = sorted({e.individual_id for e in ne_estimates})
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(individuals), 1)))
    for ind, color in zip(individuals, colors):
        ind_ests = [e for e in ne_estimates if e.individual_id == ind]
        ind_ests_sorted = sorted(ind_ests, key=lambda e: e.generations_ago)
        gens = [e.generations_ago for e in ind_ests_sorted if e.ne < 1_000_000]
        nes = [e.ne for e in ind_ests_sorted if e.ne < 1_000_000]
        if gens:
            ax.plot(gens, nes, marker="o", label=ind, color=color, linewidth=1.5)
    ax.set_xlabel("Generations Ago")
    ax.set_ylabel("Effective Population Size (Ne)")
    ax.set_title("Ne Trajectory (from ROH lengths)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(fontsize=7)
    ax.invert_xaxis()


def _chrom_sort_key(chrom: str) -> tuple[int, str]:
    """Sort chromosomes numerically where possible (chr1 < chr2 … < chrX < chrY)."""
    name = chrom.lower().replace("chr", "")
    try:
        return (int(name), "")
    except ValueError:
        return (9999, name)
