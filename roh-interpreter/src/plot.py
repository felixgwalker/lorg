"""ROH distribution histogram, FROH barplot, and Ne trajectory plot."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from src.froh_calculator import FROHResult
from src.ne_estimator import NeEstimate
from src.roh_detector import ROHSegment


def plot_roh(
    segments: list[ROHSegment],
    froh_results: list[FROHResult],
    ne_estimates: list[NeEstimate],
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


def _plot_roh_histogram(ax: plt.Axes, segments: list[ROHSegment]) -> None:
    if not segments:
        ax.text(0.5, 0.5, "No ROH detected", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("ROH Length Distribution")
        return
    lengths_mb = [s.length_bp / 1e6 for s in segments]
    bins = np.logspace(np.log10(max(0.01, min(lengths_mb))), np.log10(max(lengths_mb) + 0.1), 25)
    colors_map = {"SHORT": "#4C72B0", "MEDIUM": "#DD8452", "LONG": "#C44E52", "VERY_SHORT": "#8EBA42"}
    for cls in ["VERY_SHORT", "SHORT", "MEDIUM", "LONG"]:
        sub = [s.length_bp / 1e6 for s in segments if s.length_class == cls]
        if sub:
            ax.hist(sub, bins=bins, alpha=0.7, label=cls, color=colors_map[cls])
    ax.set_xscale("log")
    ax.set_xlabel("ROH Length (Mb)")
    ax.set_ylabel("Count")
    ax.set_title("ROH Length Distribution")
    ax.legend(fontsize=7)
    ax.axvline(0.1, color="gray", linestyle=":", linewidth=0.8)
    ax.axvline(1.0, color="orange", linestyle=":", linewidth=0.8)
    ax.axvline(10.0, color="red", linestyle=":", linewidth=0.8)


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


def _plot_ne_trajectory(ax: plt.Axes, ne_estimates: list[NeEstimate]) -> None:
    if not ne_estimates:
        ax.text(0.5, 0.5, "No Ne data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("Ne Trajectory")
        return
    individuals = sorted({e.individual_id for e in ne_estimates})
    colors = plt.cm.tab10(np.linspace(0, 1, len(individuals)))
    gen_order = [50.0, 200.0, 1000.0]
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
