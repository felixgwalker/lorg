"""Heatmap + dendrogram plot for phylogenetic distances."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform


def plot_heatmap_dendrogram(
    names: list[str],
    matrix: np.ndarray,
    output_dir: Path,
) -> tuple[Path, Path]:
    """Render a heatmap with dendrogram and save PNG + SVG."""
    n = len(names)
    if n < 2:
        raise ValueError("Need at least 2 sequences to plot.")

    condensed = squareform(matrix, checks=False)
    Z = linkage(condensed, method="average")

    fig = plt.figure(figsize=(max(10, n * 1.2), max(8, n * 1.0)))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 4], height_ratios=[4, 1],
                           wspace=0.02, hspace=0.02)

    ax_dendro_left = fig.add_subplot(gs[0, 0])
    ax_heatmap = fig.add_subplot(gs[0, 1])
    ax_dendro_top = fig.add_subplot(gs[1, 1])

    dend_left = dendrogram(
        Z,
        ax=ax_dendro_left,
        orientation="left",
        labels=names,
        leaf_font_size=9,
        color_threshold=0.7 * max(Z[:, 2]) if len(Z) > 0 else 1.0,
    )
    ax_dendro_left.set_xlabel("Distance")
    ax_dendro_left.set_title("Dendrogram")
    ax_dendro_left.invert_yaxis()

    order = dend_left["leaves"]
    ordered_names = [names[i] for i in order]
    ordered_matrix = matrix[np.ix_(order, order)]

    im = ax_heatmap.imshow(ordered_matrix, aspect="auto", cmap="YlOrRd", interpolation="nearest")
    ax_heatmap.set_xticks(range(n))
    ax_heatmap.set_xticklabels(ordered_names, rotation=45, ha="left", fontsize=8)
    ax_heatmap.xaxis.set_label_position("top")
    ax_heatmap.xaxis.tick_top()
    ax_heatmap.set_yticks(range(n))
    ax_heatmap.set_yticklabels(ordered_names, fontsize=8)
    ax_heatmap.set_title("Pairwise Distance Heatmap", pad=20)

    for i in range(n):
        for j in range(n):
            val = ordered_matrix[i, j]
            text_color = "white" if val > ordered_matrix.max() * 0.6 else "black"
            ax_heatmap.text(j, i, f"{val:.3f}", ha="center", va="center",
                            fontsize=6, color=text_color)

    plt.colorbar(im, ax=ax_heatmap, fraction=0.046, pad=0.04, label="Distance")

    dend_top = dendrogram(
        Z,
        ax=ax_dendro_top,
        orientation="bottom",
        labels=names,
        leaf_font_size=8,
        no_labels=True,
        color_threshold=0.7 * max(Z[:, 2]) if len(Z) > 0 else 1.0,
    )
    ax_dendro_top.set_visible(False)

    ax_dendro_left.set_frame_on(False)

    png_path = output_dir / "phylogenetic_heatmap.png"
    svg_path = output_dir / "phylogenetic_heatmap.svg"
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    return png_path, svg_path
