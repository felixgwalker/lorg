import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
import matplotlib.patches as patches
import numpy as np
import os


TYPE_COLORS = {
    "collinear": "#1f77b4",
    "inversion": "#d62728",
    "translocation": "#2ca02c",
}


def _build_offsets(genome):
    offsets = {}
    cumulative = 0
    for chrom in sorted(genome.keys()):
        offsets[chrom] = cumulative
        cumulative += len(genome[chrom]) + 3000
    return offsets, cumulative


def _bezier_ribbon(ax, x1_left, x1_right, y1, x2_left, x2_right, y2, color, alpha=0.35):
    verts = [
        (x1_left, y1),
        (x1_left, (y1 + y2) / 2),
        (x2_left, (y1 + y2) / 2),
        (x2_left, y2),
        (x2_right, y2),
        (x2_right, (y1 + y2) / 2),
        (x1_right, (y1 + y2) / 2),
        (x1_right, y1),
        (x1_left, y1),
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.LINETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.CLOSEPOLY,
    ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor=color, edgecolor="none", alpha=alpha)
    ax.add_patch(patch)


def plot_ribbon(blocks, genome1, genome2, output_dir):
    offsets1, total1 = _build_offsets(genome1)
    offsets2, total2 = _build_offsets(genome2)

    total = max(total1, total2)

    fig, ax = plt.subplots(figsize=(14, 6))

    bar_height = 0.06
    y_genome1 = 0.85
    y_genome2 = 0.15

    for chrom in sorted(genome1.keys()):
        off = offsets1[chrom] / total
        width = len(genome1[chrom]) / total
        rect = mpatches.FancyBboxPatch(
            (off, y_genome1 - bar_height / 2), width, bar_height,
            boxstyle="round,pad=0.002", facecolor="#aec7e8", edgecolor="#1f77b4", linewidth=1.5
        )
        ax.add_patch(rect)
        ax.text(off + width / 2, y_genome1 + bar_height / 2 + 0.01, chrom,
                ha="center", va="bottom", fontsize=8)

    for chrom in sorted(genome2.keys()):
        off = offsets2[chrom] / total
        width = len(genome2[chrom]) / total
        rect = mpatches.FancyBboxPatch(
            (off, y_genome2 - bar_height / 2), width, bar_height,
            boxstyle="round,pad=0.002", facecolor="#ffbb78", edgecolor="#ff7f0e", linewidth=1.5
        )
        ax.add_patch(rect)
        ax.text(off + width / 2, y_genome2 - bar_height / 2 - 0.02, chrom,
                ha="center", va="top", fontsize=8)

    ax.text(-0.01, y_genome1, "Genome 1", ha="right", va="center", fontsize=10, fontweight="bold")
    ax.text(-0.01, y_genome2, "Genome 2", ha="right", va="center", fontsize=10, fontweight="bold")

    for block in blocks:
        color = TYPE_COLORS.get(block["type"], "purple")
        off1 = offsets1.get(block["g1_chrom"], 0)
        off2 = offsets2.get(block["g2_chrom"], 0)

        x1_left = (off1 + block["g1_start"]) / total
        x1_right = (off1 + block["g1_end"]) / total
        x2_left = (off2 + block["g2_start"]) / total
        x2_right = (off2 + block["g2_end"]) / total

        if block["type"] == "inversion":
            x2_left, x2_right = x2_right, x2_left

        _bezier_ribbon(ax, x1_left, x1_right, y_genome1 - bar_height / 2,
                       x2_left, x2_right, y_genome2 + bar_height / 2, color)

    legend_handles = [
        mpatches.Patch(color=v, alpha=0.7, label=k.capitalize())
        for k, v in TYPE_COLORS.items()
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=9)
    ax.set_xlim(-0.05, 1.02)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Synteny Ribbon Diagram", fontsize=13, pad=12)

    plt.tight_layout()
    paths = []
    for ext in ("png", "svg"):
        out = os.path.join(output_dir, f"ribbon.{ext}")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        paths.append(out)
    plt.close(fig)
    return paths
