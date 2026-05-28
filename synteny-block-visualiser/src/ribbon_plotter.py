import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
import matplotlib.patches as patches
import numpy as np
import os


# Colours per spec: collinear=green, inversion=orange, translocation=purple
TYPE_COLORS = {
    "collinear": "#2ca02c",
    "inversion": "#ff7f0e",
    "translocation": "#9467bd",
}


def _build_offsets(genome):
    """Build cumulative chromosome offsets (normalised 0-1 later)."""
    offsets = {}
    cumulative = 0
    for chrom in sorted(genome.keys()):
        offsets[chrom] = cumulative
        cumulative += len(genome[chrom]) + 3000
    return offsets, cumulative


def _bezier_ribbon(ax, x1_left, x1_right, y1, x2_left, x2_right, y2, color, alpha=0.35):
    """Draw a filled Bezier ribbon connecting two horizontal intervals.

    The ribbon goes from [x1_left, x1_right] at height y1 to
    [x2_left, x2_right] at height y2, with cubic Bezier curves on the sides.

    A valid cubic-Bezier closed path requires the vertex list to satisfy
    matplotlib's Path conventions:
      MOVETO  + 3*CURVE4  (left side, top->bottom)
      LINETO              (bottom right)
      3*CURVE4            (right side, bottom->top)
      CLOSEPOLY           (close back to start)
    """
    mid_y = (y1 + y2) / 2

    verts = [
        # Left edge: cubic Bezier from (x1_left, y1) to (x2_left, y2)
        (x1_left, y1),          # MOVETO
        (x1_left, mid_y),       # CURVE4 cp1
        (x2_left, mid_y),       # CURVE4 cp2
        (x2_left, y2),          # CURVE4 endpoint
        # Bottom edge: straight line to right side of genome2
        (x2_right, y2),         # LINETO
        # Right edge: cubic Bezier from (x2_right, y2) back to (x1_right, y1)
        (x2_right, mid_y),      # CURVE4 cp1
        (x1_right, mid_y),      # CURVE4 cp2
        (x1_right, y1),         # CURVE4 endpoint
        # Close back to start
        (x1_left, y1),          # CLOSEPOLY
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,   # left Bezier
        Path.LINETO,                              # bottom right
        Path.CURVE4, Path.CURVE4, Path.CURVE4,   # right Bezier
        Path.CLOSEPOLY,
    ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor=color, edgecolor="none", alpha=alpha)
    ax.add_patch(patch)


def plot_ribbon(blocks, genome1, genome2, output_dir):
    """Generate a Circos-style ribbon diagram PNG and SVG.

    Top track: genome1 chromosomes; bottom track: genome2 chromosomes.
    Ribbons connect syntenic blocks, coloured by type.
    """
    offsets1, total1 = _build_offsets(genome1)
    offsets2, total2 = _build_offsets(genome2)

    total = max(total1, total2)

    fig, ax = plt.subplots(figsize=(14, 6))

    bar_height = 0.06
    y_genome1 = 0.85
    y_genome2 = 0.15

    # Draw genome1 chromosome track (top)
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

    # Draw genome2 chromosome track (bottom)
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

    # Draw ribbons
    for block in blocks:
        color = TYPE_COLORS.get(block["type"], "#9467bd")
        off1 = offsets1.get(block["g1_chrom"], 0)
        off2 = offsets2.get(block["g2_chrom"], 0)

        x1_left = (off1 + block["g1_start"]) / total
        x1_right = (off1 + block["g1_end"]) / total
        x2_left = (off2 + block["g2_start"]) / total
        x2_right = (off2 + block["g2_end"]) / total

        # Inversion: swap genome2 endpoints to show crossed ribbon
        if block["type"] == "inversion":
            x2_left, x2_right = x2_right, x2_left

        # Only draw ribbons with non-zero width
        if abs(x1_right - x1_left) < 1e-6 or abs(x2_right - x2_left) < 1e-6:
            continue

        _bezier_ribbon(ax, x1_left, x1_right, y_genome1 - bar_height / 2,
                       x2_left, x2_right, y_genome2 + bar_height / 2, color)

    legend_handles = [
        mpatches.Patch(color=v, alpha=0.7, label=k.capitalize())
        for k, v in TYPE_COLORS.items()
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=9)
    ax.set_xlim(-0.08, 1.02)
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
