import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os


# Colours per spec: collinear=green, inversion=orange, translocation=purple
TYPE_COLORS = {
    "collinear": "#2ca02c",
    "inversion": "#ff7f0e",
    "translocation": "#9467bd",
}


def _build_offsets(genome):
    """Build cumulative chromosome offsets for genome-wide coordinate space."""
    offsets = {}
    cumulative = 0
    for chrom in sorted(genome.keys()):
        offsets[chrom] = cumulative
        cumulative += len(genome[chrom]) + 5000
    return offsets, cumulative


def plot_dotplot(blocks, genome1, genome2, output_dir):
    """Generate a dot-plot PNG and SVG.

    X-axis: positions in genome1; Y-axis: positions in genome2.
    Collinear: green; inversions: orange; translocations: purple.
    """
    offsets1, total1 = _build_offsets(genome1)
    offsets2, total2 = _build_offsets(genome2)

    fig, ax = plt.subplots(figsize=(10, 10))

    # Draw chromosome boundary lines and labels for genome1 (x-axis)
    for chrom1, off1 in offsets1.items():
        chrom_len = len(genome1[chrom1])
        ax.axvline(off1, color="gray", linewidth=0.5, alpha=0.5)
        ax.axvline(off1 + chrom_len, color="gray", linewidth=0.5, alpha=0.5)
        ax.text(off1 + chrom_len / 2, total2 * 1.01, chrom1,
                ha="center", va="bottom", fontsize=7, rotation=30)

    # Draw chromosome boundary lines and labels for genome2 (y-axis)
    for chrom2, off2 in offsets2.items():
        chrom_len = len(genome2[chrom2])
        ax.axhline(off2, color="gray", linewidth=0.5, alpha=0.5)
        ax.axhline(off2 + chrom_len, color="gray", linewidth=0.5, alpha=0.5)
        ax.text(total1 * 1.01, off2 + chrom_len / 2, chrom2,
                ha="left", va="center", fontsize=7)

    for block in blocks:
        color = TYPE_COLORS.get(block["type"], "#9467bd")
        off1 = offsets1.get(block["g1_chrom"], 0)
        off2 = offsets2.get(block["g2_chrom"], 0)

        seeds = block["seeds"]
        # Seeds are (pos1, pos2) tuples from chain_builder
        xs = [off1 + s[0] for s in seeds]
        ys = [off2 + s[1] for s in seeds]

        ax.scatter(xs, ys, s=1, color=color, alpha=0.6, rasterized=True)
        if len(xs) >= 2:
            ax.plot([xs[0], xs[-1]], [ys[0], ys[-1]],
                    color=color, linewidth=1.2, alpha=0.7)

    legend_handles = [
        mpatches.Patch(color=v, label=k.capitalize())
        for k, v in TYPE_COLORS.items()
    ]
    ax.legend(handles=legend_handles, loc="upper left", fontsize=9)
    ax.set_xlabel("Genome 1 position (bp)")
    ax.set_ylabel("Genome 2 position (bp)")
    ax.set_title("Synteny Dot Plot")
    ax.set_xlim(0, total1 * 1.08)
    ax.set_ylim(0, total2 * 1.08)
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    paths = []
    for ext in ("png", "svg"):
        out = os.path.join(output_dir, f"dotplot.{ext}")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        paths.append(out)
    plt.close(fig)
    return paths
