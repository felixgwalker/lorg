import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os


def plot_edit_burden(variants, burden, chrom_lengths, output_dir, fmt="png"):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    class_counts = burden["class_counts"]
    classes = list(class_counts.keys())
    counts = [class_counts[c] for c in classes]
    weights = {"SNV": 1, "SMALL_INS": 3, "SMALL_DEL": 3,
               "LARGE_INS": 10, "LARGE_DEL": 10, "SV_INS": 50, "SV_DEL": 50}
    weighted_counts = [class_counts[c] * weights.get(c, 1) for c in classes]

    x = np.arange(len(classes))
    width = 0.35
    bars1 = ax1.bar(x - width/2, counts, width, label="Raw count", color="steelblue", alpha=0.8)
    bars2 = ax1.bar(x + width/2, weighted_counts, width, label="Weighted burden", color="coral", alpha=0.8)
    ax1.set_xticks(x)
    ax1.set_xticklabels(classes, rotation=30, ha="right")
    ax1.set_ylabel("Count / Burden")
    ax1.set_title("Edit Counts by Variant Class")
    ax1.legend()

    for bar in bars1:
        h = bar.get_height()
        if h > 0:
            ax1.text(bar.get_x() + bar.get_width()/2, h, str(int(h)),
                     ha="center", va="bottom", fontsize=7)

    chroms = list(chrom_lengths.keys())
    if not chroms:
        ax2.text(0.5, 0.5, "No chromosome data", ha="center", va="center",
                 transform=ax2.transAxes)
    else:
        for ci, chrom in enumerate(chroms):
            clen = chrom_lengths[chrom]
            chrom_vars = [v for v in variants if v.get("chrom") == chrom]
            if not chrom_vars:
                continue
            positions = [v["pos"] for v in chrom_vars]
            n_bins = max(10, clen // 500)
            counts_hist, bin_edges = np.histogram(positions, bins=n_bins, range=(0, clen))
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            density = counts_hist / max(clen / 1000, 1)
            ax2.plot(bin_centers, density + ci * density.max() * 1.5,
                     label=chrom, linewidth=1.5)

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
