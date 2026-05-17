import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os


def plot_locus_context(all_results, loci, output_dir):
    locus_names = list({r["locus_name"] for r in all_results})
    locus_names.sort()

    paths = []
    for locus_name in locus_names:
        rows = [r for r in all_results if r["locus_name"] == locus_name]
        if not rows:
            continue
        rows.sort(key=lambda r: r["start"])
        positions = [r["start"] for r in rows]
        totals = [r["total_score"] for r in rows]
        gd = [r["gene_density_score"] for r in rows]
        reg = [r["regulatory_score"] for r in rows]
        rep = [r["repeat_score"] for r in rows]
        cplx = [r["complexity_score"] for r in rows]

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True)

        ax1.plot(positions, totals, color="navy", linewidth=2, label="Total score")
        ax1.fill_between(positions, totals, alpha=0.15, color="navy")
        ax1.set_ylabel("Tolerance Score (0-100)")
        ax1.set_ylim(0, 105)
        ax1.set_title(f"Insertion Tolerance Landscape: {locus_name}")
        ax1.axhline(70, color="green", linestyle="--", linewidth=1, alpha=0.7, label="High tolerance (70)")
        ax1.axhline(40, color="orange", linestyle="--", linewidth=1, alpha=0.7, label="Low tolerance (40)")
        ax1.legend(loc="upper right", fontsize=8)
        ax1.grid(True, alpha=0.3)

        ax2.stackplot(
            positions, gd, reg, rep, cplx,
            labels=["Gene density (0-30)", "Regulatory (0-30)", "Repeat (0-20)", "Complexity (0-20)"],
            colors=["#4C72B0", "#DD8452", "#55A868", "#C44E52"],
            alpha=0.8,
        )
        ax2.set_ylabel("Component Scores")
        ax2.set_xlabel("Genomic Position (bp)")
        ax2.set_ylim(0, 105)
        ax2.legend(loc="upper right", fontsize=8)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        safe_name = locus_name.replace(":", "_").replace("/", "_")
        for ext in ("png", "svg"):
            out = os.path.join(output_dir, f"locus_context_{safe_name}.{ext}")
            fig.savefig(out, dpi=150, bbox_inches="tight")
            paths.append(out)
        plt.close(fig)

    return paths
