import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os


# Tier colours
_TIER_COLOURS = {"high": "#2ca02c", "moderate": "#ff7f0e", "low": "#d62728"}


def plot_composite_scores(locus_results, output_dir):
    """Bar chart of composite tolerance scores per locus.

    locus_results: list of dicts with composite_tolerance, tolerance_tier,
                   locus_name, chrom, start, end.
    Returns list of written file paths.
    """
    if not locus_results:
        return []

    # Sort by composite tolerance descending
    sorted_res = sorted(locus_results,
                        key=lambda r: float(r.get("composite_tolerance", 0)),
                        reverse=True)

    labels = [r.get("locus_name", f"{r['chrom']}:{r['start']}-{r['end']}")
              for r in sorted_res]
    values = [float(r.get("composite_tolerance", 0)) for r in sorted_res]
    tiers = [r.get("tolerance_tier", "moderate") for r in sorted_res]
    colours = [_TIER_COLOURS.get(t, "#7f7f7f") for t in tiers]

    fig, ax = plt.subplots(figsize=(max(6, len(labels) * 1.4), 5))
    bars = ax.bar(labels, values, color=colours, edgecolor="white", linewidth=0.8)
    ax.axhline(0.7, color="green", linestyle="--", linewidth=1, alpha=0.7,
               label="High tolerance (0.7)")
    ax.axhline(0.4, color="orange", linestyle="--", linewidth=1, alpha=0.7,
               label="Moderate threshold (0.4)")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Composite Tolerance Score (0–1)")
    ax.set_title("Large-Insert Tolerance: Per-Locus Composite Scores")
    ax.set_xlabel("Locus")

    # Legend patches for tiers
    legend_patches = [
        mpatches.Patch(color=_TIER_COLOURS["high"], label="High tolerance"),
        mpatches.Patch(color=_TIER_COLOURS["moderate"], label="Moderate tolerance"),
        mpatches.Patch(color=_TIER_COLOURS["low"], label="Low tolerance"),
    ]
    ax.legend(handles=legend_patches + [
        plt.Line2D([0], [0], color="green", linestyle="--", label="High threshold (0.7)"),
        plt.Line2D([0], [0], color="orange", linestyle="--", label="Moderate threshold (0.4)"),
    ], fontsize=8, loc="upper right")

    # Label bar values
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.01,
                f"{val:.2f}", ha="center", va="bottom", fontsize=8)

    plt.tight_layout()
    paths = []
    for ext in ("png", "svg"):
        out = os.path.join(output_dir, f"composite_tolerance_scores.{ext}")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        paths.append(out)
    plt.close(fig)
    return paths


def plot_score_components(locus_results, output_dir):
    """Grouped bar chart breaking down score components per locus."""
    if not locus_results:
        return []

    sorted_res = sorted(locus_results,
                        key=lambda r: float(r.get("composite_tolerance", 0)),
                        reverse=True)
    labels = [r.get("locus_name", f"{r['chrom']}:{r['start']}")
              for r in sorted_res]

    component_keys = [
        ("chromatin_score", "#4C72B0", "Chromatin"),
        ("sequence_complexity_score", "#55A868", "Seq complexity"),
        ("gene_density_score", "#C44E52", "Gene density"),
        ("size_penalty", "#DD8452", "Size penalty"),
    ]

    n = len(labels)
    x = list(range(n))
    bar_width = 0.18
    fig, ax = plt.subplots(figsize=(max(7, n * 1.6), 5))

    for i, (key, colour, label) in enumerate(component_keys):
        vals = [float(r.get(key, 0)) for r in sorted_res]
        offsets = [xi + (i - 1.5) * bar_width for xi in x]
        ax.bar(offsets, vals, width=bar_width, color=colour, label=label,
               edgecolor="white", linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15, ha="right")
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("Score Component (0–1)")
    ax.set_title("Tolerance Score Components per Locus")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    paths = []
    for ext in ("png", "svg"):
        out = os.path.join(output_dir, f"score_components.{ext}")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        paths.append(out)
    plt.close(fig)
    return paths


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
