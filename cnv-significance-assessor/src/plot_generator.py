"""
CNV ideogram and significance distribution visualisation.

Produces a two-panel figure:

  Left panel  — horizontal chromosome ideogram showing each CNV as a coloured
                rectangle on its chromosome track.  CNVs are colour-coded by
                significance tier: green (LIKELY_BENIGN), orange (VUS),
                red (LIKELY_PATHOGENIC).

  Right panel — stacked bar chart of CNV counts per chromosome split by tier,
                and a small legend for the colour scheme.

matplotlib.use("Agg") is set before any plt import so the figure renders
correctly in headless / server environments.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")   # must precede all other matplotlib imports

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

from src.classifier import ScoredCNV, TIER_BENIGN, TIER_PATHOGENIC, TIER_VUS
from src.config import OUTFILE_PLOT_PNG, OUTFILE_PLOT_SVG

logger = logging.getLogger(__name__)

# ── Colour scheme ─────────────────────────────────────────────────────────
_TIER_COLOURS: dict[str, str] = {
    TIER_BENIGN:     "#2ecc71",   # green
    TIER_VUS:        "#f39c12",   # amber
    TIER_PATHOGENIC: "#e74c3c",   # red
}
_TIER_ORDER = [TIER_PATHOGENIC, TIER_VUS, TIER_BENIGN]
_CHROM_TRACK_COLOUR = "#dcdcdc"
_CHROM_BORDER_COLOUR = "#aaaaaa"


def generate_ideogram_plot(
    scored: list[ScoredCNV],
    output_dir: Path,
    sample_name: str = "",
    dpi: int = 150,
) -> tuple[Path, Path]:
    """
    Generate and save the CNV ideogram figure (PNG + SVG).

    Args:
        scored:      List of ScoredCNV from classifier.classify_cnvs().
        output_dir:  Directory for output files.
        sample_name: Optional label for the figure title.
        dpi:         PNG output resolution.

    Returns:
        Tuple of (png_path, svg_path).
    """
    if not scored:
        logger.warning("No scored CNVs; skipping plot generation.")
        raise ValueError("No CNVs to plot.")

    # ── Collect chromosome extents ─────────────────────────────────────
    chrom_max: dict[str, int] = {}
    for s in scored:
        rec = s.annotated.record
        chrom_max[rec.chrom] = max(chrom_max.get(rec.chrom, 0), rec.end)

    chroms = _sort_chroms(list(chrom_max.keys()))
    n_chroms = len(chroms)
    chrom_idx = {c: i for i, c in enumerate(chroms)}

    # ── Layout ────────────────────────────────────────────────────────
    fig_height = max(4.0, n_chroms * 0.45 + 1.5)
    fig, (ax_ideogram, ax_bar) = plt.subplots(
        1, 2,
        figsize=(14, fig_height),
        gridspec_kw={"width_ratios": [3, 1]},
    )
    fig.patch.set_facecolor("white")

    title = f"CNV Significance Assessor — {sample_name}" if sample_name else "CNV Significance Assessor"
    fig.suptitle(title, fontsize=13, fontweight="bold", y=1.01)

    # ── Ideogram panel ────────────────────────────────────────────────
    _draw_ideogram(ax_ideogram, scored, chroms, chrom_idx, chrom_max)

    # ── Stacked bar panel ─────────────────────────────────────────────
    _draw_bar_panel(ax_bar, scored, chroms)

    plt.tight_layout()

    png_path = output_dir / OUTFILE_PLOT_PNG
    svg_path = output_dir / OUTFILE_PLOT_SVG
    fig.savefig(str(png_path), dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(str(svg_path), format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    logger.info("CNV ideogram saved: %s, %s", png_path, svg_path)
    return png_path, svg_path


# ── Ideogram drawing ──────────────────────────────────────────────────────

def _draw_ideogram(
    ax: plt.Axes,
    scored: list[ScoredCNV],
    chroms: list[str],
    chrom_idx: dict[str, int],
    chrom_max: dict[str, int],
) -> None:
    """Draw chromosome tracks and CNV rectangles."""
    n_chroms = len(chroms)
    track_height = 0.6
    track_pad    = 1.0      # y-spacing between tracks

    genome_max = max(chrom_max.values()) if chrom_max else 1

    for chrom in chroms:
        y = chrom_idx[chrom] * track_pad
        length = chrom_max[chrom]

        # Grey chromosome background bar
        ax.add_patch(mpatches.FancyBboxPatch(
            (0, y - track_height / 2), length, track_height,
            boxstyle="round,pad=0",
            facecolor=_CHROM_TRACK_COLOUR,
            edgecolor=_CHROM_BORDER_COLOUR,
            linewidth=0.5,
            zorder=1,
        ))

    # CNV rectangles (draw in tier order so pathogenic is on top)
    for tier in _TIER_ORDER:
        colour = _TIER_COLOURS[tier]
        for s in scored:
            if s.significance_tier != tier:
                continue
            rec = s.annotated.record
            if rec.chrom not in chrom_idx:
                continue
            y   = chrom_idx[rec.chrom] * track_pad
            ax.add_patch(mpatches.Rectangle(
                (rec.start, y - track_height / 2),
                max(rec.size, genome_max * 0.002),   # min visible width
                track_height,
                facecolor=colour,
                edgecolor="none",
                alpha=0.85,
                zorder=2,
            ))

    # Axes cosmetics
    ax.set_xlim(0, genome_max * 1.02)
    ax.set_ylim(-track_pad * 0.5, (n_chroms - 1) * track_pad + track_pad * 0.5)

    ax.set_yticks([chrom_idx[c] * track_pad for c in chroms])
    ax.set_yticklabels(chroms, fontsize=8)
    ax.invert_yaxis()

    ax.xaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, _: _bp_label(x))
    )
    ax.tick_params(axis="x", labelsize=8)
    ax.set_xlabel("Genomic position", fontsize=9)
    ax.set_title("CNV distribution by chromosome", fontsize=10, fontweight="semibold")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color="#e5e5e5", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)

    # Legend
    legend_patches = [
        mpatches.Patch(facecolor=_TIER_COLOURS[t], label=t.replace("_", " ").title())
        for t in _TIER_ORDER
    ]
    ax.legend(
        handles=legend_patches,
        loc="lower right",
        fontsize=8,
        framealpha=0.8,
        title="Significance",
        title_fontsize=8,
    )


# ── Bar panel ─────────────────────────────────────────────────────────────

def _draw_bar_panel(
    ax: plt.Axes,
    scored: list[ScoredCNV],
    chroms: list[str],
) -> None:
    """Draw stacked count bars per chromosome."""
    # Count CNVs per chrom × tier
    counts: dict[str, dict[str, int]] = {c: {t: 0 for t in _TIER_ORDER} for c in chroms}
    for s in scored:
        chrom = s.annotated.record.chrom
        if chrom in counts:
            counts[chrom][s.significance_tier] += 1

    y_positions = np.arange(len(chroms))
    bar_height  = 0.6

    left_offset = np.zeros(len(chroms))
    for tier in reversed(_TIER_ORDER):    # benign on left → pathogenic on right
        widths = np.array([counts[c][tier] for c in chroms], dtype=float)
        ax.barh(
            y_positions,
            widths,
            left=left_offset,
            height=bar_height,
            color=_TIER_COLOURS[tier],
            label=tier.replace("_", " ").title(),
            alpha=0.9,
        )
        left_offset += widths

    ax.set_yticks(y_positions)
    ax.set_yticklabels(chroms, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("CNV count", fontsize=9)
    ax.set_title("Counts by chromosome", fontsize=10, fontweight="semibold")
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    ax.tick_params(axis="x", labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ── Utility functions ─────────────────────────────────────────────────────

def _sort_chroms(chroms: list[str]) -> list[str]:
    """
    Sort chromosome names in natural order.

    Order: chr1–22 (numeric), then chrX, chrY, chrM, then anything else
    lexicographically.
    """
    def _key(c: str) -> tuple:
        c_norm = c.lower().lstrip("chr")
        try:
            return (0, int(c_norm), "")
        except ValueError:
            order = {"x": 23, "y": 24, "m": 25, "mt": 25}
            return (1, order.get(c_norm, 99), c_norm)

    return sorted(chroms, key=_key)


def _bp_label(bp: float) -> str:
    """Format a base-pair count as a human-readable string."""
    if bp >= 1_000_000:
        return f"{bp / 1_000_000:.0f} Mb"
    if bp >= 1_000:
        return f"{bp / 1_000:.0f} kb"
    return f"{int(bp)} bp"
