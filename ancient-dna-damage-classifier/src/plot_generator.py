"""
Two-panel matplotlib damage profile visualisation.

Left panel:  5' C→T substitution rates by position (blue bars = observed,
             red line = fitted geometric decay model).
Right panel: 3' G→A substitution rates by position (green bars = observed,
             orange line = fitted decay model).

The x-axis of the left panel uses 1-based positions (1 = most terminal).
The x-axis of the right panel uses negative positions (-1 = most terminal),
mirroring the mapDamage2 convention so the two panels face inward.

matplotlib.use("Agg") is called before any other matplotlib import to ensure
the non-interactive backend is used in headless environments (Windows servers,
CI pipelines).  This call must precede all other matplotlib imports.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # must precede all other matplotlib imports

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

from src.config import OUTFILE_PLOT_PNG, OUTFILE_PLOT_SVG
from src.damage_profiler import DamageProfile
from src.decay_model import ModelResult

logger = logging.getLogger(__name__)

# ── Colour palette ────────────────────────────────────────────────────────────
_COLOUR_CT_BAR  = "#4878CF"   # blue — 5' C→T observed
_COLOUR_CT_LINE = "#D65F5F"   # red  — 5' fitted model
_COLOUR_GA_BAR  = "#6ACC65"   # green — 3' G→A observed
_COLOUR_GA_LINE = "#E87D2B"   # orange — 3' fitted model
_COLOUR_GRID    = "#E5E5E5"


def generate_damage_plot(
    profile: DamageProfile,
    model: ModelResult,
    output_dir: Path,
    sample_name: str = "",
    dpi: int = 150,
) -> tuple[Path, Path]:
    """
    Generate and save the two-panel damage profile figure.

    Args:
        profile:     DamageProfile from damage_profiler.profile_damage().
        model:       ModelResult from decay_model.fit_model().
        output_dir:  Directory to write output files into.
        sample_name: Optional label for the figure title.
        dpi:         DPI for PNG output.

    Returns:
        Tuple of (png_path, svg_path) as Path objects.
    """
    n = profile.n_terminal

    # x-axis: 1-based for 5' (1..n), negative for 3' (-n..-1 with -1 most terminal)
    five_prime_positions  = np.arange(1, n + 1)
    three_prime_positions = np.arange(-n, 0)         # -n, ..., -1

    ct_rates = profile.ct_rate
    ga_rates = profile.ga_rate
    ct_fitted = model.five_prime.fitted_values
    ga_fitted = model.three_prime.fitted_values

    # Shared y-axis maximum with a 10% margin
    y_max = max(ct_rates.max(), ga_rates.max(), ct_fitted.max(), ga_fitted.max())
    y_max = y_max * 1.10 if y_max > 0 else 0.05

    fig, (ax_five, ax_three) = plt.subplots(
        1, 2, figsize=(12, 4.5), sharey=False
    )
    fig.patch.set_facecolor("white")

    title = f"Ancient DNA Damage Profile — {sample_name}" if sample_name else "Ancient DNA Damage Profile"
    fig.suptitle(title, fontsize=13, fontweight="bold", y=1.01)

    _plot_terminus_panel(
        ax=ax_five,
        observed_rates=ct_rates,
        fitted_values=ct_fitted,
        positions=five_prime_positions,
        terminus_label="5′ C→T",
        bar_color=_COLOUR_CT_BAR,
        line_color=_COLOUR_CT_LINE,
        signal_quality=model.five_prime.signal_quality,
        y_max=y_max,
        x_label="Position from 5′ end",
        r_squared=model.five_prime.r_squared,
    )

    _plot_terminus_panel(
        ax=ax_three,
        observed_rates=ga_rates,
        fitted_values=ga_fitted,
        positions=three_prime_positions,
        terminus_label="3′ G→A",
        bar_color=_COLOUR_GA_BAR,
        line_color=_COLOUR_GA_LINE,
        signal_quality=model.three_prime.signal_quality,
        y_max=y_max,
        x_label="Position from 3′ end",
        r_squared=model.three_prime.r_squared,
    )

    # Annotation footer
    fig.text(
        0.5, -0.04,
        f"Reads profiled: {profile.n_reads_passed:,}  |  "
        f"Library deamination rate: {model.library_deamination_rate:.4f}  |  "
        f"Overall signal: {model.overall_signal_quality}",
        ha="center", va="bottom", fontsize=9, color="#555555",
    )

    plt.tight_layout()

    png_path = output_dir / OUTFILE_PLOT_PNG
    svg_path = output_dir / OUTFILE_PLOT_SVG

    fig.savefig(str(png_path), dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(str(svg_path), format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    logger.info("Damage profile plot saved: %s, %s", png_path, svg_path)
    return png_path, svg_path


def _plot_terminus_panel(
    ax: plt.Axes,
    observed_rates: np.ndarray,
    fitted_values: np.ndarray,
    positions: np.ndarray,
    terminus_label: str,
    bar_color: str,
    line_color: str,
    signal_quality: str,
    y_max: float,
    x_label: str,
    r_squared: float,
) -> None:
    """
    Draw one terminus panel with observed bars and fitted model line.

    Args:
        ax:              Matplotlib Axes to draw on.
        observed_rates:  Observed substitution rates per position.
        fitted_values:   Fitted model values per position.
        positions:       x-axis values (1-based for 5', negative for 3').
        terminus_label:  Panel title.
        bar_color:       Color for observation bars.
        line_color:      Color for fitted model line.
        signal_quality:  Text annotation for signal grade.
        y_max:           Upper y-axis limit.
        x_label:         x-axis label text.
        r_squared:       R² for fitted model (shown in annotation).
    """
    bar_width = 0.7

    ax.bar(
        positions,
        observed_rates,
        width=bar_width,
        color=bar_color,
        alpha=0.75,
        label="Observed",
        zorder=2,
    )
    ax.plot(
        positions,
        fitted_values,
        color=line_color,
        linewidth=2.0,
        marker="o",
        markersize=3.5,
        label=f"Fitted model (R²={r_squared:.3f})",
        zorder=3,
    )

    ax.set_ylim(0.0, y_max)
    ax.set_xlabel(x_label, fontsize=10)
    ax.set_ylabel("Substitution frequency", fontsize=10)
    ax.set_title(terminus_label, fontsize=11, fontweight="semibold")

    ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1.0, decimals=1))
    ax.set_xticks(positions)
    ax.tick_params(axis="x", labelsize=8, rotation=45 if len(positions) > 15 else 0)
    ax.tick_params(axis="y", labelsize=9)

    ax.grid(axis="y", color=_COLOUR_GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Signal quality badge in upper-right corner
    _quality_color = {
        "strong": "#2ecc71", "moderate": "#f39c12",
        "weak": "#e74c3c", "absent": "#95a5a6",
    }.get(signal_quality, "#95a5a6")

    ax.text(
        0.97, 0.95,
        f"Signal: {signal_quality}",
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=8.5,
        color="white",
        bbox=dict(
            boxstyle="round,pad=0.3",
            facecolor=_quality_color,
            edgecolor="none",
            alpha=0.9,
        ),
    )

    ax.legend(fontsize=8.5, framealpha=0.7, loc="upper right" if positions[0] > 0 else "lower left")
