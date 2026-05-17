"""
Duplex diagram renderer for the CRISPR Base Editor Window Visualiser.

Produces a colour-coded guide RNA / target DNA duplex diagram with an
overlaid activity-window bracket and a per-position efficiency bar chart.
Outputs both PNG and SVG.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

from src.analyser import PositionInfo
from src.config import BaseEditorProfile

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Colour palette
# ---------------------------------------------------------------------------

_PALETTE = {
    "primary":    "#27ae60",   # green  — primary edit target
    "bystander":  "#e67e22",   # orange — bystander risk
    "in_window":  "#2980b9",   # blue   — in window, wrong base type
    "out_window": "#95a5a6",   # grey   — outside activity window
    "pam_guide":  "#c0a0d0",   # light purple — PAM on guide row
    "pam_target": "#8e44ad",   # purple       — PAM on target row
    "text_light": "white",
    "text_dark":  "#2c3e50",
    "border":     "#bdc3c7",
    "pair_line":  "#7f8c8d",
    "window_ann": "#e74c3c",
}


def _pos_color(p: PositionInfo) -> str:
    if p.is_primary_target:
        return _PALETTE["primary"]
    if p.is_bystander:
        return _PALETTE["bystander"]
    if p.in_window:
        return _PALETTE["in_window"]
    return _PALETTE["out_window"]


# ---------------------------------------------------------------------------
# Main rendering function
# ---------------------------------------------------------------------------

def generate_duplex_diagram(
    positions: List[PositionInfo],
    editor: BaseEditorProfile,
    pam_seq: str,
    output_prefix: str,
) -> Tuple[Path, Path]:
    """
    Render the guide–target duplex diagram and efficiency bar chart.

    Layout (top → bottom):
        - Activity-window bracket annotation
        - Guide RNA boxes (5'→3'), colour-coded by editing role
        - Pairing lines
        - Target DNA boxes (non-template strand, 3'→5' complement below)
        - Position numbers
        - PAM boxes on both rows
        - Efficiency bar chart

    Args:
        positions:      List of 20 PositionInfo objects from analyse_sequence().
        editor:         BaseEditorProfile used for analysis.
        pam_seq:        PAM sequence extracted from target_dna (e.g. 'NGG').
        output_prefix:  File path prefix; '_duplex.png' and '_duplex.svg' appended.

    Returns:
        Tuple of (png_path, svg_path).
    """
    n = 20
    gap = 1.10       # horizontal spacing between base boxes
    bw = 0.88        # box width
    bh = 0.76        # box height
    pad = 0.07       # FancyBboxPatch padding

    # Row y-coordinates
    y_guide  = 2.7
    y_target = 1.2
    y_pos    = y_target - 0.65   # position number row

    # Total figure width scaled to sequence length
    fig_w = max(20, n * gap + 6)
    fig = plt.figure(figsize=(fig_w, 7.2))

    ax_dup = fig.add_axes([0.05, 0.30, 0.90, 0.60])
    ax_eff = fig.add_axes([0.05, 0.06, 0.90, 0.18])

    # ── Protospacer boxes ────────────────────────────────────────────────
    for p in positions:
        x = (p.position - 1) * gap
        color = _pos_color(p)

        # Guide RNA box
        ax_dup.add_patch(FancyBboxPatch(
            (x - bw / 2, y_guide - bh / 2), bw, bh,
            boxstyle=f"round,pad={pad}",
            facecolor=color, edgecolor="white", linewidth=0.8, zorder=2,
        ))
        ax_dup.text(
            x, y_guide, p.guide_base,
            ha="center", va="center", fontsize=12,
            fontweight="bold", color=_PALETTE["text_light"], zorder=3,
        )

        # Pairing line
        ax_dup.plot(
            [x, x],
            [y_guide - bh / 2 - 0.06, y_target + bh / 2 + 0.06],
            color=_PALETTE["pair_line"], lw=0.9, zorder=1,
        )

        # Target DNA box (non-template strand)
        ax_dup.add_patch(FancyBboxPatch(
            (x - bw / 2, y_target - bh / 2), bw, bh,
            boxstyle=f"round,pad={pad}",
            facecolor="#ecf0f1", edgecolor=_PALETTE["border"], linewidth=0.5, zorder=2,
        ))
        ax_dup.text(
            x, y_target, p.protospacer_base,
            ha="center", va="center", fontsize=12,
            color=_PALETTE["text_dark"], zorder=3,
        )

        # Position label
        ax_dup.text(
            x, y_pos, str(p.position),
            ha="center", va="center", fontsize=7,
            color="#95a5a6",
        )

    # ── PAM boxes ────────────────────────────────────────────────────────
    pam_display = (pam_seq[:3] if pam_seq else "NGG").upper()
    pam_gap = gap * 0.88
    for j, base in enumerate(pam_display):
        xp = n * gap + j * pam_gap

        # Guide row — pale indicator (guide does not cover PAM)
        ax_dup.add_patch(FancyBboxPatch(
            (xp - bw / 2, y_guide - bh / 2), bw, bh,
            boxstyle=f"round,pad={pad}",
            facecolor=_PALETTE["pam_guide"], edgecolor="white", linewidth=0.5, zorder=2,
        ))
        ax_dup.text(
            xp, y_guide, "—",
            ha="center", va="center", fontsize=11,
            color="#9b59b6", zorder=3,
        )

        # Target row — PAM bases
        ax_dup.add_patch(FancyBboxPatch(
            (xp - bw / 2, y_target - bh / 2), bw, bh,
            boxstyle=f"round,pad={pad}",
            facecolor=_PALETTE["pam_target"], edgecolor="white", linewidth=0.5, zorder=2,
        ))
        ax_dup.text(
            xp, y_target, base,
            ha="center", va="center", fontsize=12,
            fontweight="bold", color="white", zorder=3,
        )

        ax_dup.plot(
            [xp, xp],
            [y_guide - bh / 2 - 0.06, y_target + bh / 2 + 0.06],
            color=_PALETTE["pair_line"], lw=0.9, linestyle=":", zorder=1, alpha=0.6,
        )

    pam_label_x = n * gap + (len(pam_display) - 1) * pam_gap / 2
    ax_dup.text(
        pam_label_x, y_pos, "PAM",
        ha="center", va="center", fontsize=7, color="#9b59b6", fontweight="bold",
    )

    # ── Strand labels ────────────────────────────────────────────────────
    x_left  = -1.0
    x_right = n * gap + len(pam_display) * pam_gap + 0.1

    for y, l5, l3, label in [
        (y_guide,  "5'", "3'", "Guide RNA"),
        (y_target, "3'", "5'", "Target DNA\n(non-template)"),
    ]:
        ax_dup.text(x_left, y, l5, ha="right", va="center",
                    fontsize=9, fontweight="bold", color=_PALETTE["text_dark"])
        ax_dup.text(x_right, y, l3, ha="left", va="center",
                    fontsize=9, fontweight="bold", color=_PALETTE["text_dark"])
        ax_dup.text(x_left - 0.55, y, label, ha="right", va="center",
                    fontsize=8, color="#555", style="italic")

    # ── Activity window bracket ──────────────────────────────────────────
    win_pos = [p for p in positions if p.in_window]
    if win_pos:
        wx0 = (win_pos[0].position - 1) * gap - bw / 2 - 0.12
        wx1 = (win_pos[-1].position - 1) * gap + bw / 2 + 0.12
        y_bracket = y_guide + bh / 2 + 0.28

        ax_dup.annotate(
            "",
            xy=(wx1, y_bracket), xytext=(wx0, y_bracket),
            arrowprops=dict(arrowstyle="<->", color=_PALETTE["window_ann"], lw=1.8),
        )
        ax_dup.text(
            (wx0 + wx1) / 2, y_bracket + 0.20,
            f"Activity window ({editor.name})  ·  pos {editor.window_start}–{editor.window_end}",
            ha="center", va="bottom", fontsize=8.5,
            color=_PALETTE["window_ann"], fontweight="bold",
        )

    # ── Legend ───────────────────────────────────────────────────────────
    legend_elements = [
        mpatches.Patch(facecolor=_PALETTE["primary"],    label="Primary edit target"),
        mpatches.Patch(facecolor=_PALETTE["bystander"],  label="Bystander edit risk"),
        mpatches.Patch(facecolor=_PALETTE["in_window"],  label="In window (not editable base)"),
        mpatches.Patch(facecolor=_PALETTE["out_window"], label="Outside activity window"),
        mpatches.Patch(facecolor=_PALETTE["pam_target"], label="PAM"),
    ]
    ax_dup.legend(
        handles=legend_elements,
        loc="upper right", fontsize=8,
        framealpha=0.93, ncol=3,
        bbox_to_anchor=(1.0, 1.02),
    )

    ax_dup.set_xlim(-2.2, x_right + 0.5)
    ax_dup.set_ylim(0.2, 4.6)
    ax_dup.axis("off")

    # ── Efficiency bar chart ─────────────────────────────────────────────
    xs      = [(p.position - 1) * gap for p in positions]
    pcts    = [p.absolute_efficiency * 100 for p in positions]
    bcols   = [_pos_color(p) for p in positions]
    ylim    = editor.max_absolute_efficiency * 100 * 1.20

    ax_eff.bar(xs, pcts, width=bw * 0.92, color=bcols, edgecolor="white", lw=0.4)
    ax_eff.set_xlim(ax_dup.get_xlim()[0], ax_dup.get_xlim()[1])
    ax_eff.set_ylim(0, ylim)
    ax_eff.set_xticks(xs)
    ax_eff.set_xticklabels([str(p.position) for p in positions], fontsize=7)
    ax_eff.set_ylabel("Predicted\nediting (%)", fontsize=8)
    ax_eff.set_xlabel(
        "Protospacer position  (1 = PAM-distal  →  20 = PAM-proximal)",
        fontsize=8,
    )
    ax_eff.spines["top"].set_visible(False)
    ax_eff.spines["right"].set_visible(False)
    ax_eff.tick_params(axis="y", labelsize=7)
    # Threshold guideline at 10 %
    ax_eff.axhline(10, color="#95a5a6", lw=0.7, linestyle="--", alpha=0.6,
                   label="10 % threshold")

    fig.suptitle(
        f"CRISPR Base Editor Window Visualisation  ·  {editor.name}\n"
        f"{editor.description}",
        fontsize=11, fontweight="bold", y=0.99,
    )

    png_path = Path(f"{output_prefix}_duplex.png")
    svg_path = Path(f"{output_prefix}_duplex.svg")

    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)

    logger.info("Duplex diagram saved: %s, %s", png_path, svg_path)
    return png_path, svg_path
