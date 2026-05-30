"""
Duplex diagram renderer for the CRISPR Base Editor Window Visualiser.

Produces a colour-coded guide RNA / target DNA duplex diagram with an
overlaid activity-window bracket and a per-position efficiency bar chart.
Outputs PNG, SVG, and an HTML/CSS interactive version.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

from src.analyser import PositionInfo, PamInfo
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
    pam_info: Optional[PamInfo] = None,
) -> Tuple[Path, Path, Path]:
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

    Also writes an HTML/CSS interactive version of the diagram.

    Args:
        positions:      List of 20 PositionInfo objects from analyse_sequence().
        editor:         BaseEditorProfile used for analysis.
        pam_seq:        PAM sequence extracted from target_dna (e.g. 'NGG').
        output_prefix:  File path prefix; '_duplex.png', '_duplex.svg', and
                        '_duplex.html' appended.
        pam_info:       Optional PamInfo for PAM match annotation in HTML.

    Returns:
        Tuple of (png_path, svg_path, html_path).
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

    png_path  = Path(f"{output_prefix}_duplex.png")
    svg_path  = Path(f"{output_prefix}_duplex.svg")
    html_path = Path(f"{output_prefix}_duplex.html")

    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)

    # ── HTML/CSS diagram ─────────────────────────────────────────────────
    _write_html_diagram(positions, editor, pam_seq, pam_info, html_path)

    logger.info("Duplex diagram saved: %s, %s, %s", png_path, svg_path, html_path)
    return png_path, svg_path, html_path


# ---------------------------------------------------------------------------
# HTML/CSS diagram generator
# ---------------------------------------------------------------------------

def _base_role_css_class(p: PositionInfo) -> str:
    """Return a CSS class name for the given position's editing role."""
    if p.is_primary_target:
        return "primary"
    if p.is_bystander:
        return "bystander"
    if p.in_window:
        return "in-window"
    return "out-window"


def _write_html_diagram(
    positions: List[PositionInfo],
    editor: BaseEditorProfile,
    pam_seq: str,
    pam_info: Optional[PamInfo],
    html_path: Path,
) -> None:
    """
    Write a self-contained HTML/CSS diagram showing the guide–target duplex.

    The diagram uses a horizontal flex layout with one column per protospacer
    position.  Each column shows:
      - Guide RNA base (coloured box with tooltip)
      - Pairing connector
      - Target DNA base (light box)
      - Position number
    Followed by PAM columns.
    Below the duplex a bar chart of predicted editing efficiency is rendered
    via inline SVG.

    Annotations:
      - Filled circle = editable base (primary or bystander)
      - Open circle = editable base that is a bystander
      - Shaded box (purple) = PAM
    """
    pam_display = (pam_seq[:3] if pam_seq else "NGG").upper()
    pam_match_str = ""
    if pam_info is not None:
        status = "matches" if pam_info.pam_ok else "MISMATCH"
        pam_match_str = (
            f'<p class="pam-status pam-{"ok" if pam_info.pam_ok else "mismatch"}">'
            f'PAM: <strong>{pam_info.pam_seq}</strong> '
            f'(required: <strong>{pam_info.pam_pattern}</strong>) — {status}</p>'
        )

    # Build guide / target rows
    guide_cells = []
    target_cells = []
    pos_labels = []

    for p in positions:
        cls = _base_role_css_class(p)
        tooltip_lines = [
            f"Position {p.position}",
            f"Guide: {p.guide_base} / Target: {p.protospacer_base}",
            f"Template: {p.template_base}",
        ]
        if p.in_window:
            tooltip_lines.append(f"In activity window")
        if p.is_primary_target:
            tooltip_lines.append(f"PRIMARY target — predicted {p.absolute_efficiency*100:.1f}%")
        elif p.is_bystander:
            tooltip_lines.append(f"BYSTANDER — predicted {p.absolute_efficiency*100:.1f}%")
        elif p.in_window:
            tooltip_lines.append(f"In window but not editable base")
        tooltip = "&#10;".join(tooltip_lines)

        # Circle annotation: filled = primary, open = bystander
        circle_html = ""
        if p.is_primary_target:
            circle_html = '<span class="circle filled" title="Primary edit target"></span>'
        elif p.is_bystander:
            circle_html = '<span class="circle open" title="Bystander risk"></span>'

        guide_cells.append(
            f'<div class="cell guide-cell {cls}" title="{tooltip}">'
            f'{p.guide_base}{circle_html}</div>'
        )
        target_cells.append(
            f'<div class="cell target-cell" title="{tooltip}">{p.protospacer_base}</div>'
        )
        pos_labels.append(
            f'<div class="cell pos-label">{p.position}</div>'
        )

    # PAM cells
    for base in pam_display:
        guide_cells.append(
            '<div class="cell guide-cell pam-guide" title="PAM (guide does not cover)">—</div>'
        )
        target_cells.append(
            f'<div class="cell target-cell pam-target" title="PAM base">{base}</div>'
        )
        pos_labels.append('<div class="cell pos-label pam-pos">PAM</div>')

    guide_row_html  = "\n".join(guide_cells)
    target_row_html = "\n".join(target_cells)
    pos_row_html    = "\n".join(pos_labels)

    # Inline SVG bar chart
    bar_w = 28
    bar_gap = 4
    chart_h = 80
    max_pct = editor.max_absolute_efficiency * 100 * 1.20
    chart_total_w = (bar_w + bar_gap) * (len(positions) + len(pam_display))

    bars_svg = []
    for i, p in enumerate(positions):
        pct = p.absolute_efficiency * 100
        bar_h = int((pct / max_pct) * chart_h) if max_pct > 0 else 0
        x = i * (bar_w + bar_gap)
        y = chart_h - bar_h
        color_map = {
            "primary":    "#27ae60",
            "bystander":  "#e67e22",
            "in-window":  "#2980b9",
            "out-window": "#95a5a6",
        }
        fill = color_map.get(_base_role_css_class(p), "#95a5a6")
        bars_svg.append(
            f'<rect x="{x}" y="{y}" width="{bar_w}" height="{bar_h}" '
            f'fill="{fill}" opacity="0.85">'
            f'<title>Pos {p.position}: {pct:.1f}%</title></rect>'
        )
        if pct > 0:
            bars_svg.append(
                f'<text x="{x + bar_w//2}" y="{max(y - 2, 8)}" '
                f'text-anchor="middle" font-size="8" fill="#333">{pct:.0f}%</text>'
            )

    threshold_y = int(chart_h - (10.0 / max_pct) * chart_h) if max_pct > 0 else chart_h // 2
    bars_svg.append(
        f'<line x1="0" y1="{threshold_y}" x2="{chart_total_w}" y2="{threshold_y}" '
        f'stroke="#95a5a6" stroke-width="1" stroke-dasharray="4,3" opacity="0.7">'
        f'<title>10% threshold</title></line>'
    )

    svg_chart = (
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{chart_total_w}" height="{chart_h + 20}" '
        f'style="display:block;margin-top:6px;">'
        + "\n".join(bars_svg) +
        f'<text x="{chart_total_w//2}" y="{chart_h + 16}" text-anchor="middle" '
        f'font-size="10" fill="#555">Protospacer position (1=PAM-distal → 20=PAM-proximal)</text>'
        f'</svg>'
    )

    # Legend
    legend_items = [
        ("primary", "#27ae60", "Primary edit target"),
        ("bystander", "#e67e22", "Bystander edit risk"),
        ("in-window", "#2980b9", "In window (wrong base)"),
        ("out-window", "#95a5a6", "Outside window"),
        ("pam", "#8e44ad", "PAM"),
    ]
    legend_html = "".join(
        f'<span class="legend-item"><span class="legend-swatch" '
        f'style="background:{color};"></span>{label}</span>'
        for _, color, label in legend_items
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>CRISPR Base Editor Window Visualiser — {editor.name}</title>
<style>
  body {{
    font-family: 'Segoe UI', Arial, sans-serif;
    background: #f9f9fb;
    color: #2c3e50;
    margin: 0;
    padding: 20px;
  }}
  h1 {{ font-size: 1.2em; margin-bottom: 4px; }}
  .subtitle {{ font-size: 0.85em; color: #555; margin-bottom: 14px; }}
  .duplex-wrapper {{
    overflow-x: auto;
    padding: 10px 0 6px;
  }}
  .strand-row {{
    display: flex;
    align-items: center;
    margin-bottom: 2px;
    min-width: max-content;
  }}
  .strand-label {{
    font-size: 0.75em;
    font-style: italic;
    color: #555;
    min-width: 140px;
    text-align: right;
    padding-right: 8px;
    white-space: nowrap;
  }}
  .strand-end {{
    font-size: 0.8em;
    font-weight: bold;
    padding-left: 6px;
    color: #2c3e50;
  }}
  .cell {{
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 28px;
    height: 28px;
    margin: 1px;
    font-size: 0.85em;
    font-weight: bold;
    border-radius: 4px;
    cursor: default;
    position: relative;
    user-select: none;
    flex-shrink: 0;
  }}
  .guide-cell.primary    {{ background: #27ae60; color: white; }}
  .guide-cell.bystander  {{ background: #e67e22; color: white; }}
  .guide-cell.in-window  {{ background: #2980b9; color: white; }}
  .guide-cell.out-window {{ background: #95a5a6; color: white; }}
  .guide-cell.pam-guide  {{ background: #c0a0d0; color: #7d3fa0; font-style: italic; }}
  .target-cell           {{ background: #ecf0f1; border: 1px solid #bdc3c7; color: #2c3e50; }}
  .target-cell.pam-target {{ background: #8e44ad; color: white; border: none; }}
  .pos-label             {{ font-size: 0.65em; color: #999; font-weight: normal; height: 16px; }}
  .pos-label.pam-pos     {{ color: #9b59b6; font-weight: bold; font-size: 0.60em; }}
  .connector-row {{
    display: flex;
    align-items: center;
    min-width: max-content;
  }}
  .connector {{
    display: inline-block;
    width: 30px;
    text-align: center;
    color: #7f8c8d;
    font-size: 0.7em;
    line-height: 10px;
    flex-shrink: 0;
  }}
  .connector-pam {{
    display: inline-block;
    width: 30px;
    text-align: center;
    color: #7f8c8d;
    font-size: 0.7em;
    opacity: 0.5;
    flex-shrink: 0;
  }}
  .circle {{
    display: inline-block;
    width: 7px;
    height: 7px;
    border-radius: 50%;
    position: absolute;
    top: 2px;
    right: 2px;
  }}
  .circle.filled {{ background: white; }}
  .circle.open   {{ background: transparent; border: 1.5px solid white; }}
  .legend {{
    display: flex;
    flex-wrap: wrap;
    gap: 10px;
    margin: 10px 0;
    font-size: 0.78em;
  }}
  .legend-item {{
    display: flex;
    align-items: center;
    gap: 5px;
  }}
  .legend-swatch {{
    display: inline-block;
    width: 14px;
    height: 14px;
    border-radius: 3px;
  }}
  .chart-section {{ margin-top: 8px; overflow-x: auto; }}
  .chart-label {{ font-size: 0.75em; color: #555; margin-bottom: 2px; }}
  .pam-status {{ font-size: 0.82em; margin: 6px 0; padding: 4px 8px; border-radius: 4px; }}
  .pam-ok {{ background: #d4efdf; color: #1a7a3a; }}
  .pam-mismatch {{ background: #fde8e4; color: #a93226; }}
  .window-info {{
    font-size: 0.80em;
    color: #c0392b;
    font-weight: bold;
    margin: 4px 0 8px 148px;
  }}
</style>
</head>
<body>
<h1>CRISPR Base Editor Window Visualiser &mdash; {editor.name}</h1>
<p class="subtitle">{editor.description}</p>
{pam_match_str}
<div class="legend">{legend_html}</div>
<div class="window-info">
  Activity window: positions {editor.window_start}&ndash;{editor.window_end}
  &nbsp;|&nbsp; Edit type: {editor.target_base}&rarr;{editor.product_base}
</div>
<div class="duplex-wrapper">
  <div class="strand-row">
    <span class="strand-label">Guide RNA</span>
    <span style="font-size:0.8em;font-weight:bold;padding-right:4px;">5'</span>
    {guide_row_html}
    <span class="strand-end">3'</span>
  </div>
  <div class="connector-row">
    <span class="strand-label"></span>
    <span style="width:28px;display:inline-block;"></span>
    {"".join(f'<span class="connector">|</span>' for _ in positions)}
    {"".join(f'<span class="connector-pam">&#183;</span>' for _ in pam_display)}
  </div>
  <div class="strand-row">
    <span class="strand-label">Target DNA (non-template)</span>
    <span style="font-size:0.8em;font-weight:bold;padding-right:4px;">3'</span>
    {target_row_html}
    <span class="strand-end">5'</span>
  </div>
  <div class="strand-row">
    <span class="strand-label"></span>
    <span style="width:28px;display:inline-block;"></span>
    {pos_row_html}
  </div>
</div>
<div class="chart-section">
  <p class="chart-label">Predicted editing efficiency per position (dashed line = 10% threshold):</p>
  {svg_chart}
</div>
</body>
</html>
"""

    html_path.write_text(html, encoding="utf-8")
    logger.info("HTML diagram saved: %s", html_path)
