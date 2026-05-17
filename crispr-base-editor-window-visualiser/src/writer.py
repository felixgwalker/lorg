"""
Output writers for the CRISPR Base Editor Window Visualiser.

Generates:
  - Editability CSV   — per-position annotation table
  - Outcomes CSV      — predicted edit products with frequencies
  - Bystander TXT     — plain-text warning report
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Tuple

import pandas as pd

from src.analyser import PositionInfo
from src.config import BaseEditorProfile

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Editability table
# ---------------------------------------------------------------------------

def write_editability_csv(
    positions: List[PositionInfo],
    output_prefix: str,
) -> Path:
    """
    Write the per-position editability table to CSV.

    Columns:
        position, guide_base, protospacer_base, template_base,
        in_window, editable, is_primary_target, is_bystander,
        relative_efficiency, predicted_editing_pct, bystander_risk_score

    Args:
        positions:      20-element list from analyse_sequence().
        output_prefix:  File path prefix; '_editability.csv' appended.

    Returns:
        Path to the written CSV file.
    """
    rows = [
        {
            "position":              p.position,
            "guide_base":            p.guide_base,
            "protospacer_base":      p.protospacer_base,
            "template_base":         p.template_base,
            "in_window":             p.in_window,
            "editable":              p.is_editable,
            "is_primary_target":     p.is_primary_target,
            "is_bystander":          p.is_bystander,
            "relative_efficiency":   round(p.relative_efficiency, 4),
            "predicted_editing_pct": round(p.absolute_efficiency * 100, 1),
            "bystander_risk_score":  round(p.bystander_risk_score, 4),
        }
        for p in positions
    ]

    path = Path(f"{output_prefix}_editability.csv")
    pd.DataFrame(rows).to_csv(path, index=False)
    logger.info("Editability table saved: %s", path)
    return path


# ---------------------------------------------------------------------------
# Outcome prediction table
# ---------------------------------------------------------------------------

def write_outcomes_csv(
    positions: List[PositionInfo],
    editor: BaseEditorProfile,
    output_prefix: str,
) -> Path:
    """
    Write the edit outcome prediction table to CSV.

    Only editable positions (primary target and bystanders) are included.

    Columns:
        position, original_base, edited_base, edit_type, role,
        relative_efficiency, predicted_editing_pct

    Args:
        positions:      20-element list from analyse_sequence().
        editor:         BaseEditorProfile for edit type annotation.
        output_prefix:  File path prefix; '_outcomes.csv' appended.

    Returns:
        Path to the written CSV file.
    """
    rows = [
        {
            "position":              p.position,
            "original_base":         p.protospacer_base,
            "edited_base":           editor.product_base,
            "edit_type":             f"{editor.target_base}→{editor.product_base}",
            "role":                  "primary_target" if p.is_primary_target else "bystander",
            "relative_efficiency":   round(p.relative_efficiency, 4),
            "predicted_editing_pct": round(p.absolute_efficiency * 100, 1),
        }
        for p in positions
        if p.is_editable
    ]

    path = Path(f"{output_prefix}_outcomes.csv")
    if rows:
        pd.DataFrame(rows).to_csv(path, index=False)
    else:
        # Write an empty file with headers so downstream tools don't error
        pd.DataFrame(columns=[
            "position", "original_base", "edited_base", "edit_type",
            "role", "relative_efficiency", "predicted_editing_pct",
        ]).to_csv(path, index=False)

    logger.info("Outcomes table saved: %s", path)
    return path


# ---------------------------------------------------------------------------
# Bystander warning report
# ---------------------------------------------------------------------------

def write_bystander_warnings(
    positions: List[PositionInfo],
    editor: BaseEditorProfile,
    bystander_threshold: float,
    output_prefix: str,
) -> Tuple[Path, List[str]]:
    """
    Write the bystander edit warning report to a plain-text file.

    Each bystander position is classified as HIGH (predicted editing ≥ threshold)
    or LOW.

    Args:
        positions:            20-element list from analyse_sequence().
        editor:               BaseEditorProfile for edit type annotation.
        bystander_threshold:  Absolute editing frequency above which a bystander
                              is flagged HIGH (e.g. 0.10 for 10 %).
        output_prefix:        File path prefix; '_bystander_warnings.txt' appended.

    Returns:
        Tuple of (file_path, list_of_warning_lines).
    """
    bystanders = sorted(
        [p for p in positions if p.is_bystander],
        key=lambda x: x.bystander_risk_score,
        reverse=True,
    )

    lines: List[str] = []
    lines.append("=" * 60)
    lines.append(f"Bystander Edit Warning Report  —  {editor.name}")
    lines.append(f"Edit type : {editor.target_base}→{editor.product_base}")
    lines.append(f"Window    : positions {editor.window_start}–{editor.window_end}")
    lines.append(f"Threshold : {bystander_threshold * 100:.0f}% predicted editing → HIGH")
    lines.append("=" * 60)

    if not bystanders:
        lines.append("")
        lines.append("No bystander edit risks detected within the activity window.")
    else:
        lines.append("")
        lines.append(
            f"{'Pos':>4}  {'Base':>4}  {'Edit':>7}  "
            f"{'Rel.Eff':>8}  {'Pred.%':>7}  {'Risk':>5}"
        )
        lines.append("-" * 46)
        for p in bystanders:
            pct = p.absolute_efficiency * 100
            risk = "HIGH" if p.bystander_risk_score >= bystander_threshold else "LOW"
            lines.append(
                f"{p.position:>4}  {p.protospacer_base:>4}  "
                f"{editor.target_base}→{editor.product_base}  "
                f"{p.relative_efficiency:>8.4f}  {pct:>6.1f}%  {risk:>5}"
            )

    lines.append("")

    path = Path(f"{output_prefix}_bystander_warnings.txt")
    path.write_text("\n".join(lines), encoding="utf-8")
    logger.info("Bystander warnings saved: %s", path)
    return path, lines
