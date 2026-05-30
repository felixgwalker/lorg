"""
Output writers for the CRISPR Base Editor Window Visualiser.

Generates:
  - Editability CSV       — per-position annotation table
  - Outcomes CSV          — predicted edit products with frequencies
  - Bystander TXT         — plain-text warning report
  - Per-target summary TSV — one-row-per-target summary
  - Window coordinates TSV — editing window coordinates per guide
  - Bystander TSV          — bystander edit predictions table
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd

from src.analyser import PositionInfo, PamInfo
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


# ---------------------------------------------------------------------------
# TSV summary writers
# ---------------------------------------------------------------------------

def write_target_summary_tsv(
    positions: List[PositionInfo],
    editor: BaseEditorProfile,
    pam_info: Optional[PamInfo],
    output_prefix: str,
    guide_rna: str = "",
    target_dna: str = "",
) -> Path:
    """
    Write a one-row-per-target TSV summary file.

    Columns:
        guide_rna, target_dna, editor_name, editor_class, edit_type,
        window_start, window_end, pam_seq, pam_pattern, pam_ok,
        n_editable_in_window, n_primary_targets, n_bystanders,
        primary_target_positions, bystander_positions,
        max_predicted_editing_pct

    Args:
        positions:      20-element list from analyse_sequence().
        editor:         BaseEditorProfile used for analysis.
        pam_info:       PamInfo returned by analyse_sequence(), or None.
        output_prefix:  File path prefix; '_target_summary.tsv' appended.
        guide_rna:      Original guide RNA string (for record keeping).
        target_dna:     Original target DNA string (for record keeping).

    Returns:
        Path to the written TSV file.
    """
    primary_positions   = [p.position for p in positions if p.is_primary_target]
    bystander_positions = [p.position for p in positions if p.is_bystander]
    editable_in_window  = [p for p in positions if p.is_editable]
    max_editing_pct = (
        max(p.absolute_efficiency * 100 for p in editable_in_window)
        if editable_in_window else 0.0
    )

    row = {
        "guide_rna":                guide_rna,
        "target_dna":               target_dna,
        "editor_name":              editor.name,
        "editor_class":             editor.editor_class,
        "edit_type":                f"{editor.target_base}>{editor.product_base}",
        "window_start":             editor.window_start,
        "window_end":               editor.window_end,
        "pam_seq":                  pam_info.pam_seq if pam_info else "N/A",
        "pam_pattern":              pam_info.pam_pattern if pam_info else "NGG",
        "pam_ok":                   pam_info.pam_ok if pam_info else False,
        "n_editable_in_window":     len(editable_in_window),
        "n_primary_targets":        len(primary_positions),
        "n_bystanders":             len(bystander_positions),
        "primary_target_positions": ";".join(map(str, primary_positions)),
        "bystander_positions":      ";".join(map(str, bystander_positions)),
        "max_predicted_editing_pct": round(max_editing_pct, 1),
    }

    path = Path(f"{output_prefix}_target_summary.tsv")
    pd.DataFrame([row]).to_csv(path, index=False, sep="\t")
    logger.info("Target summary TSV saved: %s", path)
    return path


def write_window_coords_tsv(
    positions: List[PositionInfo],
    editor: BaseEditorProfile,
    output_prefix: str,
    guide_rna: str = "",
) -> Path:
    """
    Write a per-guide editing window coordinates TSV.

    Columns:
        guide_rna, editor_name, window_start, window_end,
        window_length, editable_positions, editable_bases,
        n_editable

    Args:
        positions:      20-element list from analyse_sequence().
        editor:         BaseEditorProfile used for analysis.
        output_prefix:  File path prefix; '_window_coords.tsv' appended.
        guide_rna:      Original guide RNA string (for record keeping).

    Returns:
        Path to the written TSV file.
    """
    editable_pos   = [p.position for p in positions if p.is_editable]
    editable_bases = [p.protospacer_base for p in positions if p.is_editable]

    row = {
        "guide_rna":          guide_rna,
        "editor_name":        editor.name,
        "window_start":       editor.window_start,
        "window_end":         editor.window_end,
        "window_length":      editor.window_end - editor.window_start + 1,
        "editable_positions": ";".join(map(str, editable_pos)),
        "editable_bases":     ";".join(editable_bases),
        "n_editable":         len(editable_pos),
    }

    path = Path(f"{output_prefix}_window_coords.tsv")
    pd.DataFrame([row]).to_csv(path, index=False, sep="\t")
    logger.info("Window coordinates TSV saved: %s", path)
    return path


def write_bystander_tsv(
    positions: List[PositionInfo],
    editor: BaseEditorProfile,
    bystander_threshold: float,
    output_prefix: str,
) -> Path:
    """
    Write bystander edit predictions as a TSV file.

    Only bystander positions are included.

    Columns:
        position, protospacer_base, template_base, edit_type,
        relative_efficiency, predicted_editing_pct, bystander_risk_score,
        risk_level

    Args:
        positions:            20-element list from analyse_sequence().
        editor:               BaseEditorProfile for edit type annotation.
        bystander_threshold:  Absolute editing frequency above which a bystander
                              is flagged HIGH.
        output_prefix:        File path prefix; '_bystander_predictions.tsv' appended.

    Returns:
        Path to the written TSV file.
    """
    rows = []
    for p in sorted(positions, key=lambda x: x.position):
        if not p.is_bystander:
            continue
        pct = p.absolute_efficiency * 100
        risk = "HIGH" if p.bystander_risk_score >= bystander_threshold else "LOW"
        rows.append({
            "position":              p.position,
            "protospacer_base":      p.protospacer_base,
            "template_base":         p.template_base,
            "edit_type":             f"{editor.target_base}>{editor.product_base}",
            "relative_efficiency":   round(p.relative_efficiency, 4),
            "predicted_editing_pct": round(pct, 1),
            "bystander_risk_score":  round(p.bystander_risk_score, 4),
            "risk_level":            risk,
        })

    columns = [
        "position", "protospacer_base", "template_base", "edit_type",
        "relative_efficiency", "predicted_editing_pct", "bystander_risk_score",
        "risk_level",
    ]
    path = Path(f"{output_prefix}_bystander_predictions.tsv")
    if rows:
        pd.DataFrame(rows, columns=columns).to_csv(path, index=False, sep="\t")
    else:
        pd.DataFrame(columns=columns).to_csv(path, index=False, sep="\t")

    logger.info("Bystander predictions TSV saved: %s", path)
    return path
