"""
Sequence analysis module for the CRISPR Base Editor Window Visualiser.

Maps a base editor's activity window onto a guide RNA–target duplex and
classifies each protospacer position as primary target, bystander edit risk,
in-window non-editable, or outside the editing window.
"""

from __future__ import annotations

import logging
import sys
from dataclasses import dataclass
from typing import List, Optional

from src.config import BaseEditorProfile

logger = logging.getLogger(__name__)

# Single-base DNA complement lookup
_COMPLEMENT: dict[str, str] = {
    "A": "T", "T": "A", "G": "C", "C": "G", "N": "N",
    "a": "t", "t": "a", "g": "c", "c": "g", "n": "n",
}


@dataclass
class PositionInfo:
    """Per-nucleotide annotation for a single protospacer position."""

    position: int               # 1-indexed, 1 = PAM-distal, 20 = PAM-proximal
    guide_base: str             # RNA base displayed (A/U/G/C)
    protospacer_base: str       # DNA base on the non-template strand
    template_base: str          # DNA base on the template strand (complement)
    in_window: bool             # within the editor's declared activity window
    is_editable: bool           # correct base type for this editor class
    is_primary_target: bool     # designated primary edit target
    is_bystander: bool          # editable base in window but not primary target
    relative_efficiency: float  # 0–1, relative positional preference
    absolute_efficiency: float  # predicted absolute editing frequency (0–1)
    bystander_risk_score: float # absolute_efficiency when is_bystander, else 0.0


def analyse_sequence(
    guide_rna: str,
    target_dna: str,
    editor: BaseEditorProfile,
    target_position: Optional[int] = None,
) -> List[PositionInfo]:
    """
    Map the base editor activity window onto the guide–target duplex.

    The protospacer is taken as the first 20 nt of *target_dna* (non-template
    strand, 5'→3').  Position 1 is PAM-distal; position 20 is PAM-proximal.

    When multiple editable bases exist within the window the one with the
    highest relative efficiency is designated the primary target unless
    *target_position* is supplied explicitly.

    Args:
        guide_rna:        20-nt protospacer sequence, 5'→3' (DNA or RNA alphabet).
        target_dna:       Non-template strand including PAM, 5'→3' (DNA alphabet).
                          Minimum length: 20 nt.
        editor:           BaseEditorProfile describing window and efficiency curve.
        target_position:  1-indexed protospacer position of the intended primary
                          edit.  Overrides automatic selection when provided.

    Returns:
        List of 20 PositionInfo objects (index 0 = position 1).

    Raises:
        ValueError: On invalid sequence lengths or out-of-range target_position.
    """
    guide = guide_rna.upper().replace("U", "T").strip()
    target = target_dna.upper().replace("U", "T").strip()

    if len(guide) != 20:
        raise ValueError(
            f"Guide RNA must be exactly 20 nt (protospacer only); got {len(guide)} nt."
        )
    if len(target) < 20:
        raise ValueError(
            f"Target DNA must be at least 20 nt; got {len(target)} nt."
        )

    protospacer = target[:20]

    # Warn on guide/protospacer mismatches (tolerated — use protospacer)
    mismatches = [
        i + 1
        for i, (g, p) in enumerate(zip(guide, protospacer))
        if g != p
    ]
    if mismatches:
        logger.warning(
            "Guide/protospacer mismatches at position(s) %s — "
            "protospacer sequence will be used for analysis.",
            ", ".join(map(str, mismatches)),
        )

    if target_position is not None and not (1 <= target_position <= 20):
        raise ValueError(
            f"--target-position must be 1–20; got {target_position}."
        )

    # Editable positions within the window
    editable_in_window = [
        pos
        for pos in range(editor.window_start, editor.window_end + 1)
        if protospacer[pos - 1] == editor.target_base
    ]

    if not editable_in_window:
        logger.info(
            "No %s bases found in the activity window (positions %d–%d).  "
            "The guide cannot be edited by %s at this target.",
            editor.target_base, editor.window_start, editor.window_end, editor.name,
        )

    # Resolve primary target
    if target_position is not None:
        primary_pos: Optional[int] = target_position
        if primary_pos not in editable_in_window:
            logger.warning(
                "Specified --target-position %d is not an editable %s in the "
                "activity window.  It will be shown as primary target regardless.",
                primary_pos, editor.target_base,
            )
    elif editable_in_window:
        primary_pos = max(
            editable_in_window,
            key=lambda p: editor.efficiency_profile[p - 1],
        )
    else:
        primary_pos = None

    positions: List[PositionInfo] = []
    for i in range(20):
        pos = i + 1
        pb = protospacer[i]
        gb = pb.replace("T", "U")
        tb = _COMPLEMENT.get(pb, "N")

        in_window = editor.window_start <= pos <= editor.window_end
        is_editable = in_window and (pb == editor.target_base)
        is_primary = is_editable and (pos == primary_pos)
        is_bystander = is_editable and not is_primary

        rel_eff = editor.efficiency_profile[i] if in_window else 0.0
        abs_eff = rel_eff * editor.max_absolute_efficiency
        bsr = abs_eff if is_bystander else 0.0

        positions.append(PositionInfo(
            position=pos,
            guide_base=gb,
            protospacer_base=pb,
            template_base=tb,
            in_window=in_window,
            is_editable=is_editable,
            is_primary_target=is_primary,
            is_bystander=is_bystander,
            relative_efficiency=rel_eff,
            absolute_efficiency=abs_eff,
            bystander_risk_score=bsr,
        ))

    return positions
