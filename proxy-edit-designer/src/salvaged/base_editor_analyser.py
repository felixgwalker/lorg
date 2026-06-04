"""Map base editor activity window onto a guide-target duplex.

Salvaged from crispr-base-editor-window-visualiser (deleted stage1f).
Import fixed: was 'from src.config import ...' -> now relative.
"""

from __future__ import annotations
import logging
from dataclasses import dataclass
from typing import List, Optional, Tuple

from .base_editor_config import BaseEditorProfile, get_pam_requirement, pam_matches

logger = logging.getLogger(__name__)

_COMPLEMENT: dict[str, str] = {
    "A": "T", "T": "A", "G": "C", "C": "G", "N": "N",
    "a": "t", "t": "a", "g": "c", "c": "g", "n": "n",
}


@dataclass
class PamInfo:
    pam_seq: str
    pam_pattern: str
    pam_ok: bool


@dataclass
class PositionInfo:
    position: int
    guide_base: str
    protospacer_base: str
    template_base: str
    in_window: bool
    is_editable: bool
    is_primary_target: bool
    is_bystander: bool
    relative_efficiency: float
    absolute_efficiency: float
    bystander_risk_score: float


def analyse_sequence(
    guide_rna: str,
    target_dna: str,
    editor: BaseEditorProfile,
    target_position: Optional[int] = None,
) -> Tuple[List[PositionInfo], PamInfo]:
    """Map the base editor activity window onto the guide-target duplex."""
    guide = guide_rna.upper().replace("U", "T").strip()
    target = target_dna.upper().replace("U", "T").strip()
    if len(guide) != 20:
        raise ValueError(f"Guide RNA must be exactly 20 nt; got {len(guide)} nt.")
    if len(target) < 20:
        raise ValueError(f"Target DNA must be at least 20 nt; got {len(target)} nt.")
    protospacer = target[:20]
    pam_pattern = get_pam_requirement(editor.name)
    pam_len = len(pam_pattern)
    extracted_pam = target[20:20 + pam_len] if len(target) >= 20 + pam_len else target[20:]
    pam_ok = pam_matches(extracted_pam, pam_pattern) if extracted_pam else False
    pam_info = PamInfo(pam_seq=extracted_pam or "N/A", pam_pattern=pam_pattern, pam_ok=pam_ok)
    if target_position is not None and not (1 <= target_position <= 20):
        raise ValueError(f"target_position must be 1-20; got {target_position}.")
    editable_in_window = [
        pos for pos in range(editor.window_start, editor.window_end + 1)
        if protospacer[pos - 1] == editor.target_base
    ]
    if target_position is not None:
        primary_pos: Optional[int] = target_position
    elif editable_in_window:
        primary_pos = max(editable_in_window, key=lambda p: editor.efficiency_profile[p - 1])
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
        positions.append(PositionInfo(
            position=pos, guide_base=gb, protospacer_base=pb, template_base=tb,
            in_window=in_window, is_editable=is_editable,
            is_primary_target=is_primary, is_bystander=is_bystander,
            relative_efficiency=rel_eff, absolute_efficiency=abs_eff,
            bystander_risk_score=abs_eff if is_bystander else 0.0,
        ))
    return positions, pam_info
