"""Introduce silent PAM-breaking mutations in HDR template right arm.

Salvaged from hdr-template-designer (deleted stage1f).
"""

from typing import TypedDict

_PAM_SCAN_RADIUS = 10
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


class SuggestedMutation(TypedDict):
    pos: int
    ref: str
    alt: str
    rationale: str


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def suggest_pam_disruption(template: str, pam_seq: str, cut_pos: int) -> list[SuggestedMutation]:
    """Scan for NGG PAM sites within ±10 bp of cut_pos and propose disrupting mutations."""
    window_start = max(0, cut_pos - _PAM_SCAN_RADIUS)
    window_end = min(len(template), cut_pos + _PAM_SCAN_RADIUS + 3)
    suggestions: list[SuggestedMutation] = []
    seen_positions: set[int] = set()
    for i in range(window_start, window_end - 2):
        if template[i + 1] == "G" and template[i + 2] == "G":
            _add_suggestion(suggestions, seen_positions, template, i, "+")
    for i in range(window_start, window_end - 2):
        if template[i] == "C" and template[i + 1] == "C":
            _add_suggestion_rev(suggestions, seen_positions, template, i, "-")
    suggestions.sort(key=lambda s: abs(s["pos"] - cut_pos))
    return suggestions


def _add_suggestion(suggestions, seen, template, pam_start, strand):
    target = pam_start + 2
    if target in seen:
        return
    seen.add(target)
    ref_base = template[target]
    alt_base = "A" if ref_base == "G" else "G"
    suggestions.append(SuggestedMutation(
        pos=target, ref=ref_base, alt=alt_base,
        rationale=f"NGG PAM ({strand}) at {pam_start}-{pam_start+2}; break GG at pos {target}",
    ))


def _add_suggestion_rev(suggestions, seen, template, cc_start, strand):
    target = cc_start
    if target in seen:
        return
    seen.add(target)
    ref_base = template[target]
    alt_base = "T" if ref_base == "C" else "A"
    suggestions.append(SuggestedMutation(
        pos=target, ref=ref_base, alt=alt_base,
        rationale=f"CCN motif ({strand}) at {cc_start}-{cc_start+2}; break CC at pos {target}",
    ))


def disrupt_pam(right_arm: str) -> tuple[str, bool, str]:
    """Find NGG PAM in right_arm, apply silent disruption. Returns (modified_arm, was_disrupted, note)."""
    for i in range(len(right_arm) - 2):
        if right_arm[i + 1] == "G" and right_arm[i + 2] == "G":
            arm_list = list(right_arm)
            arm_list[i + 2] = "A"
            return "".join(arm_list), True, f"PAM disrupted at arm pos {i} (G->A)"
    return right_arm, False, "No NGG PAM found in right arm"
