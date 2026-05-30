"""Introduce silent PAM-breaking mutations in right arm."""

import re
from typing import TypedDict


class SuggestedMutation(TypedDict):
    """A proposed single-base change that disrupts an NGG PAM site."""
    pos: int          # 0-based position in the full template sequence
    ref: str          # reference base
    alt: str          # proposed replacement base
    rationale: str    # human-readable explanation


_PAM_SCAN_RADIUS = 10  # bp on each side of cut to scan for PAM sites


_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def suggest_pam_disruption(
    template: str,
    pam_seq: str,
    cut_pos: int,
) -> list[SuggestedMutation]:
    """Scan for NGG PAM sites within ±10 bp of *cut_pos* and propose disrupting mutations.

    Both the forward strand and reverse complement are searched, which covers all
    possible guide-RNA orientations.  For each PAM found, the function proposes
    the minimal single-base change that destroys the GG dinucleotide of the NGG
    motif, preferring a transversion at the second G (pos+2 of match) over the
    first G (pos+1) so as to minimise unintended amino-acid changes.

    Parameters
    ----------
    template : str
        The full HDR donor template sequence (uppercase).
    pam_seq : str
        PAM motif to search for.  Only ``"NGG"`` is currently meaningful; the
        ``N`` is treated as a wildcard.
    cut_pos : int
        0-based index of the cut site within *template*.

    Returns
    -------
    list[SuggestedMutation]
        One dict per discovered PAM site, sorted by proximity to *cut_pos*.
        Empty list when no PAM sites are found in the scan window.
    """
    # Determine window bounds for the scan.
    window_start = max(0, cut_pos - _PAM_SCAN_RADIUS)
    window_end = min(len(template), cut_pos + _PAM_SCAN_RADIUS + 3)  # +3 to allow full PAM at boundary

    suggestions: list[SuggestedMutation] = []
    seen_positions: set[int] = set()

    # ---- Forward strand scan ------------------------------------------------
    for i in range(window_start, window_end - 2):
        if template[i + 1] == "G" and template[i + 2] == "G":
            # NGG found at i..i+2
            # Prefer mutating the second G (i+2) G->A (transversion, disrupts GG)
            _add_suggestion(suggestions, seen_positions, template, i, strand="+")

    # ---- Reverse complement scan --------------------------------------------
    # On the reverse strand an NGG PAM appears as CCN on the forward strand.
    for i in range(window_start, window_end - 2):
        if template[i] == "C" and template[i + 1] == "C":
            # CCN on forward = NGG on reverse at this region.
            # To disrupt the PAM on the reverse strand we mutate the first C
            # on the forward strand (which is the second G on the reverse).
            _add_suggestion_rev(suggestions, seen_positions, template, i, strand="-")

    # Sort by proximity to cut site.
    suggestions.sort(key=lambda s: abs(s["pos"] - cut_pos))
    return suggestions


def _add_suggestion(
    suggestions: list,
    seen: set,
    template: str,
    pam_start: int,
    strand: str,
) -> None:
    """Propose a mutation at the second G of an NGG PAM (forward strand)."""
    # Mutate position pam_start+2 (the second G) to A.
    target = pam_start + 2
    if target in seen:
        return
    seen.add(target)
    ref_base = template[target]
    alt_base = "A" if ref_base == "G" else "G"
    suggestions.append(SuggestedMutation(
        pos=target,
        ref=ref_base,
        alt=alt_base,
        rationale=(
            f"NGG PAM ({strand} strand) at template positions "
            f"{pam_start}-{pam_start + 2}; mutate pos {target} "
            f"{ref_base}->{alt_base} to break GG dinucleotide"
        ),
    ))


def _add_suggestion_rev(
    suggestions: list,
    seen: set,
    template: str,
    cc_start: int,
    strand: str,
) -> None:
    """Propose a mutation that disrupts a CCN motif on forward (= NGG rev-complement PAM)."""
    # Mutate the first C (cc_start) to T to break CC.
    target = cc_start
    if target in seen:
        return
    seen.add(target)
    ref_base = template[target]
    alt_base = "T" if ref_base == "C" else "A"
    suggestions.append(SuggestedMutation(
        pos=target,
        ref=ref_base,
        alt=alt_base,
        rationale=(
            f"CCN motif ({strand} strand / NGG rev-comp PAM) at template positions "
            f"{cc_start}-{cc_start + 2}; mutate pos {target} "
            f"{ref_base}->{alt_base} to break CC dinucleotide"
        ),
    ))


_CODON_TABLE = {
    "TTT": "TTC", "TTC": "TTT", "TTA": "TTG", "TTG": "TTA",
    "CTT": "CTC", "CTC": "CTT", "CTA": "CTG", "CTG": "CTA",
    "ATT": "ATC", "ATC": "ATT", "ATA": "ATC", "ATG": "ATG",
    "GTT": "GTC", "GTC": "GTT", "GTA": "GTG", "GTG": "GTA",
    "TCT": "TCC", "TCC": "TCT", "TCA": "TCG", "TCG": "TCA",
    "CCT": "CCC", "CCC": "CCT", "CCA": "CCG", "CCG": "CCA",
    "ACT": "ACC", "ACC": "ACT", "ACA": "ACG", "ACG": "ACA",
    "GCT": "GCC", "GCC": "GCT", "GCA": "GCG", "GCG": "GCA",
    "TAT": "TAC", "TAC": "TAT", "TAA": "TAG", "TAG": "TAA",
    "CAT": "CAC", "CAC": "CAT", "CAA": "CAG", "CAG": "CAA",
    "AAT": "AAC", "AAC": "AAT", "AAA": "AAG", "AAG": "AAA",
    "GAT": "GAC", "GAC": "GAT", "GAA": "GAG", "GAG": "GAA",
    "TGT": "TGC", "TGC": "TGT", "TGA": "TGA", "TGG": "TGG",
    "CGT": "CGC", "CGC": "CGT", "CGA": "CGG", "CGG": "CGA",
    "AGT": "AGC", "AGC": "AGT", "AGA": "AGG", "AGG": "AGA",
    "GGT": "GGC", "GGC": "GGT", "GGA": "GGG", "GGG": "GGA",
}


def find_pam_in_arm(arm: str) -> int:
    """Return index of first NGG PAM in arm, or -1."""
    for i in range(len(arm) - 2):
        if arm[i + 1] == "G" and arm[i + 2] == "G":
            return i
    return -1


def disrupt_pam(right_arm: str) -> tuple[str, bool, str]:
    """
    Find NGG PAM in right_arm, mutate 3rd position of nearest codon.
    Returns (modified_arm, was_disrupted, note).
    """
    pam_idx = find_pam_in_arm(right_arm)
    if pam_idx < 0:
        return right_arm, False, "No NGG PAM found in right arm"

    arm_list = list(right_arm)
    # Mutate the G at position pam_idx+1 to A to break GG
    target_pos = pam_idx + 1
    # Try synonymous change at 3rd codon position nearest to PAM
    codon_start = (target_pos // 3) * 3
    if codon_start + 2 < len(arm_list):
        codon = "".join(arm_list[codon_start:codon_start + 3])
        synonymous = _CODON_TABLE.get(codon)
        if synonymous and synonymous != codon:
            arm_list[codon_start:codon_start + 3] = list(synonymous)
            modified = "".join(arm_list)
            # Verify PAM is disrupted
            if modified[pam_idx + 1:pam_idx + 3] != "GG":
                return modified, True, f"PAM disrupted at arm pos {pam_idx} via synonymous codon change"

    # Fallback: directly mutate second G of PAM
    arm_list[pam_idx + 2] = "A"
    modified = "".join(arm_list)
    return modified, True, f"PAM disrupted at arm pos {pam_idx} (direct G->A substitution)"
