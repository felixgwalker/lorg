"""Introduce silent PAM-breaking mutations in right arm."""

import re


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
