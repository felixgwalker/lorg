"""Apply SNV/insertion/deletion edits to sequences.

Salvaged from hdr-template-designer (deleted stage1f).
"""

_MAX_EDIT_BP = 20


def apply_edit(seq: str, edit_spec: dict) -> str:
    """Apply a single edit defined by edit_spec and return the mutated sequence."""
    edit_type = edit_spec.get("type", "").upper()
    if edit_type == "SNV":
        edit_type = "SNP"
    ref_bases = str(edit_spec.get("ref_bases", "")).upper()
    alt_bases = str(edit_spec.get("alt_bases", "")).upper()
    if "abs_position" in edit_spec:
        pos = int(edit_spec["abs_position"])
    elif "cut_offset" in edit_spec:
        pos = int(edit_spec["cut_offset"]) + int(edit_spec.get("position", 0))
    else:
        pos = int(edit_spec.get("position", 0))
    if pos < 0 or pos > len(seq):
        raise ValueError(f"Edit position {pos} out of range for sequence length {len(seq)}")
    max_bases = max(len(ref_bases), len(alt_bases))
    if max_bases > _MAX_EDIT_BP:
        raise ValueError(f"Edit size {max_bases} bp exceeds maximum ({_MAX_EDIT_BP} bp)")
    if edit_type in ("SNP", "MNP"):
        ref_len = len(ref_bases) if ref_bases else 1
        if ref_bases:
            actual = seq[pos:pos + ref_len].upper()
            if actual != ref_bases:
                raise ValueError(f"Reference mismatch at pos {pos}: expected '{ref_bases}', found '{actual}'")
        if not alt_bases:
            raise ValueError("alt_bases must not be empty for SNP/MNP edits")
        return seq[:pos] + alt_bases + seq[pos + ref_len:]
    elif edit_type == "INSERTION":
        if not alt_bases:
            raise ValueError("alt_bases must not be empty for insertion edits")
        return seq[:pos] + alt_bases + seq[pos:]
    elif edit_type == "DELETION":
        del_len = len(ref_bases) if ref_bases else 1
        if ref_bases:
            actual = seq[pos:pos + del_len].upper()
            if actual != ref_bases:
                raise ValueError(f"Reference mismatch at pos {pos}: expected '{ref_bases}', found '{actual}'")
        return seq[:pos] + seq[pos + del_len:]
    else:
        raise ValueError(f"Unknown edit type '{edit_spec.get('type')}'. Supported: SNP, SNV, insertion, deletion")


def build_template(ref_seq: str, cut_pos: int, left_arm_len: int, right_arm_len: int,
                   edit_type: str, ref_allele: str, alt_allele: str) -> dict:
    """Extract arms around cut site, apply edit, return template components."""
    left_start = max(0, cut_pos - left_arm_len)
    right_end = min(len(ref_seq), cut_pos + right_arm_len)
    region = ref_seq[left_start:right_end]
    local_cut = cut_pos - left_start
    edit_type_l = edit_type.lower()
    if edit_type_l == "snv":
        edited = region[:local_cut] + alt_allele.upper() + region[local_cut + 1:]
        new_local_cut = local_cut
    elif edit_type_l == "ins":
        edited = region[:local_cut] + alt_allele.upper() + region[local_cut:]
        new_local_cut = local_cut + len(alt_allele)
    elif edit_type_l == "del":
        edited = region[:local_cut] + region[local_cut + len(ref_allele):]
        new_local_cut = local_cut
    else:
        raise ValueError(f"Unknown edit type: {edit_type}")
    left_arm = edited[:new_local_cut]
    right_arm = edited[new_local_cut + (1 if edit_type_l == "snv" else 0):]
    return {
        "left_arm": left_arm, "edit_seq": alt_allele.upper(), "right_arm": right_arm,
        "left_arm_len": len(left_arm), "right_arm_len": len(right_arm),
        "full_template": left_arm + alt_allele.upper() + right_arm,
        "cut_pos": cut_pos, "left_start": left_start,
    }
