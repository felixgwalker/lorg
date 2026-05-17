"""Apply SNV/insertion/deletion edits to sequences."""


def apply_edit(seq: str, pos: int, edit_type: str, ref_allele: str, alt_allele: str) -> tuple[str, int]:
    """
    Apply an edit and return (new_seq, new_cut_pos).
    edit_type: 'snv', 'ins', 'del'
    pos: 0-based position of the edit in seq.
    """
    edit_type = edit_type.lower()
    if edit_type == "snv":
        if seq[pos].upper() != ref_allele.upper():
            raise ValueError(
                f"Reference mismatch at pos {pos}: expected {ref_allele}, found {seq[pos]}"
            )
        new_seq = seq[:pos] + alt_allele.upper() + seq[pos + 1:]
        return new_seq, pos
    elif edit_type == "ins":
        new_seq = seq[:pos] + alt_allele.upper() + seq[pos:]
        return new_seq, pos + len(alt_allele)
    elif edit_type == "del":
        del_len = len(ref_allele)
        new_seq = seq[:pos] + seq[pos + del_len:]
        return new_seq, pos
    else:
        raise ValueError(f"Unknown edit type: {edit_type}")


def build_template(ref_seq: str, cut_pos: int, left_arm_len: int, right_arm_len: int,
                   edit_type: str, ref_allele: str, alt_allele: str) -> dict:
    """Extract arms around cut site, apply edit, return template components."""
    left_start = max(0, cut_pos - left_arm_len)
    right_end = min(len(ref_seq), cut_pos + right_arm_len)

    region = ref_seq[left_start:right_end]
    local_cut = cut_pos - left_start

    edited_region, new_local_cut = apply_edit(region, local_cut, edit_type, ref_allele, alt_allele)

    left_arm = edited_region[:new_local_cut]
    right_arm = edited_region[new_local_cut + (1 if edit_type == "snv" else 0):]
    edit_seq = alt_allele.upper()

    return {
        "left_arm": left_arm,
        "edit_seq": edit_seq,
        "right_arm": right_arm,
        "left_arm_len": len(left_arm),
        "right_arm_len": len(right_arm),
        "full_template": left_arm + edit_seq + right_arm,
        "cut_pos": cut_pos,
        "left_start": left_start,
    }
