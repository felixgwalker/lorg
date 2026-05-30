"""Apply SNV/insertion/deletion edits to sequences."""

_MAX_EDIT_BP = 20


def apply_edit(seq: str, edit_spec: dict) -> str:
    """Apply a single edit defined by *edit_spec* and return the mutated sequence.

    Parameters
    ----------
    seq : str
        Full template sequence (e.g. the region extracted around the cut site).
    edit_spec : dict
        Must contain:

        * ``type``      – one of ``"SNP"``, ``"insertion"``, or ``"deletion"``
          (case-insensitive; ``"SNV"`` is also accepted as a synonym for ``"SNP"``).
        * ``position``  – int, position of the edit **relative to the cut site**.
          Positive values are downstream (right of cut), negative are upstream
          (left of cut).  The caller must translate this to an absolute index in
          *seq* using the ``cut_offset`` key, or supply the already-resolved
          absolute index directly via ``abs_position``.  When ``abs_position`` is
          present it takes precedence over ``position``.
        * ``ref_bases`` – str, expected reference bases at the edit locus
          (used for validation; may be empty string for pure insertions).
        * ``alt_bases`` – str, replacement bases
          (may be empty string for pure deletions).

    Returns
    -------
    str
        Mutated sequence string.

    Raises
    ------
    ValueError
        On reference mismatch, out-of-range position, oversized edit (> 20 bp),
        or unknown edit type.
    """
    edit_type = edit_spec.get("type", "").upper()
    # Normalise synonyms
    if edit_type == "SNV":
        edit_type = "SNP"

    ref_bases = str(edit_spec.get("ref_bases", "")).upper()
    alt_bases = str(edit_spec.get("alt_bases", "")).upper()

    # Resolve absolute position in seq.
    if "abs_position" in edit_spec:
        pos = int(edit_spec["abs_position"])
    elif "cut_offset" in edit_spec:
        pos = int(edit_spec["cut_offset"]) + int(edit_spec.get("position", 0))
    else:
        # Treat 'position' as already absolute (legacy / simple callers).
        pos = int(edit_spec.get("position", 0))

    if pos < 0 or pos > len(seq):
        raise ValueError(f"Edit position {pos} is out of range for sequence of length {len(seq)}")

    # Validate edit size.
    max_bases = max(len(ref_bases), len(alt_bases))
    if max_bases > _MAX_EDIT_BP:
        raise ValueError(
            f"Edit size {max_bases} bp exceeds maximum allowed ({_MAX_EDIT_BP} bp)"
        )

    if edit_type in ("SNP", "MNP"):
        # Single or multi-nucleotide substitution.
        ref_len = len(ref_bases) if ref_bases else 1
        if ref_bases:
            actual = seq[pos: pos + ref_len].upper()
            if actual != ref_bases:
                raise ValueError(
                    f"Reference mismatch at pos {pos}: expected '{ref_bases}', found '{actual}'"
                )
        if not alt_bases:
            raise ValueError("alt_bases must not be empty for SNP/MNP edits")
        new_seq = seq[:pos] + alt_bases + seq[pos + ref_len:]
        return new_seq

    elif edit_type == "INSERTION":
        # Insert alt_bases *before* pos (i.e. between pos-1 and pos).
        if not alt_bases:
            raise ValueError("alt_bases must not be empty for insertion edits")
        new_seq = seq[:pos] + alt_bases + seq[pos:]
        return new_seq

    elif edit_type == "DELETION":
        del_len = len(ref_bases) if ref_bases else 1
        if ref_bases:
            actual = seq[pos: pos + del_len].upper()
            if actual != ref_bases:
                raise ValueError(
                    f"Reference mismatch at pos {pos}: expected '{ref_bases}', found '{actual}'"
                )
        new_seq = seq[:pos] + seq[pos + del_len:]
        return new_seq

    else:
        raise ValueError(
            f"Unknown edit type '{edit_spec.get('type')}'. "
            "Supported: SNP, SNV, insertion, deletion"
        )


def _apply_edit_positional(seq: str, pos: int, edit_type: str, ref_allele: str, alt_allele: str) -> tuple[str, int]:
    """Internal helper: positional API used by build_template.

    Returns (new_seq, new_cut_pos).
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

    edited_region, new_local_cut = _apply_edit_positional(region, local_cut, edit_type, ref_allele, alt_allele)

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
