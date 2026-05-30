"""Write HDR template outputs: FASTA, annotated sequence, QC CSV, arm variants CSV."""

from pathlib import Path
import pandas as pd


def write_template_fasta(template: dict, output_dir: Path, name: str = "hdr_template") -> Path:
    out_path = output_dir / "hdr_template.fa"
    full_seq = template["full_template"]
    with open(out_path, "w") as fh:
        fh.write(f">{name} left_arm={template['left_arm_len']}bp edit={template['edit_seq']} "
                 f"right_arm={template['right_arm_len']}bp\n")
        for i in range(0, len(full_seq), 60):
            fh.write(full_seq[i:i + 60] + "\n")
    return out_path


def write_annotated_sequence(
    template: dict,
    output_dir: Path,
    pam_mutations: list[dict] | None = None,
    line_width: int = 60,
) -> Path:
    """Write an annotated plain-text view of the HDR template sequence.

    The output file shows the template sequence with a ruler, a boundary
    annotation line indicating where the left arm ends, where the edit sits,
    and where the right arm begins.  Optional PAM-disruption mutations are
    flagged with caret (^) markers below the sequence.

    Parameters
    ----------
    template : dict
        Template dict as returned by ``build_template``, with keys
        ``left_arm``, ``edit_seq``, ``right_arm``, ``full_template``,
        ``left_arm_len``, ``right_arm_len``.
    output_dir : Path
        Directory to write ``annotated_template.txt`` into.
    pam_mutations : list[dict] | None
        Optional list of PAM-disruption mutations (from
        :func:`~src.pam_disruptor.suggest_pam_disruption`), each with keys
        ``pos``, ``ref``, ``alt``, ``rationale``.
    line_width : int
        Number of sequence bases per line (default 60).

    Returns
    -------
    Path
        Path to the written file.
    """
    out_path = output_dir / "annotated_template.txt"

    left_arm = template["left_arm"]
    edit_seq = template["edit_seq"]
    right_arm = template["right_arm"]
    full_seq = template["full_template"]
    left_len = len(left_arm)
    edit_len = len(edit_seq)

    # Build annotation characters for each position:
    #   L = left arm, E = edit, R = right arm
    ann_chars = ["L"] * left_len + ["E"] * edit_len + ["R"] * len(right_arm)

    # Build caret mask for PAM mutation positions.
    caret_mask: dict[int, str] = {}
    if pam_mutations:
        for mut in pam_mutations:
            caret_mask[mut["pos"]] = "^"

    lines: list[str] = []
    lines.append("=" * 72)
    lines.append("HDR Donor Template — Annotated Sequence")
    lines.append("=" * 72)
    lines.append(
        f"Total length : {len(full_seq)} bp"
    )
    lines.append(
        f"Left arm     : positions 1–{left_len} ({left_len} bp)"
    )
    edit_start = left_len + 1
    edit_end = left_len + edit_len
    lines.append(
        f"Edit         : positions {edit_start}–{edit_end} ({edit_len} bp)  [{edit_seq}]"
    )
    right_start = edit_end + 1
    lines.append(
        f"Right arm    : positions {right_start}–{len(full_seq)} ({len(right_arm)} bp)"
    )
    lines.append("")

    # Legend
    lines.append("Annotation key:  L = left arm   E = edit   R = right arm")
    if pam_mutations:
        lines.append("                 ^ = PAM-disruption mutation site")
    lines.append("")

    # Sequence blocks
    for block_start in range(0, len(full_seq), line_width):
        block_end = min(block_start + line_width, len(full_seq))
        seq_block = full_seq[block_start:block_end]
        ann_block = "".join(ann_chars[block_start:block_end])

        # Ruler: 1-based positions at the start and end of the block.
        ruler = f"{block_start + 1:>6}  "
        lines.append(ruler + seq_block)
        lines.append(" " * len(ruler) + ann_block)

        # Caret line (only if there are carets in this block).
        caret_line = ""
        has_caret = False
        for i in range(block_start, block_end):
            if i in caret_mask:
                caret_line += caret_mask[i]
                has_caret = True
            else:
                caret_line += " "
        if has_caret:
            lines.append(" " * len(ruler) + caret_line)

        lines.append("")

    # PAM mutation annotations
    if pam_mutations:
        lines.append("-" * 72)
        lines.append("PAM-disruption mutations:")
        for mut in pam_mutations:
            lines.append(
                f"  pos {mut['pos'] + 1:>5} (0-based {mut['pos']:>5})  "
                f"{mut['ref']}->{mut['alt']}  {mut['rationale']}"
            )
        lines.append("")

    lines.append("=" * 72)

    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    return out_path


def write_qc_report(qc_checks: list[dict], output_dir: Path) -> Path:
    out_path = output_dir / "qc_report.csv"
    df = pd.DataFrame(qc_checks)
    df.to_csv(out_path, index=False)
    return out_path


def write_arm_variants(variants: list[dict], output_dir: Path) -> Path:
    out_path = output_dir / "arm_variants.csv"
    df = pd.DataFrame(variants)
    df.to_csv(out_path, index=False)
    return out_path
