"""Write TSV/CSV outputs for guide RNA off-target scorer."""

from pathlib import Path
import pandas as pd


def write_off_target_table(hits: list[dict], output_dir: Path) -> Path:
    """Write all off-target hits ranked by CFD score to a TSV file.

    Columns: guide_name, guide_seq, chrom, pos, strand, target_seq,
             mismatches, mismatch_positions, pam, has_pam, cfd_score
    """
    out_path = output_dir / "off_target_sites.tsv"
    cols = [
        "guide_name", "guide_seq", "chrom", "pos", "strand",
        "target_seq", "mismatches", "mismatch_positions",
        "pam", "has_pam", "cfd_score",
    ]

    if not hits:
        df = pd.DataFrame(columns=cols)
    else:
        df = pd.DataFrame(hits)
        # Stringify mismatch_positions list for TSV readability
        if "mismatch_positions" in df.columns:
            df["mismatch_positions"] = df["mismatch_positions"].apply(
                lambda v: ";".join(str(x) for x in v) if isinstance(v, list) else str(v)
            )
        else:
            df["mismatch_positions"] = ""
        # Keep only known columns, in order
        df = df[[c for c in cols if c in df.columns]]
        df = df.sort_values("cfd_score", ascending=False).reset_index(drop=True)

    df.to_csv(out_path, index=False, sep="\t")
    return out_path


def write_specificity_summary(hits: list[dict], guides: list[tuple[str, str]], output_dir: Path) -> Path:
    """Write per-guide specificity summary to a TSV file.

    Columns: guide_name, guide_seq, on_target_sites, off_target_1mm,
             off_target_2mm, off_target_3mm, total_off_targets,
             max_off_target_cfd, specificity_score
    """
    out_path = output_dir / "specificity_summary.tsv"
    rows = []
    for guide_name, guide_seq in guides:
        guide_hits = [h for h in hits if h["guide_name"] == guide_name]
        on_target = [h for h in guide_hits if h["mismatches"] == 0 and h.get("has_pam")]
        off_1mm = [h for h in guide_hits if h["mismatches"] == 1]
        off_2mm = [h for h in guide_hits if h["mismatches"] == 2]
        off_3mm = [h for h in guide_hits if h["mismatches"] == 3]
        cfd_scores = [h["cfd_score"] for h in guide_hits if h["mismatches"] > 0]
        max_off_cfd = max(cfd_scores) if cfd_scores else 0.0
        rows.append({
            "guide_name": guide_name,
            "guide_seq": guide_seq,
            "on_target_sites": len(on_target),
            "off_target_1mm": len(off_1mm),
            "off_target_2mm": len(off_2mm),
            "off_target_3mm": len(off_3mm),
            "total_off_targets": len(off_1mm) + len(off_2mm) + len(off_3mm),
            "max_off_target_cfd": round(max_off_cfd, 4),
            "specificity_score": round(1.0 - max_off_cfd, 4),
        })
    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False, sep="\t")
    return out_path
