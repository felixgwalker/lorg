"""Write CSV outputs for guide RNA off-target scorer."""

from pathlib import Path
import pandas as pd


def write_off_target_table(hits: list[dict], output_dir: Path) -> Path:
    out_path = output_dir / "off_target_sites.csv"
    if not hits:
        df = pd.DataFrame(columns=[
            "guide_name", "guide_seq", "chrom", "pos", "strand",
            "target_seq", "mismatches", "pam", "has_pam", "cfd_score",
        ])
    else:
        df = pd.DataFrame(hits)
        cols = [
            "guide_name", "guide_seq", "chrom", "pos", "strand",
            "target_seq", "mismatches", "pam", "has_pam", "cfd_score",
        ]
        df = df[[c for c in cols if c in df.columns]]
        df = df.sort_values("cfd_score", ascending=False).reset_index(drop=True)
    df.to_csv(out_path, index=False)
    return out_path


def write_specificity_summary(hits: list[dict], guides: list[tuple[str, str]], output_dir: Path) -> Path:
    out_path = output_dir / "specificity_summary.csv"
    rows = []
    for guide_name, guide_seq in guides:
        guide_hits = [h for h in hits if h["guide_name"] == guide_name]
        on_target = [h for h in guide_hits if h["mismatches"] == 0 and h.get("has_pam")]
        off_1mm = [h for h in guide_hits if h["mismatches"] == 1]
        off_2mm = [h for h in guide_hits if h["mismatches"] == 2]
        off_3mm = [h for h in guide_hits if h["mismatches"] == 3]
        off_pam = [h for h in guide_hits if h["mismatches"] == 0 and h.get("has_pam")]
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
    df.to_csv(out_path, index=False)
    return out_path
