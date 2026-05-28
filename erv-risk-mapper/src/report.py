"""Write BED file and CSV report for ERV loci."""

from pathlib import Path
import pandas as pd


def write_bed(hits: list[dict], output_dir: Path) -> Path:
    """Write a 6-column BED file (chrom start end name score strand).

    score is composite_risk scaled to the BED 0–1000 integer range.
    strand is '.' (unknown) for all ERV elements.
    """
    out_path = output_dir / "erv_loci.bed"
    with open(out_path, "w") as fh:
        fh.write("# track name=ERV_loci\n")
        for h in hits:
            name = f"{h['family']}|{h['ltr_type']}|{h['risk_tier']}"
            # BED score: integer 0–1000
            score = int(min(1000, round(h["composite_risk"] * 1000)))
            fh.write(
                f"{h['chrom']}\t{h['start']}\t{h['end']}\t{name}\t{score}\t.\n"
            )
    return out_path


def write_csv_report(hits: list[dict], output_dir: Path) -> Path:
    """Write a CSV report with all hit fields including composite_risk and risk_tier."""
    out_path = output_dir / "erv_risk_scores.csv"
    cols = [
        "chrom", "start", "end", "family", "ltr_type", "kmer_hits",
        "has_poly_a", "longest_orf_bp", "gc_content",
        "age_mya",
        "orf_score", "ltr_score", "age_score", "family_score",
        "composite_risk", "risk_tier",
    ]
    df = pd.DataFrame(hits)
    # Keep only columns that exist (graceful if a field is absent)
    df = df[[c for c in cols if c in df.columns]]
    df = df.sort_values("composite_risk", ascending=False).reset_index(drop=True)
    df.to_csv(out_path, index=False)
    return out_path
