"""Write BED file and CSV report for ERV loci."""

from pathlib import Path
import pandas as pd


def write_bed(hits: list[dict], output_dir: Path) -> Path:
    out_path = output_dir / "erv_loci.bed"
    with open(out_path, "w") as fh:
        fh.write("# track name=ERV_loci\n")
        for h in hits:
            # BED: chrom, start, end, name, score, strand
            name = f"{h['family']}|{h['ltr_type']}|{h['risk_tier']}"
            score = min(1000, h["risk_score"] * 125)
            fh.write(f"{h['chrom']}\t{h['start']}\t{h['end']}\t{name}\t{score}\t.\n")
    return out_path


def write_csv_report(hits: list[dict], output_dir: Path) -> Path:
    out_path = output_dir / "erv_risk_scores.csv"
    cols = [
        "chrom", "start", "end", "family", "ltr_type", "kmer_hits",
        "has_poly_a", "longest_orf_bp", "gc_content",
        "score_orf", "score_ltr", "score_gc", "risk_score", "risk_tier",
        "age_mya",
    ]
    df = pd.DataFrame(hits)
    df = df[[c for c in cols if c in df.columns]]
    df = df.sort_values("risk_score", ascending=False).reset_index(drop=True)
    df.to_csv(out_path, index=False)
    return out_path
