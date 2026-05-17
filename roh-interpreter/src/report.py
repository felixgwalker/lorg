"""Write ROH BED catalog and FROH CSV."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.froh_calculator import FROHResult
from src.roh_detector import ROHSegment


def write_roh_bed(segments: list[ROHSegment], output_dir: Path) -> Path:
    """Write ROH segments as a BED file."""
    rows = [
        {
            "chrom": s.chrom,
            "start": s.start_pos,
            "end": s.end_pos,
            "name": f"{s.individual_id}_{s.length_class}",
            "score": int(s.mean_homozygosity * 1000),
            "strand": ".",
            "individual_id": s.individual_id,
            "n_snps": s.n_snps,
            "length_bp": s.length_bp,
            "length_class": s.length_class,
            "mean_homozygosity": s.mean_homozygosity,
        }
        for s in segments
    ]
    df = pd.DataFrame(rows)
    out = output_dir / "roh_catalog.bed"
    df.to_csv(out, sep="\t", index=False)
    return out


def write_froh_csv(froh_results: list[FROHResult], output_dir: Path) -> Path:
    """Write FROH per-individual summary CSV."""
    rows = [
        {
            "individual_id": fr.individual_id,
            "n_roh": fr.n_roh,
            "total_roh_bp": fr.total_roh_bp,
            "froh": fr.froh,
            "n_short": fr.n_short,
            "n_medium": fr.n_medium,
            "n_long": fr.n_long,
            "bp_short": fr.bp_short,
            "bp_medium": fr.bp_medium,
            "bp_long": fr.bp_long,
        }
        for fr in froh_results
    ]
    df = pd.DataFrame(rows)
    out = output_dir / "froh_summary.csv"
    df.to_csv(out, index=False)
    return out
