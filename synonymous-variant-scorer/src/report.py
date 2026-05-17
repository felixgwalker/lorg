"""Write per-mechanism score CSV and composite impact index CSV."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.aggregator import ScoredVariant


def write_mechanism_csv(scored: list[ScoredVariant], output_dir: Path) -> Path:
    """Write per-mechanism scores CSV."""
    rows = [
        {
            "variant_id": s.variant_id,
            "chrom": s.chrom,
            "pos": s.pos,
            "ref_codon": s.ref_codon,
            "alt_codon": s.alt_codon,
            "gene": s.gene,
            "splicing_score": s.splicing_score,
            "codon_usage_score": s.codon_usage_score,
            "mrna_stability_score": s.mrna_stability_score,
            "folding_score": s.folding_score,
        }
        for s in scored
    ]
    df = pd.DataFrame(rows)
    out = output_dir / "mechanism_scores.csv"
    df.to_csv(out, index=False)
    return out


def write_composite_csv(scored: list[ScoredVariant], output_dir: Path) -> Path:
    """Write composite impact index CSV."""
    rows = [
        {
            "variant_id": s.variant_id,
            "chrom": s.chrom,
            "pos": s.pos,
            "ref_codon": s.ref_codon,
            "alt_codon": s.alt_codon,
            "gene": s.gene,
            "transcript": s.transcript,
            "composite_score": s.composite_score,
            "impact_tier": s.impact_tier,
        }
        for s in scored
    ]
    df = pd.DataFrame(rows)
    df.sort_values("composite_score", ascending=False, inplace=True)
    out = output_dir / "composite_impact.csv"
    df.to_csv(out, index=False)
    return out
