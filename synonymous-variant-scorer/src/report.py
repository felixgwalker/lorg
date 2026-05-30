"""Write per-mechanism score CSV and composite impact index CSV.

Output files
------------
mechanism_scores.csv
    One row per variant with all four mechanism scores.
composite_impact.csv
    One row per variant with the composite score and impact tier, sorted
    descending by composite score.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.aggregator import ScoredVariant


# Column order for mechanism CSV
_MECH_COLS = [
    "variant_id",
    "chrom",
    "pos",
    "ref_codon",
    "alt_codon",
    "gene",
    "transcript",
    "splicing_score",
    "codon_usage_score",
    "mrna_stability_score",
    "folding_score",
    "composite_score",
    "impact_tier",
]

# Column order for composite CSV
_COMP_COLS = [
    "variant_id",
    "chrom",
    "pos",
    "ref_codon",
    "alt_codon",
    "gene",
    "transcript",
    "composite_score",
    "impact_tier",
    "splicing_score",
    "codon_usage_score",
    "mrna_stability_score",
    "folding_score",
]


def write_mechanism_csv(scored: list[ScoredVariant], output_dir: Path) -> Path:
    """Write per-mechanism scores CSV.

    Includes all four mechanism scores plus the composite index and impact
    tier for convenience.

    Parameters
    ----------
    scored:
        List of ``ScoredVariant`` dataclass instances.
    output_dir:
        Directory where the file will be written.

    Returns
    -------
    Path
        Absolute path to the written CSV file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = [
        {
            "variant_id": s.variant_id,
            "chrom": s.chrom,
            "pos": s.pos,
            "ref_codon": s.ref_codon,
            "alt_codon": s.alt_codon,
            "gene": s.gene,
            "transcript": s.transcript,
            "splicing_score": s.splicing_score,
            "codon_usage_score": s.codon_usage_score,
            "mrna_stability_score": s.mrna_stability_score,
            "folding_score": s.folding_score,
            "composite_score": s.composite_score,
            "impact_tier": s.impact_tier,
        }
        for s in scored
    ]
    df = pd.DataFrame(rows, columns=_MECH_COLS)
    out = output_dir / "mechanism_scores.csv"
    df.to_csv(out, index=False)
    return out


def write_composite_csv(scored: list[ScoredVariant], output_dir: Path) -> Path:
    """Write composite impact index CSV, sorted by descending composite score.

    Parameters
    ----------
    scored:
        List of ``ScoredVariant`` dataclass instances.
    output_dir:
        Directory where the file will be written.

    Returns
    -------
    Path
        Absolute path to the written CSV file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
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
            "splicing_score": s.splicing_score,
            "codon_usage_score": s.codon_usage_score,
            "mrna_stability_score": s.mrna_stability_score,
            "folding_score": s.folding_score,
        }
        for s in scored
    ]
    df = pd.DataFrame(rows, columns=_COMP_COLS)
    df.sort_values("composite_score", ascending=False, inplace=True)
    out = output_dir / "composite_impact.csv"
    df.to_csv(out, index=False)
    return out
