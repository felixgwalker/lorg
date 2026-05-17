"""Weighted composite score aggregator (equal weights across four mechanisms)."""

from __future__ import annotations

from dataclasses import dataclass


MECHANISM_WEIGHTS = {
    "splicing": 0.25,
    "codon_usage": 0.25,
    "mrna_stability": 0.25,
    "folding": 0.25,
}


@dataclass
class ScoredVariant:
    """Per-variant mechanism scores and composite impact index."""
    variant_id: str
    chrom: str
    pos: int
    ref_codon: str
    alt_codon: str
    gene: str
    transcript: str
    splicing_score: float
    codon_usage_score: float
    mrna_stability_score: float
    folding_score: float
    composite_score: float
    impact_tier: str


def aggregate_scores(
    variant_id: str,
    chrom: str,
    pos: int,
    ref_codon: str,
    alt_codon: str,
    gene: str,
    transcript: str,
    splicing: float,
    codon_usage: float,
    mrna_stability: float,
    folding: float,
) -> ScoredVariant:
    """Compute equal-weight composite score and assign impact tier."""
    composite = (
        MECHANISM_WEIGHTS["splicing"] * splicing
        + MECHANISM_WEIGHTS["codon_usage"] * codon_usage
        + MECHANISM_WEIGHTS["mrna_stability"] * mrna_stability
        + MECHANISM_WEIGHTS["folding"] * folding
    )
    composite = round(min(1.0, max(0.0, composite)), 4)
    tier = _assign_tier(composite)
    return ScoredVariant(
        variant_id=variant_id,
        chrom=chrom,
        pos=pos,
        ref_codon=ref_codon,
        alt_codon=alt_codon,
        gene=gene,
        transcript=transcript,
        splicing_score=round(splicing, 4),
        codon_usage_score=round(codon_usage, 4),
        mrna_stability_score=round(mrna_stability, 4),
        folding_score=round(folding, 4),
        composite_score=composite,
        impact_tier=tier,
    )


def _assign_tier(score: float) -> str:
    if score >= 0.6:
        return "HIGH"
    if score >= 0.3:
        return "MODERATE"
    return "LOW"
