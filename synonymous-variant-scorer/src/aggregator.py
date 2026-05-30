"""Weighted composite score aggregator for the Synonymous Variant Scorer.

Weights (per specification):
  - Splicing disruption : 0.30
  - mRNA stability      : 0.25
  - Codon usage bias    : 0.30
  - Co-translational folding : 0.15

Impact tiers:
  - benign           : composite < 0.2
  - uncertain        : 0.2 <= composite <= 0.5
  - likely_functional: composite > 0.5
"""

from __future__ import annotations

from dataclasses import dataclass


MECHANISM_WEIGHTS: dict[str, float] = {
    "splicing":    0.30,
    "codon_usage": 0.30,
    "mrna_stability": 0.25,
    "folding":     0.15,
}

# Sanity check weights sum to 1.0
_weight_sum = sum(MECHANISM_WEIGHTS.values())
assert abs(_weight_sum - 1.0) < 1e-9, (
    f"Mechanism weights sum to {_weight_sum}, expected 1.0"
)


@dataclass
class ScoredVariant:
    """Per-variant mechanism scores and composite functional impact index."""
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
    """Compute weighted composite score and assign impact tier.

    Parameters
    ----------
    variant_id, chrom, pos, ref_codon, alt_codon, gene, transcript:
        Variant identity fields.
    splicing:
        Splicing disruption score in [0, 1].
    codon_usage:
        Codon usage bias score in [0, 1].
    mrna_stability:
        mRNA stability impact score in [0, 1].
    folding:
        Co-translational folding risk score in [0, 1].

    Returns
    -------
    ScoredVariant
        Dataclass with all component scores and the composite index.
    """
    composite = (
        MECHANISM_WEIGHTS["splicing"]      * splicing
        + MECHANISM_WEIGHTS["codon_usage"] * codon_usage
        + MECHANISM_WEIGHTS["mrna_stability"] * mrna_stability
        + MECHANISM_WEIGHTS["folding"]     * folding
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
    """Assign an impact tier based on the composite score.

    Thresholds (per specification):
      - benign           : score < 0.2
      - uncertain        : 0.2 <= score <= 0.5
      - likely_functional: score > 0.5
    """
    if score > 0.5:
        return "likely_functional"
    if score >= 0.2:
        return "uncertain"
    return "benign"
