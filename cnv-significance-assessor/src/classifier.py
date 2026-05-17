"""
Rule-based CNV significance classifier.

Scoring model
─────────────
Each CNV receives an integer total score computed from four components:

  size_score        (0–3)   based on CNV size in bp
  gene_score        (0–3)   based on number of overlapping gene bodies
  dosage_score      (0–3)   based on max pHaplo/pLI (deletions) or pTriplo
                            (duplications) across overlapping genes
  pop_modifier      (0–3)   subtracted; based on population allele frequency

  total_score = size_score + gene_score + dosage_score − pop_modifier

Classification tiers (after overrides)
───────────────────────────────────────
  LIKELY_BENIGN       total_score ≤ SCORE_BENIGN_MAX  OR  pop_freq > POP_FREQ_CUTOFF
  VUS                 score between BENIGN and PATHOGENIC
  LIKELY_PATHOGENIC   total_score ≥ SCORE_PATHOGENIC_MIN

Override rules (applied after score thresholding):
  1. If pop_freq > DEFAULT_POP_FREQ_CUTOFF → force LIKELY_BENIGN regardless of score.
  2. If any overlapping gene has dosage sensitivity ≥ HIGH_DS_THRESHOLD AND
     pop_freq is absent or ≤ DEFAULT_POP_FREQ_CUTOFF → raise floor to VUS.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

from src.annotator import AnnotatedCNV
from src.config import (
    DEFAULT_POP_FREQ_CUTOFF,
    DOSAGE_SCORE_BREAKS,
    DOSAGE_SCORE_MAX,
    GENE_COUNT_SCORE_BREAKS,
    GENE_COUNT_SCORE_MAX,
    HIGH_DS_THRESHOLD,
    POP_FREQ_MOD_BREAKS,
    POP_FREQ_MOD_MAX,
    SCORE_BENIGN_MAX,
    SCORE_PATHOGENIC_MIN,
    SIZE_SCORE_BREAKS,
    SIZE_SCORE_MAX,
)

logger = logging.getLogger(__name__)

# Human-readable tier labels
TIER_BENIGN     = "LIKELY_BENIGN"
TIER_VUS        = "VUS"
TIER_PATHOGENIC = "LIKELY_PATHOGENIC"


@dataclass
class ScoredCNV:
    """A single CNV with its component scores and significance classification."""

    annotated: AnnotatedCNV

    size_score:   int
    gene_score:   int
    dosage_score: int
    pop_modifier: int
    total_score:  int

    significance_tier: str   # LIKELY_BENIGN | VUS | LIKELY_PATHOGENIC
    classification_reason: str


@dataclass
class ClassificationSummary:
    """Counts and fractions of CNVs in each significance tier."""

    n_likely_benign:     int
    n_vus:               int
    n_likely_pathogenic: int
    n_total:             int

    fraction_benign:     float
    fraction_vus:        float
    fraction_pathogenic: float


def classify_cnvs(
    annotated: list[AnnotatedCNV],
    pop_freq_cutoff: float = DEFAULT_POP_FREQ_CUTOFF,
) -> tuple[list[ScoredCNV], ClassificationSummary]:
    """
    Classify each annotated CNV and return per-CNV results and a summary.

    Args:
        annotated:       List from annotator.annotate_cnvs().
        pop_freq_cutoff: CNVs with pop_frequency > this → forced LIKELY_BENIGN.

    Returns:
        Tuple of (list[ScoredCNV], ClassificationSummary).
    """
    scored: list[ScoredCNV] = [
        _score_one(ann, pop_freq_cutoff) for ann in annotated
    ]

    n_benign  = sum(1 for s in scored if s.significance_tier == TIER_BENIGN)
    n_vus     = sum(1 for s in scored if s.significance_tier == TIER_VUS)
    n_path    = sum(1 for s in scored if s.significance_tier == TIER_PATHOGENIC)
    n_total   = len(scored)

    def _frac(n: int) -> float:
        return round(n / n_total, 4) if n_total > 0 else 0.0

    summary = ClassificationSummary(
        n_likely_benign=n_benign,
        n_vus=n_vus,
        n_likely_pathogenic=n_path,
        n_total=n_total,
        fraction_benign=_frac(n_benign),
        fraction_vus=_frac(n_vus),
        fraction_pathogenic=_frac(n_path),
    )

    logger.info(
        "Classification: %d LIKELY_BENIGN (%.1f%%), %d VUS (%.1f%%), "
        "%d LIKELY_PATHOGENIC (%.1f%%)  [total: %d]",
        n_benign, summary.fraction_benign * 100,
        n_vus,    summary.fraction_vus * 100,
        n_path,   summary.fraction_pathogenic * 100,
        n_total,
    )
    return scored, summary


# ── Per-CNV scoring ───────────────────────────────────────────────────────

def _score_one(ann: AnnotatedCNV, pop_freq_cutoff: float) -> ScoredCNV:
    """Compute component scores and apply classification logic for one CNV."""
    cnv = ann.record

    size_score   = _size_score(cnv.size)
    gene_score   = _gene_count_score(ann.n_genes)
    dosage_score = _dosage_score(ann.max_dosage_sensitivity)
    pop_modifier = _pop_modifier(ann.pop_frequency)
    total        = size_score + gene_score + dosage_score - pop_modifier

    pop_freq = ann.pop_frequency

    # ── Override 1: common variant → always LIKELY_BENIGN ────────────────
    if pop_freq is not None and pop_freq > pop_freq_cutoff:
        tier   = TIER_BENIGN
        reason = (
            f"Common variant: population frequency {pop_freq:.2%} exceeds "
            f"cutoff {pop_freq_cutoff:.2%}"
        )
        return ScoredCNV(
            annotated=ann,
            size_score=size_score, gene_score=gene_score,
            dosage_score=dosage_score, pop_modifier=pop_modifier,
            total_score=total,
            significance_tier=tier,
            classification_reason=reason,
        )

    # ── Score-based tier ─────────────────────────────────────────────────
    if total <= SCORE_BENIGN_MAX:
        tier   = TIER_BENIGN
        reason = f"Low aggregate score ({total}): size={size_score}, genes={gene_score}, dosage={dosage_score}, pop_mod=−{pop_modifier}"
    elif total >= SCORE_PATHOGENIC_MIN:
        tier   = TIER_PATHOGENIC
        reason = f"High aggregate score ({total}): size={size_score}, genes={gene_score}, dosage={dosage_score}, pop_mod=−{pop_modifier}"
    else:
        tier   = TIER_VUS
        reason = f"Intermediate score ({total}): size={size_score}, genes={gene_score}, dosage={dosage_score}, pop_mod=−{pop_modifier}"

    # ── Override 2: high-dosage-sensitivity gene raises floor to VUS ─────
    if (
        tier == TIER_BENIGN
        and ann.max_dosage_sensitivity >= HIGH_DS_THRESHOLD
        and (pop_freq is None or pop_freq <= pop_freq_cutoff)
    ):
        tier   = TIER_VUS
        reason += (
            f"; reclassified to VUS: overlapping gene has {ann.dosage_metric} "
            f"= {ann.max_dosage_sensitivity:.3f} ≥ {HIGH_DS_THRESHOLD}"
        )

    return ScoredCNV(
        annotated=ann,
        size_score=size_score, gene_score=gene_score,
        dosage_score=dosage_score, pop_modifier=pop_modifier,
        total_score=total,
        significance_tier=tier,
        classification_reason=reason,
    )


# ── Score component functions ─────────────────────────────────────────────

def _size_score(size_bp: int) -> int:
    """Map CNV size in bp to a 0–3 integer score."""
    for upper, score in SIZE_SCORE_BREAKS:
        if size_bp <= upper:
            return score
    return SIZE_SCORE_MAX


def _gene_count_score(n_genes: int) -> int:
    """Map number of overlapping genes to a 0–3 integer score."""
    for upper, score in GENE_COUNT_SCORE_BREAKS:
        if n_genes <= upper:
            return score
    return GENE_COUNT_SCORE_MAX


def _dosage_score(max_ds: float) -> int:
    """
    Map maximum dosage sensitivity score to a 0–3 integer score.

    The score is assigned to the interval [break_lower, break_upper) where
    break_lower is the lower bound of the matching DOSAGE_SCORE_BREAKS entry.
    """
    # DOSAGE_SCORE_BREAKS is ordered by lower bound; iterate to find the bucket.
    score = 0
    for lower, s in DOSAGE_SCORE_BREAKS:
        if max_ds >= lower:
            score = s
        else:
            break

    # Final bucket: values at or above the last breakpoint get the max score
    if max_ds >= 0.9:
        return DOSAGE_SCORE_MAX
    return score


def _pop_modifier(pop_freq: float | None) -> int:
    """
    Map population frequency to a modifier that is *subtracted* from the score.

    Higher frequency → larger subtraction.  Returns 0 if frequency is unknown.
    """
    if pop_freq is None:
        return 0

    modifier = 0
    for lower, mod in POP_FREQ_MOD_BREAKS:
        if pop_freq >= lower:
            modifier = mod
        else:
            break

    if pop_freq >= 0.01:
        return POP_FREQ_MOD_MAX
    return modifier
