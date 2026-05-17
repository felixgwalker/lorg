"""
Layer 1: Ancient DNA Quality

Scores the input ancient DNA dataset on four components:
  - Coverage depth   (mean read depth across genome)
  - Coverage breadth (fraction of genome at ≥5× depth)
  - Fragment length  (mean read length; ancient DNA degrades into short fragments)
  - Contamination    (estimated fraction of reads from modern sources)

Score range: 0–100 (higher = better quality / more feasible).
"""

import math
from ..config import ADQ_WEIGHTS


def _score_coverage_depth(depth: float) -> float:
    """
    Saturating function: returns near-100 above ~50×, diminishes rapidly below 10×.
    Formula: 100 × (1 − e^(−depth/15))
    Reference points: 5× → 28, 10× → 49, 20× → 74, 30× → 86, 47× → 96.
    """
    return round(min(100.0, 100.0 * (1.0 - math.exp(-depth / 15.0))), 1)


def _score_coverage_breadth(fraction: float) -> float:
    """Linear mapping of fraction-of-genome-covered to 0–100."""
    return round(100.0 * max(0.0, min(1.0, fraction)), 1)


def _score_fragment_length(mean_bp: float) -> float:
    """
    Piecewise score reflecting ancient DNA preservation.
    Ancient fragments are typically 50–150 bp; shorter = poorer preservation.
    Longer fragments improve mappability and phasing confidence.
    """
    if mean_bp < 25:
        return 5.0
    elif mean_bp < 60:
        return round(5.0 + 50.0 * (mean_bp - 25) / 35, 1)
    elif mean_bp < 120:
        return round(55.0 + 30.0 * (mean_bp - 60) / 60, 1)
    elif mean_bp < 200:
        return round(85.0 + 10.0 * (mean_bp - 120) / 80, 1)
    else:
        return 95.0


def _score_contamination(rate: float) -> float:
    """
    Exponential penalty for modern contamination.
    Formula: 100 × e^(−6 × rate)
    Reference points: 0% → 100, 1% → 94, 2% → 89, 5% → 74, 10% → 55, 20% → 30.
    """
    return round(max(0.0, 100.0 * math.exp(-6.0 * rate)), 1)


def score_adq(metrics: dict) -> dict:
    """
    Compute the Ancient DNA Quality layer score.

    Expected metrics keys:
        mean_coverage_depth     (float) Mean read depth across the genome
        coverage_breadth_5x     (float) Fraction of genome covered at ≥5×
        mean_fragment_length_bp (float) Mean fragment length in base pairs
        contamination_rate      (float) Estimated contamination fraction (0–1)
        damage_ct_5prime        (float) 5′ C→T misincorporation rate (informational)
    """
    components = {
        "coverage_depth":   _score_coverage_depth(metrics["mean_coverage_depth"]),
        "coverage_breadth": _score_coverage_breadth(metrics["coverage_breadth_5x"]),
        "fragment_length":  _score_fragment_length(metrics["mean_fragment_length_bp"]),
        "contamination":    _score_contamination(metrics["contamination_rate"]),
    }

    score = round(sum(components[k] * ADQ_WEIGHTS[k] for k in components), 1)

    flags: list[str] = []
    if metrics["contamination_rate"] > 0.05:
        flags.append("HIGH_CONTAMINATION")
    if metrics["mean_coverage_depth"] < 10:
        flags.append("LOW_COVERAGE_DEPTH")
    if metrics["mean_fragment_length_bp"] < 40:
        flags.append("VERY_SHORT_FRAGMENTS")
    if metrics.get("damage_ct_5prime", 0) > 0.50:
        flags.append("SEVERE_DEAMINATION_DAMAGE")

    return {
        "score": score,
        "grade": _grade(score),
        "components": components,
        "interpretation": _interpret(score, metrics),
        "flags": flags,
    }


def _interpret(score: float, m: dict) -> str:
    depth = m["mean_coverage_depth"]
    contamination_pct = m["contamination_rate"] * 100
    frag = m["mean_fragment_length_bp"]
    if score >= 80:
        return (
            f"Excellent ancient DNA quality. Mean coverage of {depth:.1f}× "
            f"and {contamination_pct:.1f}% contamination support high-confidence "
            f"variant calling. Fragment length ({frag:.0f} bp) is within the "
            f"acceptable range for well-preserved ancient specimens."
        )
    elif score >= 60:
        return (
            f"Good ancient DNA quality. Coverage ({depth:.1f}×) and contamination "
            f"({contamination_pct:.1f}%) are within acceptable ranges. Fragment "
            f"length ({frag:.0f} bp) may limit mappability in repetitive regions."
        )
    elif score >= 40:
        return (
            f"Marginal ancient DNA quality. Significant limitations in at least "
            f"one metric: coverage {depth:.1f}×, contamination {contamination_pct:.1f}%, "
            f"mean fragment {frag:.0f} bp. Downstream analyses should be interpreted cautiously."
        )
    else:
        return (
            f"Poor ancient DNA quality. Coverage ({depth:.1f}×), contamination "
            f"({contamination_pct:.1f}%), or fragment preservation ({frag:.0f} bp) "
            f"fall below thresholds for reliable genome reconstruction."
        )


def _grade(score: float) -> str:
    if score >= 80: return "A"
    if score >= 65: return "B"
    if score >= 50: return "C"
    if score >= 35: return "D"
    return "F"
