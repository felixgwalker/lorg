"""
Layer 2: Genomic Completeness

Scores how completely the extinct genome has been recovered relative to the
proxy reference. Three components are evaluated:
  - Overall breadth   (fraction of genome at ≥5× depth)
  - Coding coverage   (fraction of annotated coding regions recovered)
  - Regulatory coverage (fraction of annotated regulatory elements recovered)

Score range: 0–100 (higher = more complete recovery / more feasible).
"""

from ..config import GC_WEIGHTS


def score_gc(metrics: dict) -> dict:
    """
    Compute the Genomic Completeness layer score.

    Expected metrics keys:
        coverage_fraction_5x                (float) Overall genome fraction at ≥5×
        coding_regions_covered_fraction     (float) Coding regions fraction recovered
        regulatory_regions_covered_fraction (float) Regulatory regions fraction recovered
    """
    components = {
        "overall_breadth":     round(100.0 * metrics["coverage_fraction_5x"], 1),
        "coding_coverage":     round(100.0 * metrics["coding_regions_covered_fraction"], 1),
        "regulatory_coverage": round(100.0 * metrics["regulatory_regions_covered_fraction"], 1),
    }

    score = round(sum(components[k] * GC_WEIGHTS[k] for k in components), 1)

    flags: list[str] = []
    if metrics["coverage_fraction_5x"] < 0.80:
        flags.append("LOW_OVERALL_COVERAGE_BREADTH")
    if metrics["coding_regions_covered_fraction"] < 0.90:
        flags.append("INCOMPLETE_CODING_RECOVERY")
    if metrics["regulatory_regions_covered_fraction"] < 0.80:
        flags.append("INCOMPLETE_REGULATORY_RECOVERY")

    return {
        "score": score,
        "grade": _grade(score),
        "components": components,
        "interpretation": _interpret(score, metrics),
        "flags": flags,
    }


def _interpret(score: float, m: dict) -> str:
    breadth_pct   = m["coverage_fraction_5x"] * 100
    coding_pct    = m["coding_regions_covered_fraction"] * 100
    reg_pct       = m["regulatory_regions_covered_fraction"] * 100

    if score >= 85:
        return (
            f"Near-complete genome recovery. {breadth_pct:.1f}% of the extinct genome "
            f"is covered at ≥5× depth. Coding ({coding_pct:.1f}%) and regulatory "
            f"({reg_pct:.1f}%) regions are well-represented, enabling confident "
            f"divergence analysis across all functional partitions."
        )
    elif score >= 70:
        return (
            f"Good genomic completeness. Overall coverage breadth of {breadth_pct:.1f}% "
            f"supports most analyses. Gaps in regulatory regions ({reg_pct:.1f}%) may "
            f"introduce uncertainty in regulatory divergence estimates."
        )
    elif score >= 50:
        return (
            f"Partial genome recovery ({breadth_pct:.1f}% overall). Missing regions "
            f"limit confidence in whole-genome divergence estimates. Targeted enrichment "
            f"of coding ({coding_pct:.1f}%) and regulatory ({reg_pct:.1f}%) regions "
            f"may mitigate this."
        )
    else:
        return (
            f"Incomplete genome recovery ({breadth_pct:.1f}% overall). Substantial gaps "
            f"in coding and regulatory coverage will propagate errors to all downstream layers."
        )


def _grade(score: float) -> str:
    if score >= 80: return "A"
    if score >= 65: return "B"
    if score >= 50: return "C"
    if score >= 35: return "D"
    return "F"
