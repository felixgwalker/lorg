"""
Layer 4: Edit Burden Estimation

Estimates the practical genome-editing workload required to convert the proxy
genome toward the extinct genome, accounting for SNPs, indels, and structural
variants across functional regions.

Scoring uses a log-scale penalty calibrated to current and near-term CRISPR
throughput capacity (~10,000 targeted edits per project cycle):
  ~1,000 edits  → ~86  (within current capability)
  ~10,000 edits → ~72  (near-term feasible)
  ~100,000 edits → ~58  (decade-scale horizon)
  ~1,000,000 edits → ~44  (requires transformative technology advances)
  ~10,000,000 edits → ~30 (currently infeasible at scale)

Score range: 0–100 (higher = lower edit burden / more feasible).
"""

import math
from ..config import EB_WEIGHTS


def _count_score(n: int, scale: float) -> float:
    """
    Log-scale scoring: score = max(0, 100 − scale × log10(n)).
    scale controls how steeply the score decays with edit count.
    """
    if n <= 0:
        return 100.0
    return round(max(0.0, 100.0 - scale * math.log10(n)), 1)


def score_edit_burden(metrics: dict) -> dict:
    """
    Compute the Edit Burden layer score.

    Expected metrics keys:
        total_minimum_edits  (int)   Minimum functional edits required
        coding_edits         (int)   Edits in coding regions
        regulatory_edits     (int)   Edits in regulatory elements
        functional_svs       (int)   Structural variants requiring resolution
        crispr_efficiency_estimate  (float) Estimated on-target CRISPR efficiency
        current_max_crispr_throughput (int) Current max edits per project cycle
    """
    total    = metrics["total_minimum_edits"]
    coding   = metrics["coding_edits"]
    reg      = metrics["regulatory_edits"]
    svs      = metrics["functional_svs"]
    efficiency = metrics.get("crispr_efficiency_estimate", 1.0)

    # Effective counts accounting for CRISPR efficiency (more attempts needed)
    effective_total  = int(total  / efficiency) if efficiency > 0 else total
    effective_coding = int(coding / efficiency) if efficiency > 0 else coding
    effective_reg    = int(reg    / efficiency) if efficiency > 0 else reg

    components = {
        "total_edit_count":       _count_score(effective_total,  14.0),
        "coding_edit_burden":     _count_score(effective_coding, 16.0),
        "regulatory_edit_burden": _count_score(effective_reg,    16.0),
        "structural_complexity":  _count_score(svs,              12.0),
    }

    score = round(sum(components[k] * EB_WEIGHTS[k] for k in components), 1)

    throughput = metrics.get("current_max_crispr_throughput", 10000)
    years_at_current = round(total / throughput, 1)

    flags: list[str] = []
    if total > throughput * 10:
        flags.append("EDIT_COUNT_EXCEEDS_CURRENT_CAPABILITY")
    if coding > 50000:
        flags.append("HIGH_CODING_EDIT_BURDEN")
    if reg > 100000:
        flags.append("HIGH_REGULATORY_EDIT_BURDEN")
    if svs > 2000:
        flags.append("COMPLEX_STRUCTURAL_VARIANT_LOAD")

    return {
        "score": score,
        "grade": _grade(score),
        "components": components,
        "interpretation": _interpret(score, metrics, years_at_current),
        "flags": flags,
        "context": {
            "total_minimum_edits":  total,
            "coding_edits":         coding,
            "regulatory_edits":     reg,
            "functional_svs":       svs,
            "years_at_current_throughput": years_at_current,
            "crispr_efficiency":    efficiency,
        },
    }


def _interpret(score: float, m: dict, years: float) -> str:
    total = m["total_minimum_edits"]
    throughput = m.get("current_max_crispr_throughput", 10000)
    coding = m["coding_edits"]
    reg = m["regulatory_edits"]
    svs = m["functional_svs"]

    if score >= 70:
        return (
            f"Tractable edit burden. {total:,} functional edits required — "
            f"within or near current CRISPR throughput ({throughput:,}/project). "
            f"Coding ({coding:,}) and regulatory ({reg:,}) edits are manageable."
        )
    elif score >= 50:
        return (
            f"Significant edit burden. {total:,} edits required, approximately "
            f"{total // throughput}× current throughput. Achievable on a 5–15 year "
            f"horizon assuming continued improvement in multiplexed editing technologies."
        )
    elif score >= 35:
        return (
            f"Heavy edit burden. {total:,} functional edits ({coding:,} coding, "
            f"{reg:,} regulatory, {svs:,} SVs) would require ~{years:.0f} years at "
            f"current CRISPR capacity. Requires substantial advances in high-throughput "
            f"genome editing."
        )
    else:
        return (
            f"Prohibitive edit burden under current technology. {total:,} minimum "
            f"functional edits — {total // throughput}× current CRISPR throughput "
            f"({throughput:,}/project). At current rates: ~{years:.0f} years. "
            f"Feasibility requires an order-of-magnitude advance in genome editing "
            f"throughput, potentially via base editing arrays or synthetic chromosome assembly."
        )


def _grade(score: float) -> str:
    if score >= 80: return "A"
    if score >= 65: return "B"
    if score >= 50: return "C"
    if score >= 35: return "D"
    return "F"
