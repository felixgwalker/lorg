"""
Layer 5: Ethical and Ecological Flagging

Generates structured flags and a composite score across four sub-dimensions:
  - Habitat availability and quality
  - Ecological role (is the niche currently unfilled?)
  - Welfare considerations for surrogate/resurrected animals
  - Regulatory and conservation conflicts

Unlike the technical layers, this layer incorporates qualitative flags. The
numeric score reflects the overall ethical/ecological *tractability* — how
straightforward a responsible reintroduction programme would be if technical
barriers were resolved.

Score range: 0–100 (higher = fewer barriers / more tractable).
"""

from ..config import EE_WEIGHTS

_HABITAT_QUALITY_SCORES: dict[str, float] = {
    "excellent": 95.0,
    "good":      75.0,
    "moderate":  55.0,
    "poor":      25.0,
    "none":       0.0,
}

_ROLE_SCORES: dict[str, float] = {
    "apex_predator": 75.0,
    "keystone":      80.0,
    "ecosystem_engineer": 70.0,
    "specialist":    60.0,
    "generalist":    55.0,
    "unknown":       45.0,
}


def score_ethical_ecological(metrics: dict) -> dict:
    """
    Compute the Ethical/Ecological layer score.

    Expected metrics keys:
        habitat_available               (bool)   Whether habitat exists
        habitat_quality                 (str)    "excellent"|"good"|"moderate"|"poor"|"none"
        ecological_role                 (str)    Role category
        ecological_role_currently_filled (bool)  Is niche already occupied?
        welfare_flags                   (list)   Welfare concern strings
        regulatory_flags                (list)   Regulatory conflict strings
        conservation_flags              (list)   Conservation conflict strings
    """
    # ── Habitat ──────────────────────────────────────────────────────────────
    if not metrics["habitat_available"]:
        habitat_score = 0.0
    else:
        habitat_score = _HABITAT_QUALITY_SCORES.get(
            metrics.get("habitat_quality", "unknown"), 40.0
        )

    # ── Ecological role ───────────────────────────────────────────────────────
    role_base = _ROLE_SCORES.get(metrics.get("ecological_role", "unknown"), 45.0)
    if metrics.get("ecological_role_currently_filled", True):
        # Niche is already filled — reduced justification for reintroduction
        ecological_score = role_base * 0.40
    else:
        ecological_score = role_base

    # ── Welfare ───────────────────────────────────────────────────────────────
    n_welfare = len(metrics.get("welfare_flags", []))
    welfare_score = max(0.0, 100.0 - 25.0 * n_welfare)

    # ── Regulatory / conservation ─────────────────────────────────────────────
    n_reg  = len(metrics.get("regulatory_flags", []))
    n_cons = len(metrics.get("conservation_flags", []))
    reg_cons_score = max(0.0, 100.0 - 20.0 * n_reg - 15.0 * n_cons)

    components = {
        "habitat":                 round(habitat_score, 1),
        "ecological_role":         round(ecological_score, 1),
        "welfare":                 round(welfare_score, 1),
        "regulatory_conservation": round(reg_cons_score, 1),
    }

    score = round(sum(components[k] * EE_WEIGHTS[k] for k in components), 1)

    all_flags = (
        metrics.get("welfare_flags", [])
        + metrics.get("regulatory_flags", [])
        + metrics.get("conservation_flags", [])
    )

    return {
        "score": score,
        "grade": _grade(score),
        "components": components,
        "interpretation": _interpret(score, metrics),
        "flags": all_flags,
        "welfare_flags":       metrics.get("welfare_flags", []),
        "regulatory_flags":    metrics.get("regulatory_flags", []),
        "conservation_flags":  metrics.get("conservation_flags", []),
        "habitat_details": {
            "available":    metrics["habitat_available"],
            "quality":      metrics.get("habitat_quality", "unknown"),
            "area_km2":     metrics.get("habitat_area_km2"),
        },
    }


def _interpret(score: float, m: dict) -> str:
    hq = m.get("habitat_quality", "unknown")
    role = m.get("ecological_role", "unknown").replace("_", " ")
    filled = m.get("ecological_role_currently_filled", True)
    n_welfare = len(m.get("welfare_flags", []))
    n_cons = len(m.get("conservation_flags", []))
    area = m.get("habitat_area_km2")

    habitat_str = f"{area:,} km²" if area else "area unquantified"
    niche_str = "unfilled since extinction" if not filled else "currently occupied"

    if score >= 70:
        return (
            f"Favourable ethical/ecological profile. Suitable habitat available "
            f"({hq}, {habitat_str}). The {role} niche is {niche_str}. "
            f"Welfare and regulatory concerns are minimal."
        )
    elif score >= 50:
        return (
            f"Mixed ethical/ecological profile. Positive factors: suitable habitat "
            f"({hq}, {habitat_str}), {role} niche {niche_str}. "
            f"Concerns: {n_welfare} welfare flag(s) and {n_cons} conservation "
            f"conflict(s) require mitigation before a responsible programme could proceed."
        )
    elif score >= 35:
        return (
            f"Significant ethical/ecological barriers. Habitat quality is {hq}; "
            f"{role} role is {niche_str}. {n_welfare} welfare concern(s) and "
            f"{n_cons} conservation conflict(s) require substantive resolution. "
            f"Community and regulatory engagement is prerequisite."
        )
    else:
        return (
            f"Major ethical/ecological obstacles. Habitat {'unavailable' if not m['habitat_available'] else hq}. "
            f"The {role} niche is {niche_str}. Multiple welfare and "
            f"conservation concerns ({n_welfare + n_cons} total) represent "
            f"barriers that go beyond technical feasibility."
        )


def _grade(score: float) -> str:
    if score >= 80: return "A"
    if score >= 65: return "B"
    if score >= 50: return "C"
    if score >= 35: return "D"
    return "F"
