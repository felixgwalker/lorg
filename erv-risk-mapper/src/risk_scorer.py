"""Activity risk scoring for ERV elements.

Each hit receives normalised component scores (0.0–1.0) and a composite risk
score (0.0–1.0) with a three-tier label (low / moderate / high).

Component scores
----------------
orf_score    : longest_orf_bp / 9000  (pol gene ~9 kb)
ltr_score    : 1.0 full | 0.5 partial | 0.1 solo
age_score    : 1 / (1 + age_mya / 10)  (younger = riskier)
family_score : HERV-K=1.0 | HERV-W=0.8 | HERV-H=0.6 | HERV-E=0.5 | HERV-L=0.3

composite_risk = (orf_score + ltr_score + age_score + family_score) / 4

risk_tier: low (<0.3) | moderate (0.3–0.6) | high (>=0.6)
"""

FAMILY_SCORES: dict[str, float] = {
    "HERV-K": 1.0,
    "HERV-W": 0.8,
    "HERV-H": 0.6,
    "HERV-E": 0.5,
    "HERV-L": 0.3,
}

LTR_SCORES: dict[str, float] = {
    "full": 1.0,
    "partial": 0.5,
    "solo": 0.1,
}

POL_GENE_BP = 9000  # normalisation denominator for ORF score


def _orf_score(longest_orf_bp: int) -> float:
    return min(1.0, longest_orf_bp / POL_GENE_BP)


def _ltr_score(ltr_type: str) -> float:
    return LTR_SCORES.get(ltr_type, 0.1)


def _age_score(age_mya: float) -> float:
    """Younger ERVs are riskier: score decreases as age increases."""
    return 1.0 / (1.0 + age_mya / 10.0)


def _family_score(family: str) -> float:
    return FAMILY_SCORES.get(family, 0.3)


def _risk_tier(composite_risk: float) -> str:
    if composite_risk >= 0.6:
        return "high"
    elif composite_risk >= 0.3:
        return "moderate"
    return "low"


def score_hit(hit: dict) -> dict:
    """Add normalised risk scoring fields to a hit dict and return it."""
    orf_sc = _orf_score(hit.get("longest_orf_bp", 0))
    ltr_sc = _ltr_score(hit.get("ltr_type", "solo"))
    age_sc = _age_score(hit.get("age_mya", 25.0))
    fam_sc = _family_score(hit.get("family", "HERV-L"))

    composite = (orf_sc + ltr_sc + age_sc + fam_sc) / 4.0

    hit["orf_score"] = round(orf_sc, 4)
    hit["ltr_score"] = round(ltr_sc, 4)
    hit["age_score"] = round(age_sc, 4)
    hit["family_score"] = round(fam_sc, 4)
    hit["composite_risk"] = round(composite, 4)
    hit["risk_tier"] = _risk_tier(composite)
    return hit


def score_all(hits: list[dict]) -> list[dict]:
    """Score every hit and return the updated list."""
    return [score_hit(h) for h in hits]
