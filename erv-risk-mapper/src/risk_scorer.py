"""Activity risk scoring for ERV elements (0-8 scale, tiered)."""


def score_orf(longest_orf_bp: int) -> int:
    """ORF completeness -> 0-3 points. Threshold: >500 aa equivalent = 1500 bp."""
    if longest_orf_bp >= 1500:
        return 3
    elif longest_orf_bp >= 750:
        return 2
    elif longest_orf_bp >= 300:
        return 1
    return 0


def score_ltr(ltr_type: str) -> int:
    """LTR integrity -> 0-3 points."""
    if ltr_type == "full":
        return 3
    elif ltr_type == "partial":
        return 1
    return 0


def score_gc(gc_content: float) -> int:
    """GC similarity to active elements (0.45-0.55 optimal) -> 0-2 points."""
    if 0.45 <= gc_content <= 0.55:
        return 2
    elif 0.40 <= gc_content <= 0.60:
        return 1
    return 0


def risk_tier(total_score: int) -> str:
    if total_score > 5:
        return "high"
    elif total_score >= 3:
        return "moderate"
    return "low"


def score_hit(hit: dict) -> dict:
    """Add risk scoring fields to a hit dict and return it."""
    orf_pts = score_orf(hit.get("longest_orf_bp", 0))
    ltr_pts = score_ltr(hit.get("ltr_type", "solo"))
    gc_pts = score_gc(hit.get("gc_content", 0.5))
    total = orf_pts + ltr_pts + gc_pts
    hit["score_orf"] = orf_pts
    hit["score_ltr"] = ltr_pts
    hit["score_gc"] = gc_pts
    hit["risk_score"] = total
    hit["risk_tier"] = risk_tier(total)
    return hit


def score_all(hits: list[dict]) -> list[dict]:
    return [score_hit(h) for h in hits]
