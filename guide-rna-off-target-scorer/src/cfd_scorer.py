"""CFD scoring for CRISPR off-target sites."""

KMER_LEN = 20

# Position weights: linear ramp, positions 1-12 from PAM-proximal end have highest weight.
# Index 0 = PAM-distal end, index 19 = PAM-proximal end.
def _build_position_weights() -> list[float]:
    weights = []
    for i in range(KMER_LEN):
        pam_proximal_rank = i  # 0=most distal, 19=most proximal
        if pam_proximal_rank < 8:
            w = 0.1 + (pam_proximal_rank / 7) * 0.3  # 0.1 to 0.4 for distal
        else:
            w = 0.4 + ((pam_proximal_rank - 8) / 11) * 0.6  # 0.4 to 1.0 for proximal
        weights.append(round(w, 4))
    return weights


POSITION_WEIGHTS = _build_position_weights()
MISMATCH_PENALTY = 0.6


def cfd_score(guide: str, target: str) -> float:
    """Compute CFD score between guide and target (both 20nt)."""
    guide = guide.upper()
    target = target.upper()
    score = 1.0
    for i in range(KMER_LEN):
        if guide[i] != target[i]:
            w = POSITION_WEIGHTS[i]
            score *= (1.0 - MISMATCH_PENALTY * w)
    return round(score, 6)


def score_hits(hits: list[dict]) -> list[dict]:
    """Add cfd_score to each hit dict in-place and return."""
    for hit in hits:
        mm = hit["mismatches"]
        if mm == 0:
            hit["cfd_score"] = 1.0
        else:
            hit["cfd_score"] = cfd_score(hit["guide_seq"], hit["target_seq"])
        # Penalize missing PAM
        if not hit.get("has_pam", False) and mm > 0:
            hit["cfd_score"] *= 0.1
    return hits
