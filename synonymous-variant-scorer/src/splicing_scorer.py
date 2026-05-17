"""Rule-based positional splice site proximity scoring."""

from __future__ import annotations

from src.vcf_parser import SynonymousVariant


DONOR_POSITIONS = {1, 2, 3, 4, 5}
ACCEPTOR_POSITIONS = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20}

DONOR_SCORE_MAP = {1: 1.0, 2: 0.9, 3: 0.8, 4: 0.6, 5: 0.4}
ACCEPTOR_SCORE_MAP = {d: max(0.0, 1.0 - (d - 1) * 0.05) for d in ACCEPTOR_POSITIONS}


def score_splicing(variant: SynonymousVariant) -> float:
    """Return splicing disruption score in [0, 1] based on proximity to splice sites."""
    donor_dist = variant.distance_to_donor
    acceptor_dist = variant.distance_to_acceptor

    donor_score = DONOR_SCORE_MAP.get(donor_dist, 0.0)
    acceptor_score = ACCEPTOR_SCORE_MAP.get(acceptor_dist, 0.0)

    exon_frac = _exon_position_score(variant.exon_pos, variant.exon_length)

    raw = max(donor_score, acceptor_score, exon_frac)
    return min(1.0, max(0.0, raw))


def _exon_position_score(exon_pos: int, exon_length: int) -> float:
    """Slight elevation for variants in the first/last 6 bp of an exon (ESE/ESS regions)."""
    if exon_length <= 0:
        return 0.0
    from_start = exon_pos
    from_end = exon_length - exon_pos
    margin = min(from_start, from_end)
    if margin <= 6:
        return 0.3 * (1 - margin / 6)
    return 0.0
