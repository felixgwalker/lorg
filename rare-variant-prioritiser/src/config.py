"""Configuration for Rare Variant Prioritiser."""

MAX_AF = 0.01
ULTRA_RARE_AF = 0.0001

MIN_CADD_SCORE = 20.0
HIGH_CADD_SCORE = 30.0

TIER_1_SCORE_THRESHOLD = 0.75
TIER_2_SCORE_THRESHOLD = 0.50

EVIDENCE_WEIGHTS: dict[str, float] = {
    "ultra_rare": 0.25,
    "rare": 0.10,
    "conserved": 0.15,
    "high_cadd": 0.20,
    "in_panel": 0.10,
    "hpo_match": 0.20,
    "constrained_gene": 0.10,
    "clinvar_pathogenic": 0.30,
}

HIGH_IMPACT_CONSEQUENCES = [
    "stop_gained", "frameshift_variant", "splice_donor_variant",
    "splice_acceptor_variant", "start_lost", "stop_lost",
]
