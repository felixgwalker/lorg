"""Configuration for Conservation Priority Ranker."""

CRITICAL_NE_THRESHOLD = 50
HIGH_NE_THRESHOLD = 100

CRITICAL_INBREEDING_THRESHOLD = 0.25
UNIQUE_ALLELE_WEIGHT = 0.25

SCORE_WEIGHTS: dict[str, float] = {
    "inbreeding": 0.30,
    "ne_size": 0.25,
    "adaptive_diversity": 0.20,
    "unique_alleles": 0.15,
    "threat_status": 0.10,
}
