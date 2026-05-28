# Weights for the composite tolerance score (must sum to 1.0)
SCORE_WEIGHTS = {
    "chromatin_score": 0.30,
    "sequence_complexity_score": 0.25,
    "gene_density_score": 0.30,
    "size_penalty": 0.15,
}

# Insert size threshold above which the log-linear penalty kicks in (bp)
SIZE_PENALTY_THRESHOLD_BP = 5000

# Tolerance tier boundaries (composite score 0–1)
TIER_HIGH = 0.7
TIER_MODERATE = 0.4

# GC content range considered "open / permissive" chromatin
GC_OPEN_MIN = 0.35
GC_OPEN_MAX = 0.55
