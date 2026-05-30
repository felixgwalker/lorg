"""Configuration for Guide RNA GC Optimiser."""

GC_OPTIMAL_MIN = 0.30
GC_OPTIMAL_MAX = 0.70
SEED_GC_OPTIMAL_MIN = 0.40
SEED_GC_OPTIMAL_MAX = 0.60
SEED_REGION_LENGTH = 12

HOMOPOLYMER_RUN_MAX = 3
POLY_T_RUN_MAX = 3

FEATURE_WEIGHTS: dict[str, float] = {
    "total_gc": 0.35,
    "seed_gc": 0.35,
    "homopolymer_penalty": 0.15,
    "poly_t_penalty": 0.15,
}
