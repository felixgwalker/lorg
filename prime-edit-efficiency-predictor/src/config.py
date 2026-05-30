"""Configuration for Prime Edit Efficiency Predictor."""

FEATURE_WEIGHTS: dict[str, float] = {
    "pbs_gc": 0.20,
    "rt_length": 0.15,
    "rt_gc": 0.15,
    "nick_distance": 0.20,
    "spacer_tm": 0.15,
    "rt_mfe": 0.15,
}

EFFICIENCY_BANDS: list[tuple[float, str]] = [
    (0.20, "low"),
    (0.50, "moderate"),
    (0.75, "high"),
]
EFFICIENCY_BAND_TOP = "very_high"

PBS_GC_OPTIMAL_MIN = 0.40
PBS_GC_OPTIMAL_MAX = 0.60
RT_OPTIMAL_LENGTH = 12
NICK_DISTANCE_OPTIMAL_MIN = 40
NICK_DISTANCE_OPTIMAL_MAX = 70
