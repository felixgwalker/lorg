"""Configuration for Guide RNA Specificity Ranker."""

CFD_MISMATCH_PENALTIES: dict[str, dict[int, float]] = {
    "rA:dC": {1: 0.0, 2: 0.0, 3: 0.057, 4: 0.0, 5: 0.0, 6: 0.0, 7: 0.0,
              8: 0.0, 9: 0.0, 10: 0.0, 11: 0.0, 12: 0.0, 13: 0.0, 14: 0.0,
              15: 0.0, 16: 0.0, 17: 0.0, 18: 0.0, 19: 0.0, 20: 0.0},
}

SPECIFICITY_BANDS: list[tuple[float, str]] = [
    (0.50, "poor"),
    (0.75, "moderate"),
    (0.90, "good"),
]
SPECIFICITY_BAND_TOP = "excellent"

MAX_MISMATCHES = 3
SEED_REGION_LENGTH = 12
