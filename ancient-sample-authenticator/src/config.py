"""Configuration for Ancient Sample Authenticator."""

MIN_MEAN_FRAGMENT_LENGTH = 35
MAX_MEAN_FRAGMENT_LENGTH = 200

CT_RATE_5PRIME_MIN = 0.05
CONTAMINATION_MAX = 0.05

MIN_COVERAGE = 0.5
MIN_ENDOGENOUS_FRACTION = 0.01

AUTHENTIC_SCORE_THRESHOLD = 0.75
LIKELY_AUTHENTIC_SCORE_THRESHOLD = 0.50

CRITERION_WEIGHTS: dict[str, float] = {
    "fragment_length": 0.20,
    "deamination_damage": 0.30,
    "contamination": 0.25,
    "endogenous_fraction": 0.15,
    "coverage": 0.10,
}
