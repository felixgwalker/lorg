"""Configuration for Base Edit Outcome Predictor."""

EDITING_WINDOWS: dict[str, dict[str, int]] = {
    "CBE3":   {"start": 4, "end": 8},
    "BE4max": {"start": 4, "end": 8},
    "ABE8e":  {"start": 4, "end": 8},
    "ABEmax": {"start": 4, "end": 7},
    "CBE4-max": {"start": 4, "end": 8},
}

TARGET_BASES: dict[str, str] = {
    "CBE3":   "C",
    "BE4max": "C",
    "ABE8e":  "A",
    "ABEmax": "A",
    "CBE4-max": "C",
}

INDEL_FREQUENCY: dict[str, float] = {
    "CBE3":   0.06,
    "BE4max": 0.04,
    "ABE8e":  0.02,
    "ABEmax": 0.03,
}

POSITION_EFFICIENCY: dict[int, float] = {
    4: 0.5, 5: 0.85, 6: 0.90, 7: 0.80, 8: 0.55,
}
