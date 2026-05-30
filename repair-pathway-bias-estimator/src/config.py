"""Configuration for Repair Pathway Bias Estimator."""

CELL_TYPE_HDR_BIAS: dict[str, float] = {
    "HEK293":     0.25,
    "HeLa":       0.20,
    "iPSC":       0.30,
    "primary_T":  0.10,
    "neuron":     0.05,
    "hepatocyte": 0.08,
    "default":    0.15,
}

MH_MIN_LENGTH = 2
MH_MAX_LENGTH = 20
MH_SEARCH_WINDOW = 30
MH_SCORE_GC_BONUS = 0.5

NHEJ_BASE_PROBABILITY = 0.60
