"""Configuration for Constraint Region Detector."""

LOEUF_THRESHOLD = 0.35
PLI_THRESHOLD = 0.9
Z_SCORE_THRESHOLD = 3.09
OE_RATIO_THRESHOLD = 0.35

HIGHLY_CONSTRAINED_LOEUF = 0.2
HIGHLY_CONSTRAINED_PLI = 0.99

GNOMAD_CONSTRAINT_COLUMNS = [
    "gene", "transcript", "loeuf", "pLI", "mis_z", "oe_lof", "oe_lof_upper",
]

DEFAULT_CONSTRAINT_METRIC = "LOEUF"
