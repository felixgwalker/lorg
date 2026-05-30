"""Configuration for Metabolic Pathway Balancer."""

DEFAULT_OBJECTIVE = "maximise_yield"
FBA_SOLVER = "scipy"

TOXIC_THRESHOLD_MM = 10.0
NADH_NADPH_BALANCE_TARGET = 1.0

MIN_FLUX_RATIO_FOR_BOTTLENECK = 0.1
