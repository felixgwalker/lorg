"""Configuration for Demographic History Inferencer."""

GENERATION_TIME = 30.0
MUTATION_RATE = 1.25e-8

DEFAULT_MODELS = ["constant", "exponential_growth", "two_epoch", "three_epoch"]

OPTIMISATION_REPS = 5
BOOTSTRAP_REPS = 100

NE_LOWER_BOUND = 100
NE_UPPER_BOUND = 1e7

TIME_LOWER_BOUND = 10
TIME_UPPER_BOUND = 1e6
