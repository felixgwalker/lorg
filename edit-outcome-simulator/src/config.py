"""Configuration for Edit Outcome Simulator."""

MH_SCORE_WEIGHTS: dict[str, float] = {
    "length_sq": 1.0,
    "gc_bonus":  0.5,
}

MAX_DELETION_SIZE = 30
INSERTION_TEMPLATE_POSITION = 1

CAS_CUT_OFFSETS: dict[str, int] = {
    "SpCas9":   3,
    "SaCas9":   3,
    "Cas12a":   1,
    "Cas12b":   1,
}

N_SIMULATIONS_DEFAULT = 10_000
FRAMESHIFT_WINDOW = 100
