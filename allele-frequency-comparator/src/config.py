"""Configuration for Allele Frequency Comparator."""

GNOMAD_POPULATIONS = ["afr", "amr", "asj", "eas", "fin", "nfe", "sas", "oth"]

FOLD_CHANGE_THRESHOLD = 5.0
MIN_AF_FOR_COMPARISON = 1e-5
POPULATION_SPECIFIC_AF_THRESHOLD = 0.01

FISHER_ALPHA = 0.05

FST_ESTIMATOR = "weir_cockerham"
