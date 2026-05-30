"""Configuration for Compound Heterozygosity Detector."""

MAX_AF_COMPHET = 0.01
MIN_GENOTYPE_QUALITY = 20
MIN_DEPTH = 10

HIGH_IMPACT_CONSEQUENCES = [
    "stop_gained", "frameshift_variant", "splice_donor_variant",
    "splice_acceptor_variant", "missense_variant", "start_lost",
]

TRIO_ROLES = ["proband", "father", "mother"]

PHASE_BY_TRANSMISSION = True
