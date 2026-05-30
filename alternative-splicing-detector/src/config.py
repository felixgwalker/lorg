"""Configuration for Alternative Splicing Detector."""

MIN_READS_PER_EVENT = 10
MIN_DELTA_PSI = 0.1
FDR_THRESHOLD = 0.05

EVENT_TYPES = ["exon_skipping", "intron_retention", "alt_5_prime", "alt_3_prime"]

MIN_COVERAGE_PER_JUNCTION = 5
PSI_RANGE = (0.05, 0.95)
