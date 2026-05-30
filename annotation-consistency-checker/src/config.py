"""Configuration for Annotation Consistency Checker."""

REQUIRED_FIELDS = ["gene_id", "transcript_id"]
ALLOWED_FEATURE_TYPES = ["gene", "transcript", "exon", "CDS", "start_codon", "stop_codon", "UTR"]

MAX_ACCEPTABLE_ISSUES_FRACTION = 0.01
OVERLAP_TOLERANCE_BP = 0
