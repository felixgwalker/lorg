"""Configuration for Structural Variant Prioritiser."""

SV_TYPES = ["DEL", "DUP", "INV", "BND", "INS", "CNV"]

MIN_SV_SIZE = 50
MAX_AF = 0.01

HAPLOINSUFFICIENCY_THRESHOLD = 0.9
TRIPLOSENSITIVITY_THRESHOLD = 0.9

DECIPHER_MATCH_OVERLAP = 0.5
CLINVAR_SV_MATCH_OVERLAP = 0.5

TIER_1_SCORE_THRESHOLD = 0.75
TIER_2_SCORE_THRESHOLD = 0.50

SCORE_WEIGHTS: dict[str, float] = {
    "rare": 0.20,
    "gene_overlap": 0.25,
    "hi_score": 0.20,
    "exon_disruption": 0.15,
    "decipher_match": 0.10,
    "clinvar_match": 0.10,
}
