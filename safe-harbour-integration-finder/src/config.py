"""Configuration for Safe Harbour Integration Finder."""

KNOWN_SAFE_HARBOURS: dict[str, dict] = {
    "AAVS1": {"species": "human", "chromosome": "chr19", "position": 55_115_263,
               "gene": "PPP1R12C", "description": "Most widely used human safe harbour"},
    "H11":   {"species": "human", "chromosome": "chr22", "position": 39_000_000,
               "gene": "H11", "description": "H11 locus (chr22)"},
    "Rosa26":{"species": "mouse", "chromosome": "chr6",  "position": 113_073_564,
               "gene": "ROSA26", "description": "Mouse Rosa26 safe harbour"},
    "CCR5":  {"species": "human", "chromosome": "chr3",  "position": 46_414_474,
               "gene": "CCR5", "description": "CCR5 disruption as safe harbour"},
}

ONCOGENE_DISTANCE_MIN = 1_000_000
REGULATORY_DISTANCE_MIN = 50_000
REPEAT_DENSITY_MAX = 0.30
