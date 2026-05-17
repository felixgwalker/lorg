"""
Global configuration constants for HGT Risk Assessor.

All scientifically significant thresholds and weight profiles are defined here
so they can be reviewed and updated in one place without hunting through the
codebase.  Values are intentionally plain Python dicts/lists so they can be
overridden from a JSON/TOML config file in a future version.
"""

# ---------------------------------------------------------------------------
# Score bands for the three-layer HGT Risk Index
# ---------------------------------------------------------------------------
# Each tuple: (upper_bound_exclusive, label_string)
# Anything >= the last threshold receives SCORE_BAND_TOP.
# These are the bands referred to in the biosafety narrative.

SCORE_BANDS: list[tuple[float, str]] = [
    (0.25, "low"),
    (0.50, "moderate"),
    (0.75, "high"),
]
SCORE_BAND_TOP = "very_high"

# ---------------------------------------------------------------------------
# Flat aggregation thresholds (used by the original aggregator.py)
# Kept here so they are configurable alongside the three-layer bands.
# These map to RiskLevel enum values: Low, Medium, High, Critical.
# ---------------------------------------------------------------------------
FLAT_THRESHOLDS: list[tuple[float, str]] = [
    (0.25, "Low"),
    (0.50, "Medium"),
    (0.75, "High"),
]
FLAT_TOP = "Critical"

# Human-readable descriptions used in reports and explanations
SCORE_BAND_DESCRIPTIONS: dict[str, str] = {
    "low":       "Minimal sequence-level HGT indicators. Standard biosafety procedures apply.",
    "moderate":  "Some indicators present. Expert review recommended before scale-up or release.",
    "high":      "Multiple significant indicators. Formal contained use risk assessment required.",
    "very_high": "Strong HGT risk signals. Do not proceed without biosafety officer review.",
}

# ---------------------------------------------------------------------------
# Three-layer weight profiles
# ---------------------------------------------------------------------------
# Keys must be: transfer_opportunity, establishment, consequence
# Values must sum to 1.0 within each profile.
# Add further named profiles here; the CLI exposes them via --weight-profile.

WEIGHT_PROFILES: dict[str, dict[str, float]] = {
    "default": {
        "transfer_opportunity": 0.40,
        "establishment":        0.35,
        "consequence":          0.25,
    },
    "environmental": {
        # Environment release: transfer and persistence in open ecosystem matter most.
        "transfer_opportunity": 0.45,
        "establishment":        0.40,
        "consequence":          0.15,
    },
    "clinical_amr": {
        # Clinical AMR: functional payload (resistance genes) is the primary concern.
        "transfer_opportunity": 0.30,
        "establishment":        0.25,
        "consequence":          0.45,
    },
}

DEFAULT_WEIGHT_PROFILE = "default"

# ---------------------------------------------------------------------------
# Within-layer feature weights
# ---------------------------------------------------------------------------
# Feature weights within each layer must sum to 1.0.
# Features listed here are the canonical identifiers; signal modules and
# feature extractors must use these exact names.

LAYER_FEATURE_WEIGHTS: dict[str, dict[str, float]] = {
    "transfer_opportunity": {
        # IS element and integron signals are from the existing BLAST pipeline.
        # Conjugative element signal (BLASTX) is the strongest single transfer indicator.
        # Plasmid context and transposase proximity are placeholders for future integration.
        "is_element_match":      0.25,
        "integron_association":  0.20,
        "conjugative_element":   0.25,
        "plasmid_context":       0.15,
        "transposase_proximity": 0.10,
        "repeat_density":        0.05,
    },
    "establishment": {
        # GC deviation is reused from the existing gc_content signal.
        # Codon usage distance is computable from sequence alone.
        # Taxonomic distance requires optional donor/recipient inputs.
        "gc_deviation":          0.35,
        "codon_usage_distance":  0.30,
        "taxonomic_distance":    0.20,
        "promoter_plausibility": 0.10,
        "sequence_complexity":   0.05,
    },
    "consequence": {
        # AMR and virulence are placeholders pending CARD/VFDB integration.
        # Prophage context reused from existing signal.
        # Gene completeness and payload count are computable from sequence.
        "prophage_context":      0.15,
        "amr_content":           0.35,
        "virulence_flags":       0.25,
        "gene_completeness":     0.15,
        "payload_count":         0.10,
    },
}
