"""Configuration for CRISPR Knockout Designer."""

EARLY_EXON_FRACTION = 0.33
PAM_MOTIF_DEFAULT = "NGG"
SPACER_LENGTH = 20

DOENCH_WEIGHTS: dict[str, float] = {
    "gc_content": 0.20,
    "seed_gc": 0.15,
    "a20_penalty": 0.10,
    "g_clamp": 0.10,
    "poly_t": 0.10,
    "thermodynamic": 0.35,
}

FRAMESHIFT_PREFERRED_EXON_LIMIT = 3
N_TOP_GUIDES_DEFAULT = 5
