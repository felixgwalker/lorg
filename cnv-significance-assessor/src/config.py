"""
Configuration constants for the CNV Significance Assessor.

All thresholds, defaults, and classification parameters live here.
"""

# ── Input parsing defaults ─────────────────────────────────────────────────
DEFAULT_MIN_CNV_SIZE: int = 1_000             # bp; smaller CNVs are excluded
DEFAULT_OVERLAP_FRACTION: float = 0.1         # min fraction of gene/CNV overlap to count
DEFAULT_POP_FREQ_CUTOFF: float = 0.01         # above this frequency → likely benign override

# ── Size scoring breakpoints ──────────────────────────────────────────────
# List of (upper_limit_bp, score); first entry whose upper limit >= size wins.
# Anything above the last entry gets SIZE_SCORE_MAX.
SIZE_SCORE_BREAKS: list[tuple[int, int]] = [
    (50_000,     0),   # < 50 kb       → 0
    (500_000,    1),   # 50–500 kb     → 1
    (5_000_000,  2),   # 500 kb–5 Mb   → 2
]
SIZE_SCORE_MAX: int = 3                        # > 5 Mb

# ── Gene count scoring breakpoints ───────────────────────────────────────
# Same convention: (upper_limit_inclusive, score).
GENE_COUNT_SCORE_BREAKS: list[tuple[int, int]] = [
    (0,  0),   # 0 genes      → 0
    (2,  1),   # 1–2 genes    → 1
    (9,  2),   # 3–9 genes    → 2
]
GENE_COUNT_SCORE_MAX: int = 3                  # ≥ 10 genes

# ── Dosage sensitivity scoring breakpoints ────────────────────────────────
# Applied to maximum pHaplo/pLI (deletions) or pTriplo (duplications) across
# overlapping genes.  Score assigned to the interval [break, next_break).
DOSAGE_SCORE_BREAKS: list[tuple[float, int]] = [
    (0.0,  0),   # 0.0–0.3  → 0
    (0.3,  1),   # 0.3–0.7  → 1
    (0.7,  2),   # 0.7–high → 2
]
DOSAGE_SCORE_MAX: int = 3                      # ≥ 0.9 (high-confidence sensitive gene)

# Threshold for "high dosage sensitivity" — any gene above this in a rare CNV
# triggers a minimum VUS classification regardless of total score.
HIGH_DS_THRESHOLD: float = 0.9

# ── Population frequency modifier ────────────────────────────────────────
# Subtracted from the total score.  More common → larger subtraction.
POP_FREQ_MOD_BREAKS: list[tuple[float, int]] = [
    (0.0001, 0),   # < 0.01%   → subtract 0
    (0.001,  1),   # 0.01–0.1% → subtract 1
    (0.01,   2),   # 0.1–1%    → subtract 2
]
POP_FREQ_MOD_MAX: int = 3                      # > 1% → subtract 3

# ── Classification tier cutoffs ───────────────────────────────────────────
SCORE_BENIGN_MAX: int = 1           # total_score ≤ this  → LIKELY_BENIGN
SCORE_PATHOGENIC_MIN: int = 5       # total_score ≥ this  → LIKELY_PATHOGENIC
                                    # between              → VUS

# ── GFF3 feature sets ─────────────────────────────────────────────────────
GFF3_GENE_FEATURES: set[str] = {
    "gene", "pseudogene", "lncRNA_gene", "miRNA_gene", "ncRNA_gene", "snoRNA_gene",
}
GFF3_REGULATORY_FEATURES: set[str] = {
    "enhancer", "promoter", "regulatory_region", "CTCF_binding_site",
    "DNase_I_hypersensitive_site", "TF_binding_site", "insulator",
    "open_chromatin_region",
}

# ── Output filenames ──────────────────────────────────────────────────────
OUTFILE_ANNOTATED_CSV: str          = "cnv_annotated.csv"
OUTFILE_SIGNIFICANCE_TXT: str       = "significance_summary.txt"
OUTFILE_SIGNIFICANCE_JSON: str      = "significance_summary.json"
OUTFILE_GENE_IMPACT_CSV: str        = "gene_impact_report.csv"
OUTFILE_PLOT_PNG: str               = "cnv_ideogram.png"
OUTFILE_PLOT_SVG: str               = "cnv_ideogram.svg"
