"""
Scoring weights and thresholds for the Genomic Resurrection Scorer.

All weights within a layer must sum to 1.0.
Layer weights in LAYER_WEIGHTS must sum to 1.0.
"""

# Contribution of each layer to the overall feasibility index
LAYER_WEIGHTS: dict[str, float] = {
    "ancient_dna_quality":   0.20,
    "genomic_completeness":  0.25,
    "divergence":            0.20,
    "edit_burden":           0.25,
    "ethical_ecological":    0.10,
}

# Grade thresholds: (minimum_score, grade_letter, grade_label)
GRADE_THRESHOLDS: list[tuple[float, str, str]] = [
    (80, "A", "Highly Feasible"),
    (65, "B", "Feasible with Significant Effort"),
    (50, "C", "Marginally Feasible"),
    (35, "D", "Extremely Challenging"),
    ( 0, "F", "Currently Infeasible"),
]

# ── Layer 1: Ancient DNA Quality ──────────────────────────────────────────────
ADQ_WEIGHTS: dict[str, float] = {
    "coverage_depth":   0.30,
    "coverage_breadth": 0.25,
    "fragment_length":  0.20,
    "contamination":    0.25,
}

# ── Layer 2: Genomic Completeness ────────────────────────────────────────────
GC_WEIGHTS: dict[str, float] = {
    "overall_breadth":     0.40,
    "coding_coverage":     0.35,
    "regulatory_coverage": 0.25,
}

# ── Layer 3: Divergence ───────────────────────────────────────────────────────
DIV_WEIGHTS: dict[str, float] = {
    "coding_divergence":      0.40,
    "regulatory_divergence":  0.35,
    "genome_wide_divergence": 0.25,
}

# ── Layer 4: Edit Burden ──────────────────────────────────────────────────────
EB_WEIGHTS: dict[str, float] = {
    "total_edit_count":     0.35,
    "coding_edit_burden":   0.30,
    "regulatory_edit_burden": 0.20,
    "structural_complexity": 0.15,
}

# ── Layer 5: Ethical / Ecological ─────────────────────────────────────────────
EE_WEIGHTS: dict[str, float] = {
    "habitat":                  0.30,
    "ecological_role":          0.25,
    "welfare":                  0.25,
    "regulatory_conservation":  0.20,
}
