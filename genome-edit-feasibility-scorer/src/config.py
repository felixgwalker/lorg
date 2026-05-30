"""Configuration for Genome Edit Feasibility Scorer."""

FEASIBILITY_WEIGHTS: dict[str, float] = {
    "pam_density":         0.20,
    "gc_content":          0.15,
    "chromatin_access":    0.25,
    "essentiality_risk":   0.20,
    "delivery_suitability":0.20,
}

FEASIBILITY_BANDS: list[tuple[float, str]] = [
    (0.30, "unfeasible"),
    (0.55, "challenging"),
    (0.75, "feasible"),
]
FEASIBILITY_BAND_TOP = "highly_feasible"

EDIT_TYPE_WEIGHTS: dict[str, dict[str, float]] = {
    "knockout":   {"pam_density": 0.30, "chromatin_access": 0.20},
    "base-edit":  {"pam_density": 0.25, "chromatin_access": 0.25},
    "prime-edit": {"pam_density": 0.20, "chromatin_access": 0.30},
}
