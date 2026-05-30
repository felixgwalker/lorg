"""Configuration for pegRNA Optimiser."""

PBS_RANGE = (8, 15)
RT_RANGE = (10, 16)

SCORING_WEIGHTS: dict[str, float] = {
    "pbs_gc": 0.25,
    "rt_gc": 0.20,
    "rt_mfe": 0.25,
    "spacer_score": 0.20,
    "synthesis_penalty": 0.10,
}

PARETO_OBJECTIVES = ["efficiency_score", "synthesis_complexity"]
