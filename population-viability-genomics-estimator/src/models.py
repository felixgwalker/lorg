"""Data models for Population Viability Genomics Estimator."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class ViabilityOutcome(Enum):
    VIABLE = "viable"
    QUASI_VIABLE = "quasi_viable"
    VULNERABLE = "vulnerable"
    CRITICAL = "critical"


@dataclass
class GenomicViabilityMetrics:
    population_id: str
    census_size: int | None
    effective_size: float
    inbreeding_coefficient: float
    genetic_load_estimate: float | None
    adaptive_diversity_score: float
    mutation_accumulation_rate: float | None


@dataclass
class ViabilityProjection:
    time_horizon_years: int
    probability_persistence: float
    expected_ne_at_horizon: float
    expected_inbreeding_at_horizon: float
    outcome: ViabilityOutcome


@dataclass
class PipelineResult:
    metrics: GenomicViabilityMetrics | None
    projections: list[ViabilityProjection]
    minimum_viable_population_genomic: int | None
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
