"""Data models for Conservation Priority Ranker."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class PriorityTier(Enum):
    CRITICAL = "critical"
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


@dataclass
class PopulationMetrics:
    population_id: str
    effective_population_size: float
    inbreeding_coefficient: float
    adaptive_diversity_score: float
    unique_allele_fraction: float
    threat_status: str | None
    population_size_census: int | None


@dataclass
class ConservationPriorityScore:
    population: PopulationMetrics
    composite_score: float
    tier: PriorityTier
    urgency_years: int | None
    recommended_actions: list[str] = field(default_factory=list)


@dataclass
class PipelineResult:
    ranked_populations: list[ConservationPriorityScore]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
