"""Data models for Genetic Rescue Candidate Selector."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class CompatibilityRisk(Enum):
    LOW = "low"
    MODERATE = "moderate"
    HIGH = "high"


@dataclass
class DonorPopulation:
    population_id: str
    n_samples: int
    mean_heterozygosity: float
    kinship_to_recipient: float
    geographic_distance_km: float | None
    ecotype_match: bool


@dataclass
class RescueCandidate:
    donor: DonorPopulation
    rescue_score: float
    expected_heterozygosity_gain: float
    inbreeding_depression_relief: float
    compatibility_risk: CompatibilityRisk
    outbreeding_depression_risk: float
    recommendation: str


@dataclass
class PipelineResult:
    candidates: list[RescueCandidate]
    recipient_population_id: str
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
