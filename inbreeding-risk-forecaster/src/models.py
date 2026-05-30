"""Data models for Inbreeding Risk Forecaster."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class InbreedingRiskLevel(Enum):
    CRITICAL = "critical"
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"


@dataclass
class InbreedingStats:
    population_id: str
    n_samples: int
    mean_froh: float
    mean_fis: float
    mean_fst_to_outgroup: float | None
    mean_roh_length_mb: float
    n_roh_per_individual: float
    proportion_genome_in_roh: float


@dataclass
class InbreedingRiskForecast:
    population_id: str
    current_f: float
    projected_f_10gen: float
    projected_f_50gen: float
    effective_population_size: float
    risk_level: InbreedingRiskLevel
    generations_to_critical: int | None
    recommendation: str


@dataclass
class PipelineResult:
    stats: InbreedingStats | None
    forecast: InbreedingRiskForecast | None
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
