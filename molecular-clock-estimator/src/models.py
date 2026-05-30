"""Data models for Molecular Clock Estimator."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class ClockModel(Enum):
    STRICT = "strict"
    RELAXED_LOGNORMAL = "relaxed_lognormal"
    RELAXED_EXPONENTIAL = "relaxed_exponential"


@dataclass
class CalibrationPoint:
    node_id: str
    min_age_mya: float
    max_age_mya: float
    prior: str
    source: str


@dataclass
class ClockEstimate:
    model: ClockModel
    substitution_rate: float
    rate_ci_lower: float
    rate_ci_upper: float
    coefficient_of_variation: float | None
    log_marginal_likelihood: float | None


@dataclass
class PipelineResult:
    estimates: list[ClockEstimate]
    best_model: ClockModel | None
    calibration_points: list[CalibrationPoint]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
