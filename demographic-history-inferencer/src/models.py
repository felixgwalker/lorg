"""Data models for Demographic History Inferencer."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class DemographicModel(Enum):
    CONSTANT = "constant"
    EXPONENTIAL_GROWTH = "exponential_growth"
    TWO_EPOCH = "two_epoch"
    THREE_EPOCH = "three_epoch"
    ISOLATION_WITH_MIGRATION = "isolation_with_migration"


@dataclass
class EpochParameters:
    ne: float
    start_time_gen: float
    end_time_gen: float | None
    growth_rate: float = 0.0


@dataclass
class DemographicFit:
    model: DemographicModel
    epochs: list[EpochParameters]
    log_likelihood: float
    aic: float
    bic: float
    n_params: int
    bootstrap_ci: dict[str, tuple[float, float]] = field(default_factory=dict)


@dataclass
class PipelineResult:
    best_model: DemographicFit | None
    all_fits: list[DemographicFit]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
