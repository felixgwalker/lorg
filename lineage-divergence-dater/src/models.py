"""Data models for Lineage Divergence Dater."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class DatingMethod(Enum):
    MOLECULAR_CLOCK = "molecular_clock"
    BAYESIAN_DATED_TIPS = "bayesian_dated_tips"
    PAIRWISE_DISTANCE = "pairwise_distance"


@dataclass
class DivergenceDate:
    lineage_a: str
    lineage_b: str
    divergence_mya: float
    ci_lower_mya: float
    ci_upper_mya: float
    method: DatingMethod
    n_loci: int
    calibrated: bool


@dataclass
class PipelineResult:
    divergence_dates: list[DivergenceDate]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
