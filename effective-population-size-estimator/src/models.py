"""Data models for Effective Population Size Estimator."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class NeMethod(Enum):
    WATTERSON = "watterson"
    TAJIMA = "tajima"
    LD_BASED = "ld_based"
    PSMC = "psmc"


@dataclass
class SiteFrequencySpectrum:
    counts: list[int]
    n_samples: int
    n_sites: int
    folded: bool


@dataclass
class NeEstimate:
    method: NeMethod
    ne: float
    ne_lower_ci: float
    ne_upper_ci: float
    theta: float | None
    generation_time: float


@dataclass
class TemporalNeProfile:
    time_points_gen: list[float]
    ne_values: list[float]
    method: NeMethod


@dataclass
class PipelineResult:
    estimates: list[NeEstimate]
    sfs: SiteFrequencySpectrum | None
    temporal_profile: TemporalNeProfile | None
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
