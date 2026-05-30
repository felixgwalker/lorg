"""Data models for Introgression Detector."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class IntrogressionTest(Enum):
    D_STATISTIC = "D_statistic"
    F4_RATIO = "f4_ratio"
    RND_MIN = "RND_min"
    DFOIL = "Dfoil"


@dataclass
class PopulationConfig:
    p1: str
    p2: str
    p3: str
    outgroup: str


@dataclass
class IntrogressionResult:
    test: IntrogressionTest
    d_statistic: float | None
    f4_ratio: float | None
    z_score: float
    p_value: float
    n_abba: int | None
    n_baba: int | None
    significant: bool
    population_config: PopulationConfig


@dataclass
class IntrogressionSegment:
    chrom: str
    start: int
    end: int
    recipient_population: str
    donor_population: str
    length_kb: float
    mean_d: float
    confidence: float


@dataclass
class PipelineResult:
    genome_wide_result: IntrogressionResult | None
    segments: list[IntrogressionSegment]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
