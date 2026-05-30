"""Data models for Selection Sweep Detector."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class SweepTest(Enum):
    IHS = "iHS"
    XP_EHH = "XP-EHH"
    CLR = "CLR"
    TAJIMAS_D = "tajimas_d"


class SweepType(Enum):
    HARD = "hard"
    SOFT = "soft"
    INCOMPLETE = "incomplete"
    NONE = "none"


@dataclass
class SweepWindow:
    chrom: str
    start: int
    end: int
    population: str
    ihs_score: float | None
    xp_ehh_score: float | None
    clr_score: float | None
    tajimas_d: float | None
    composite_score: float
    is_outlier: bool
    sweep_type: SweepType


@dataclass
class SweepRegion:
    chrom: str
    start: int
    end: int
    length_kb: float
    population: str
    peak_score: float
    sweep_type: SweepType
    candidate_genes: list[str] = field(default_factory=list)


@dataclass
class PipelineResult:
    windows: list[SweepWindow]
    sweep_regions: list[SweepRegion]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
