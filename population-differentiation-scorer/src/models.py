"""Data models for Population Differentiation Scorer."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class DifferentiationMetric(Enum):
    FST = "Fst"
    GST = "Gst"
    JOST_D = "Jost_D"
    PHI_ST = "Phi_st"


@dataclass
class PairwiseDifferentiation:
    pop_a: str
    pop_b: str
    fst: float
    gst: float | None
    jost_d: float | None
    n_snps: int
    p_value: float | None
    ci_lower: float | None
    ci_upper: float | None


@dataclass
class WindowDifferentiation:
    chrom: str
    start: int
    end: int
    pop_a: str
    pop_b: str
    fst: float
    n_snps: int
    is_outlier: bool


@dataclass
class PipelineResult:
    genome_wide: list[PairwiseDifferentiation]
    windows: list[WindowDifferentiation]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
