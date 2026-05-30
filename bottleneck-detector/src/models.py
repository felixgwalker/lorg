"""Data models for Bottleneck Detector."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class BottleneckTest(Enum):
    TAJIMAS_D = "tajimas_d"
    HEW = "HEW"
    MODE_SHIFT = "mode_shift"
    M_RATIO = "m_ratio"


class BottleneckSignal(Enum):
    STRONG = "strong"
    MODERATE = "moderate"
    WEAK = "weak"
    NONE = "none"


@dataclass
class SFSShape:
    n_singletons: int
    n_doubletons: int
    total_variants: int
    l_shape_ratio: float
    mode_shift_detected: bool


@dataclass
class BottleneckResult:
    test: BottleneckTest
    statistic: float
    p_value: float | None
    signal: BottleneckSignal
    tajimas_d: float | None
    m_ratio: float | None


@dataclass
class PipelineResult:
    sfs_shape: SFSShape | None
    results: list[BottleneckResult]
    combined_signal: BottleneckSignal
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
