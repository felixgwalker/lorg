"""Data models for DNA Fragmentation Profiler."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class FragmentationPattern(Enum):
    ANCIENT = "ancient"
    MODERN_LIKE = "modern_like"
    DEGRADED = "degraded"
    UNKNOWN = "unknown"


@dataclass
class FragmentLengthDistribution:
    mean_length: float
    median_length: float
    mode_length: float
    sd_length: float
    fraction_under_100bp: float
    histogram: list[int] = field(default_factory=list)
    bin_edges: list[int] = field(default_factory=list)


@dataclass
class DeaminationProfile:
    ct_freq_5prime: list[float]
    ct_freq_3prime: list[float]
    ga_freq_5prime: list[float]
    ga_freq_3prime: list[float]
    ct_rate_5prime_first_base: float
    ct_rate_3prime_first_base: float


@dataclass
class PipelineResult:
    fragment_distribution: FragmentLengthDistribution | None
    deamination_profile: DeaminationProfile | None
    fragmentation_pattern: FragmentationPattern
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
