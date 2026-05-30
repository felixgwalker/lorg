"""Data models for Adaptive Diversity Scorer."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class AdaptiveDiversityClass(Enum):
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    CRITICALLY_LOW = "critically_low"


@dataclass
class AdaptiveLocus:
    locus_id: str
    chrom: str
    pos: int
    function: str | None
    is_putatively_adaptive: bool
    fst_outlier: bool
    environmental_correlation: float | None


@dataclass
class AdaptiveDiversityScore:
    population_id: str
    n_adaptive_loci: int
    adaptive_heterozygosity: float
    neutral_heterozygosity: float
    adaptive_to_neutral_ratio: float
    diversity_class: AdaptiveDiversityClass
    climate_adaptive_score: float | None


@dataclass
class PipelineResult:
    score: AdaptiveDiversityScore | None
    adaptive_loci: list[AdaptiveLocus]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
