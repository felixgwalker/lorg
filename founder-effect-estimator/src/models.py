"""Data models for Founder Effect Estimator."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class FounderSignature(Enum):
    STRONG = "strong"
    MODERATE = "moderate"
    WEAK = "weak"
    ABSENT = "absent"


@dataclass
class DiversityStats:
    population: str
    n_samples: int
    theta_w: float
    theta_pi: float
    tajimas_d: float
    heterozygosity: float
    n_private_variants: int


@dataclass
class FounderEffectEstimate:
    study_population: str
    reference_population: str
    diversity_ratio: float
    private_variant_fraction: float
    haplotype_block_length_kb: float
    signature: FounderSignature
    estimated_founder_size: int | None
    estimated_generations_since_founding: int | None


@dataclass
class PipelineResult:
    study_stats: DiversityStats | None
    reference_stats: DiversityStats | None
    estimate: FounderEffectEstimate | None
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
