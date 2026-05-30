"""Data models for Genomic Diversity Index."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class DiversityMetric(Enum):
    THETA_W = "theta_w"
    THETA_PI = "theta_pi"
    TAJIMAS_D = "tajimas_d"
    HO = "Ho"
    HE = "He"
    FIS = "Fis"


@dataclass
class WindowDiversity:
    chrom: str
    start: int
    end: int
    n_snps: int
    theta_w: float
    theta_pi: float
    tajimas_d: float
    ho: float | None
    he: float | None
    fis: float | None


@dataclass
class PopulationDiversityIndex:
    population: str
    n_samples: int
    genome_wide_theta_w: float
    genome_wide_theta_pi: float
    mean_tajimas_d: float
    mean_ho: float | None
    mean_he: float | None
    mean_fis: float | None
    n_windows: int


@dataclass
class PipelineResult:
    population_index: PopulationDiversityIndex | None
    windows: list[WindowDiversity]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
