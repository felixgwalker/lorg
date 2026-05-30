"""Data models for Palaeogenomic Coverage Assessor."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class CoverageClass(Enum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"
    INSUFFICIENT = "insufficient"


@dataclass
class ChromosomeCoverage:
    chrom: str
    length: int
    mean_depth: float
    fraction_covered: float
    breadth_1x: float
    breadth_5x: float
    breadth_10x: float


@dataclass
class CoverageAssessment:
    total_reads: int
    mapped_reads: int
    mapping_rate: float
    mean_depth: float
    median_depth: float
    fraction_covered_1x: float
    fraction_covered_5x: float
    duplication_rate: float
    endogenous_fraction: float
    coverage_class: CoverageClass
    per_chromosome: list[ChromosomeCoverage] = field(default_factory=list)


@dataclass
class PipelineResult:
    assessment: CoverageAssessment | None
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
