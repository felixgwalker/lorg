"""Data models for Constraint Region Detector."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class ConstraintMetric(Enum):
    LOEUF = "LOEUF"
    PLI = "pLI"
    Z_SCORE = "z_score"
    OE_RATIO = "oe_ratio"


class ConstraintLevel(Enum):
    HIGHLY_CONSTRAINED = "highly_constrained"
    CONSTRAINED = "constrained"
    UNCONSTRAINED = "unconstrained"


@dataclass
class ConstrainedRegion:
    chrom: str
    start: int
    end: int
    gene: str
    transcript_id: str
    loeuf: float | None
    pli: float | None
    z_score: float | None
    oe_ratio: float | None
    constraint_level: ConstraintLevel


@dataclass
class VariantConstraintOverlap:
    chrom: str
    pos: int
    ref: str
    alt: str
    overlapping_region: ConstrainedRegion | None
    is_constrained: bool
    primary_metric: ConstraintMetric
    metric_value: float | None


@dataclass
class PipelineResult:
    overlaps: list[VariantConstraintOverlap]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
