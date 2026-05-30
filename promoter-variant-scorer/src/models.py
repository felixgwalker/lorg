"""Data models for Promoter Variant Scorer."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class TFBSEffect(Enum):
    DISRUPTION = "disruption"
    CREATION = "creation"
    NEUTRAL = "neutral"


@dataclass
class PromoterVariant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    distance_to_tss: int
    strand: str


@dataclass
class PWMHit:
    tf_name: str
    matrix_id: str
    start: int
    end: int
    strand: str
    ref_score: float
    ic_content: float


@dataclass
class TFBSDisruption:
    variant: PromoterVariant
    tf_hit: PWMHit
    alt_score: float
    delta_score: float
    effect: TFBSEffect
    p_value: float | None


@dataclass
class PipelineResult:
    disruptions: list[TFBSDisruption]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
