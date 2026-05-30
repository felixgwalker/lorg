"""Data models for Promoter Strength Estimator."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class PromoterStrengthClass(Enum):
    STRONG = "strong"
    MODERATE = "moderate"
    WEAK = "weak"
    SILENT = "silent"


@dataclass
class PromoterFeatures:
    tata_box_score: float | None
    inr_score: float | None
    cpg_island: bool
    gc_content: float
    tfbs_count: int
    h3k4me3_signal: float | None
    pol2_signal: float | None


@dataclass
class PromoterStrength:
    gene_id: str
    chrom: str
    tss_pos: int
    strand: str
    features: PromoterFeatures
    composite_strength_score: float
    strength_class: PromoterStrengthClass
    predicted_expression_level: str | None


@dataclass
class PipelineResult:
    promoter_strengths: list[PromoterStrength]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
