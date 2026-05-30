"""Data models for Enhancer Target Linker."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class LinkMethod(Enum):
    ACTIVITY_BY_CONTACT = "activity_by_contact"
    CORRELATION = "correlation"
    DISTANCE = "distance"
    HI_C = "hi_c"


@dataclass
class Enhancer:
    chrom: str
    start: int
    end: int
    element_id: str
    activity_score: float
    cell_type: str | None


@dataclass
class EnhancerTargetLink:
    enhancer: Enhancer
    target_gene: str
    target_tss_chrom: str
    target_tss_pos: int
    distance_bp: int
    link_score: float
    method: LinkMethod
    correlation: float | None
    contact_frequency: float | None


@dataclass
class PipelineResult:
    links: list[EnhancerTargetLink]
    n_enhancers: int
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
