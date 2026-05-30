"""Data models for Guide RNA GC Optimiser."""

from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class GuideGCScore:
    guide_id: str
    spacer: str
    total_gc: float
    seed_gc: float
    homopolymer_penalty: float
    poly_t_penalty: float
    composite_score: float
    passed_gc_filter: bool


@dataclass
class PipelineResult:
    scored_guides: list[GuideGCScore]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
