"""Data models for Cas Variant Selector."""

from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class CasVariant:
    name: str
    pam: str
    size_aa: int
    editing_types: list[str]
    aav_compatible: bool


@dataclass
class CasVariantRanking:
    variant: CasVariant
    pam_site_count: int
    goal_compatibility: float
    delivery_score: float
    composite_score: float
    notes: list[str] = field(default_factory=list)


@dataclass
class PipelineResult:
    rankings: list[CasVariantRanking]
    editing_goal: str
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
