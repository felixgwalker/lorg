"""Data models for Guide RNA Specificity Ranker."""

from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class OffTargetSite:
    chromosome: str
    position: int
    strand: str
    mismatches: int
    cfd_score: float
    sequence: str


@dataclass
class GuideSpecificityResult:
    guide_id: str
    spacer: str
    specificity_score: float
    specificity_band: str
    n_offtargets: int
    offtarget_sites: list[OffTargetSite] = field(default_factory=list)


@dataclass
class PipelineResult:
    ranked_guides: list[GuideSpecificityResult]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
