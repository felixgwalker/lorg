"""Data models for CRISPR Knockout Designer."""

from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class KnockoutGuide:
    guide_id: str
    spacer: str
    pam: str
    position: int
    strand: str
    on_target_score: float
    predicted_frameshift_pct: float
    feature_scores: dict[str, float] = field(default_factory=dict)


@dataclass
class KnockoutDesign:
    gene_id: str
    guides: list[KnockoutGuide]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
