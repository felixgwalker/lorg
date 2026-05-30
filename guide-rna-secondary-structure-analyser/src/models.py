"""Data models for Guide RNA Secondary Structure Analyser."""

from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class GuideStructure:
    guide_id: str
    spacer: str
    full_guide_sequence: str
    mfe_kcal_mol: float
    dot_bracket: str
    seed_accessibility: float
    scaffold_duplex_flag: bool
    passed_structure_filter: bool


@dataclass
class PipelineResult:
    structures: list[GuideStructure]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
