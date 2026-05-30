"""Data models for CRISPR Array Designer."""

from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class SpacerElement:
    spacer_id: str
    target_id: str
    spacer_sequence: str
    pam: str
    strand: str
    gc_fraction: float
    uniqueness_score: float


@dataclass
class ArrayDesign:
    cas_system: str
    spacers: list[SpacerElement]
    direct_repeat: str
    array_sequence: str
    n_spacers: int


@dataclass
class PipelineResult:
    array_design: ArrayDesign | None
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
