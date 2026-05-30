"""Data models for Multiplex Edit Planner."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class CompatibilityLevel(Enum):
    COMPATIBLE = "compatible"
    CAUTION = "caution"
    INCOMPATIBLE = "incompatible"


@dataclass
class EditTarget:
    target_id: str
    guide_spacer: str
    chromosome: str | None
    position: int | None
    edit_type: str


@dataclass
class EditCompatibility:
    target_a: str
    target_b: str
    cfd_cross_reactivity: float
    window_overlap: bool
    translocation_risk: bool
    compatibility: CompatibilityLevel


@dataclass
class MultiplexPlan:
    targets: list[EditTarget]
    compatibility_matrix: list[EditCompatibility]
    batches: list[list[str]]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
