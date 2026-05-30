"""Data models for Genome Edit Feasibility Scorer."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class FeasibilityBand(Enum):
    UNFEASIBLE = "unfeasible"
    CHALLENGING = "challenging"
    FEASIBLE = "feasible"
    HIGHLY_FEASIBLE = "highly_feasible"


@dataclass
class FeasibilityComponent:
    name: str
    score: float
    weight: float
    notes: str = ""


@dataclass
class FeasibilityScore:
    composite_score: float
    band: FeasibilityBand
    components: list[FeasibilityComponent]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
