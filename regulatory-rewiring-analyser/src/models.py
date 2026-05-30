"""Data models for Regulatory Rewiring Analyser."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class RewiringType(Enum):
    GAINED = "gained"
    LOST = "lost"
    CONSERVED = "conserved"
    RELOCATED = "relocated"


@dataclass
class RegulatoryElement:
    element_id: str
    species: str
    chrom: str
    start: int
    end: int
    element_type: str
    target_gene: str | None


@dataclass
class RegulatoryRewiring:
    element_a: RegulatoryElement
    element_b: RegulatoryElement | None
    ortholog_gene: str
    rewiring_type: RewiringType
    sequence_conservation: float | None
    position_conservation: bool
    tfbs_delta: float | None


@dataclass
class PipelineResult:
    rewirings: list[RegulatoryRewiring]
    n_elements_assessed: int
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
