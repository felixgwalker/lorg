"""Data models for Enhancer Conservation Analyser."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class ConservationLevel(Enum):
    HIGHLY_CONSERVED = "highly_conserved"
    MODERATELY_CONSERVED = "moderately_conserved"
    LINEAGE_SPECIFIC = "lineage_specific"
    UNKNOWN = "unknown"


@dataclass
class EnhancerConservation:
    element_id: str
    chrom: str
    start: int
    end: int
    mean_phastcons: float | None
    mean_phylop: float | None
    sequence_conservation: float | None
    n_species_conserved: int
    conservation_level: ConservationLevel
    functional_conservation: bool | None


@dataclass
class PipelineResult:
    conserved_enhancers: list[EnhancerConservation]
    n_elements_assessed: int
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
