"""Data models for Assembly Quality Assessor."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class AssemblyQualityClass(Enum):
    REFERENCE_QUALITY = "reference_quality"
    CHROMOSOME_LEVEL = "chromosome_level"
    SCAFFOLD_LEVEL = "scaffold_level"
    CONTIG_LEVEL = "contig_level"


@dataclass
class AssemblyStats:
    n_sequences: int
    total_length_bp: int
    n50: int
    n90: int
    l50: int
    l90: int
    largest_sequence_bp: int
    n_gaps: int
    gap_length_bp: int
    gc_content: float
    n_ambiguous_bases: int
    quality_class: AssemblyQualityClass


@dataclass
class PipelineResult:
    stats: AssemblyStats | None
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
