"""Data models for Contig Scaffolding Helper."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class EvidenceType(Enum):
    PAIRED_READS = "paired_reads"
    HI_C = "hi_c"
    REFERENCE_GUIDED = "reference_guided"
    OPTICAL_MAPPING = "optical_mapping"


@dataclass
class ContigLink:
    contig_a: str
    contig_b: str
    orientation: str
    gap_estimate_bp: int | None
    link_support: int
    evidence_type: EvidenceType


@dataclass
class Scaffold:
    scaffold_id: str
    contigs: list[str]
    orientations: list[str]
    gaps: list[int]
    total_length_bp: int
    n_contigs: int


@dataclass
class PipelineResult:
    scaffolds: list[Scaffold]
    unplaced_contigs: list[str]
    n_contigs_placed: int
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
