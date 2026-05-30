"""Data models for Assembly Gap Analyser."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class GapType(Enum):
    CONTIG = "contig"
    SCAFFOLD = "scaffold"
    CENTROMERE = "centromere"
    TELOMERE = "telomere"
    UNKNOWN = "unknown"


@dataclass
class AssemblyGap:
    chrom: str
    start: int
    end: int
    length_bp: int
    gap_type: GapType
    flanking_gene_a: str | None
    flanking_gene_b: str | None
    missing_genes: list[str] = field(default_factory=list)


@dataclass
class GapSummary:
    n_gaps: int
    total_gap_length_bp: int
    mean_gap_length_bp: float
    gaps_in_gene_bodies: int
    gaps_between_syntenous_genes: int


@dataclass
class PipelineResult:
    gaps: list[AssemblyGap]
    summary: GapSummary | None
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
