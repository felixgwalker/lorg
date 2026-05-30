"""Data models for Genome Rearrangement Mapper."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class RearrangementType(Enum):
    INVERSION = "inversion"
    TRANSLOCATION = "translocation"
    FUSION = "fusion"
    FISSION = "fission"
    TRANSPOSITION = "transposition"


@dataclass
class SyntenyBreakpoint:
    chrom_a: str
    pos_a: int
    chrom_b: str
    pos_b: int
    species_a: str
    species_b: str
    breakpoint_type: RearrangementType
    flanking_gene_a: str | None
    flanking_gene_b: str | None


@dataclass
class ChromosomeRearrangement:
    rearrangement_id: str
    rearrangement_type: RearrangementType
    species_a: str
    species_b: str
    region_a: str
    region_b: str
    size_kb: float | None
    breakpoints: list[SyntenyBreakpoint] = field(default_factory=list)


@dataclass
class PipelineResult:
    rearrangements: list[ChromosomeRearrangement]
    n_breakpoints: int
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
