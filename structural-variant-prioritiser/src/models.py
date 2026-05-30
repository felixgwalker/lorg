"""Data models for Structural Variant Prioritiser."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class SVType(Enum):
    DEL = "DEL"
    DUP = "DUP"
    INV = "INV"
    BND = "BND"
    INS = "INS"
    CNV = "CNV"


class SVPriorityTier(Enum):
    TIER_1 = "tier_1"
    TIER_2 = "tier_2"
    TIER_3 = "tier_3"
    UNRANKED = "unranked"


@dataclass
class StructuralVariant:
    chrom: str
    start: int
    end: int
    sv_type: SVType
    size: int | None
    af: float | None
    qual: float | None
    sample_id: str


@dataclass
class SVGeneOverlap:
    sv: StructuralVariant
    gene: str
    overlap_type: str
    hi_score: float | None
    ts_score: float | None
    exon_overlap: bool


@dataclass
class SVPriorityScore:
    sv: StructuralVariant
    gene_overlaps: list[SVGeneOverlap]
    composite_score: float
    tier: SVPriorityTier
    decipher_match: bool
    clinvar_sv_match: bool


@dataclass
class PipelineResult:
    prioritised_svs: list[SVPriorityScore]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
