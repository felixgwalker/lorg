"""Data models for Missense Impact Scorer."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class ImpactClass(Enum):
    BENIGN = "benign"
    LIKELY_BENIGN = "likely_benign"
    UNCERTAIN = "uncertain"
    LIKELY_PATHOGENIC = "likely_pathogenic"
    PATHOGENIC = "pathogenic"


class ConservationTool(Enum):
    PHYLOP = "phyloP"
    GERP = "GERP"
    PHASTCONS = "phastCons"


@dataclass
class MissenseVariant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    transcript: str
    hgvs_p: str
    ref_aa: str
    alt_aa: str
    codon_pos: int


@dataclass
class ImpactScore:
    variant: MissenseVariant
    conservation_score: float
    blosum62_score: float
    physicochemical_delta: float
    structural_context: str | None
    composite_score: float
    impact_class: ImpactClass


@dataclass
class PipelineResult:
    scored_variants: list[ImpactScore]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
