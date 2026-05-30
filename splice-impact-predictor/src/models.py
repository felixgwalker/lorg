"""Data models for Splice Impact Predictor."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class SiteType(Enum):
    DONOR = "donor"
    ACCEPTOR = "acceptor"
    BRANCHPOINT = "branchpoint"


class SpliceEffect(Enum):
    DISRUPTION = "disruption"
    CREATION = "creation"
    WEAKENING = "weakening"
    NEUTRAL = "neutral"


@dataclass
class SpliceSite:
    chrom: str
    pos: int
    site_type: SiteType
    exon_id: str
    transcript_id: str
    strand: str
    consensus_seq: str


@dataclass
class SpliceVariant:
    chrom: str
    pos: int
    ref: str
    alt: str
    nearest_site: SpliceSite
    distance_to_site: int


@dataclass
class SpliceScore:
    variant: SpliceVariant
    ref_score: float
    alt_score: float
    delta_score: float
    effect: SpliceEffect
    confidence: float


@dataclass
class PipelineResult:
    splice_scores: list[SpliceScore]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
