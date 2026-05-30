"""Data models for UTR Variant Analyser."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class UTRType(Enum):
    FIVE_PRIME = "5prime"
    THREE_PRIME = "3prime"


class UTREffect(Enum):
    UORF_CREATION = "uORF_creation"
    UORF_DISRUPTION = "uORF_disruption"
    KOZAK_CHANGE = "kozak_change"
    POLYA_SIGNAL_DISRUPTION = "polyA_signal_disruption"
    IRES_DISRUPTION = "IRES_disruption"
    NEUTRAL = "neutral"


@dataclass
class UTRVariant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    transcript_id: str
    utr_type: UTRType
    distance_to_cds: int


@dataclass
class uORFChange:
    variant: UTRVariant
    ref_uorf_count: int
    alt_uorf_count: int
    delta_uorf_count: int
    affected_uorf_start: int | None
    kozak_context_ref: str
    kozak_context_alt: str
    kozak_strength_ref: float
    kozak_strength_alt: float


@dataclass
class UTRAnalysisResult:
    variant: UTRVariant
    effects: list[UTREffect]
    uorf_change: uORFChange | None
    polya_signal_disrupted: bool
    impact_summary: str


@dataclass
class PipelineResult:
    results: list[UTRAnalysisResult]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
