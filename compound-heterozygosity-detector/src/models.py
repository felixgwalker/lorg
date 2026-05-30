"""Data models for Compound Heterozygosity Detector."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class PhaseStatus(Enum):
    PHASED = "phased"
    INFERRED_TRIO = "inferred_trio"
    UNPHASED = "unphased"


class CompHetConfidence(Enum):
    CONFIRMED = "confirmed"
    LIKELY = "likely"
    POSSIBLE = "possible"


@dataclass
class PhasedVariant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    consequence: str
    sample_id: str
    haplotype: int | None
    gt: str
    af: float | None


@dataclass
class FamilyMember:
    sample_id: str
    role: str
    affected: bool


@dataclass
class CompHetPair:
    variant_a: PhasedVariant
    variant_b: PhasedVariant
    gene: str
    phase_status: PhaseStatus
    confidence: CompHetConfidence
    trio_consistent: bool | None


@dataclass
class PipelineResult:
    comp_het_pairs: list[CompHetPair]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
