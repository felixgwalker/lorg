"""Data models for Rare Variant Prioritiser."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class EvidenceCode(Enum):
    ULTRA_RARE = "ultra_rare"
    RARE = "rare"
    CONSERVED = "conserved"
    HIGH_CADD = "high_cadd"
    IN_PANEL = "in_panel"
    HPO_MATCH = "hpo_match"
    CONSTRAINED_GENE = "constrained_gene"
    CLINVAR_PATHOGENIC = "clinvar_pathogenic"


class PriorityTier(Enum):
    TIER_1 = "tier_1"
    TIER_2 = "tier_2"
    TIER_3 = "tier_3"
    UNRANKED = "unranked"


@dataclass
class RareVariant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    consequence: str
    gnomad_af: float | None
    cadd_score: float | None
    clinvar_sig: str | None


@dataclass
class PriorityScore:
    variant: RareVariant
    evidence_codes: list[EvidenceCode]
    composite_score: float
    tier: PriorityTier
    hpo_matched_genes: list[str] = field(default_factory=list)


@dataclass
class PipelineResult:
    prioritised_variants: list[PriorityScore]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
