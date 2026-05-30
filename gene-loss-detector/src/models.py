"""Data models for Gene Loss Detector."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class LossType(Enum):
    COMPLETE_DELETION = "complete_deletion"
    PSEUDOGENISATION = "pseudogenisation"
    TRUNCATION = "truncation"
    EXON_LOSS = "exon_loss"


class LossConfidence(Enum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


@dataclass
class GeneLoss:
    gene_id: str
    species_with_gene: list[str]
    species_without_gene: list[str]
    loss_type: LossType
    confidence: LossConfidence
    branch_of_loss: str | None
    inactivating_mutations: list[str] = field(default_factory=list)


@dataclass
class PipelineResult:
    gene_losses: list[GeneLoss]
    n_genes_assessed: int
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
