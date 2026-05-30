"""Data models for Lineage-Specific Gene Finder."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class LSGOrigin(Enum):
    DE_NOVO = "de_novo"
    ORPHAN = "orphan"
    HORIZONTAL_TRANSFER = "horizontal_transfer"
    RAPID_DIVERGENCE = "rapid_divergence"
    UNKNOWN = "unknown"


@dataclass
class LineageSpecificGene:
    gene_id: str
    species: str
    lineage: str
    inferred_origin: LSGOrigin
    age_estimate_mya: float | None
    expression_evidence: bool
    n_exons: int | None
    functional_annotation: str | None


@dataclass
class PipelineResult:
    lineage_specific_genes: list[LineageSpecificGene]
    n_genes_screened: int
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
