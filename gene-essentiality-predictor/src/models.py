"""Data models for Gene Essentiality Predictor."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class EssentialityClass(Enum):
    ESSENTIAL = "essential"
    CONTEXT_DEPENDENT = "context_dependent"
    NON_ESSENTIAL = "non_essential"
    UNKNOWN = "unknown"


@dataclass
class EssentialityEvidence:
    source: str
    score: float | None
    is_essential: bool | None


@dataclass
class GeneEssentialityPrediction:
    gene_id: str
    gene_symbol: str | None
    composite_score: float
    essentiality_class: EssentialityClass
    loeuf: float | None
    pli: float | None
    crispr_fitness_score: float | None
    rnai_score: float | None
    DepMap_essential_lines: int | None
    evidence: list[EssentialityEvidence] = field(default_factory=list)


@dataclass
class PipelineResult:
    predictions: list[GeneEssentialityPrediction]
    n_essential: int
    n_context_dependent: int
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
