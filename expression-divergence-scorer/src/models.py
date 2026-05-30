"""Data models for Expression Divergence Scorer."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class DivergenceMetric(Enum):
    LOG_FOLD_CHANGE = "log_fold_change"
    TAU = "tau"
    JSI = "JSI"
    EUCLIDEAN = "euclidean"


@dataclass
class GeneExpressionProfile:
    gene_id: str
    species: str
    tpm_values: list[float]
    tissue_labels: list[str]
    tau: float | None


@dataclass
class ExpressionDivergenceScore:
    gene_id: str
    ortholog_a: str
    ortholog_b: str
    species_a: str
    species_b: str
    metric: DivergenceMetric
    score: float
    significant: bool
    p_value: float | None


@dataclass
class PipelineResult:
    divergence_scores: list[ExpressionDivergenceScore]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
