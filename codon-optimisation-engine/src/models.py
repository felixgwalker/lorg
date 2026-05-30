"""Data models for Codon Optimisation Engine."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class OptimisationStrategy(Enum):
    MOST_FREQUENT = "most_frequent"
    HARMONISED = "harmonised"
    RANDOM_WEIGHTED = "random_weighted"
    CAI_MAXIMISED = "CAI_maximised"


@dataclass
class CodonUsageTable:
    organism: str
    n_codons: int
    codon_frequencies: dict[str, float]
    cai_weights: dict[str, float]


@dataclass
class OptimisedSequence:
    gene_id: str
    original_protein: str
    original_dna: str | None
    optimised_dna: str
    cai_before: float | None
    cai_after: float
    gc_content_before: float | None
    gc_content_after: float
    strategy: OptimisationStrategy
    n_codons_changed: int


@dataclass
class PipelineResult:
    optimised_sequences: list[OptimisedSequence]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
