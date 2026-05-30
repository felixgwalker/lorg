"""Data models for Ancestral Gene Content Reconstructor."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class AncestralState(Enum):
    PRESENT = "present"
    ABSENT = "absent"
    UNCERTAIN = "uncertain"


@dataclass
class AncestralNode:
    node_id: str
    label: str | None
    is_leaf: bool


@dataclass
class AncestralGeneState:
    gene_id: str
    node: AncestralNode
    state: AncestralState
    posterior_probability: float
    gains_on_branch: int
    losses_on_branch: int


@dataclass
class AncestralGenome:
    node: AncestralNode
    gene_states: list[AncestralGeneState]
    estimated_gene_count: int
    confidence: float


@dataclass
class PipelineResult:
    ancestral_genomes: list[AncestralGenome]
    n_genes_assessed: int
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
