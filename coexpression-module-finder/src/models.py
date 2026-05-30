"""Data models for Coexpression Module Finder."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class ModuleMethod(Enum):
    WGCNA = "WGCNA"
    CLIQUE_BASED = "clique_based"
    KMEANS = "kmeans"


@dataclass
class CoexpressionModule:
    module_id: str
    colour_label: str
    genes: list[str]
    n_genes: int
    hub_gene: str | None
    eigengene_correlation: float | None
    enriched_go_terms: list[str] = field(default_factory=list)


@dataclass
class ModuleTrait:
    module_id: str
    trait: str
    correlation: float
    p_value: float
    significant: bool


@dataclass
class PipelineResult:
    modules: list[CoexpressionModule]
    module_trait_associations: list[ModuleTrait]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
