"""Data models for Metabolic Pathway Balancer."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class BottleneckType(Enum):
    FLUX_LIMITING = "flux_limiting"
    TOXIC_INTERMEDIATE = "toxic_intermediate"
    COFACTOR_IMBALANCE = "cofactor_imbalance"
    COMPETING_PATHWAY = "competing_pathway"


@dataclass
class Reaction:
    reaction_id: str
    enzyme: str
    substrates: list[str]
    products: list[str]
    stoichiometry: dict[str, float]
    kcat: float | None
    km: float | None


@dataclass
class PathwayBottleneck:
    reaction: Reaction
    bottleneck_type: BottleneckType
    flux_ratio: float
    recommendation: str


@dataclass
class BalancedPathway:
    pathway_id: str
    reactions: list[Reaction]
    bottlenecks: list[PathwayBottleneck]
    predicted_yield: float | None
    predicted_productivity: float | None
    cofactor_balance: dict[str, float]


@dataclass
class PipelineResult:
    balanced_pathway: BalancedPathway | None
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
