"""Data models for Gene Regulatory Network Builder."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class EdgeType(Enum):
    ACTIVATION = "activation"
    REPRESSION = "repression"
    UNKNOWN = "unknown"


class NetworkMethod(Enum):
    GENIE3 = "GENIE3"
    CORRELATION = "correlation"
    ARACNE = "ARACNE"


@dataclass
class RegEdge:
    regulator: str
    target: str
    weight: float
    edge_type: EdgeType
    method: NetworkMethod
    p_value: float | None


@dataclass
class NetworkNode:
    gene_id: str
    is_tf: bool
    in_degree: int
    out_degree: int
    hub_score: float


@dataclass
class PipelineResult:
    edges: list[RegEdge]
    nodes: list[NetworkNode]
    n_tf_regulators: int
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
