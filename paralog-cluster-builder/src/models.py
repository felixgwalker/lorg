"""Data models for Paralog Cluster Builder."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class ParalogType(Enum):
    SEGMENTAL_DUPLICATION = "segmental_duplication"
    TANDEM = "tandem"
    DISPERSED = "dispersed"
    RETROPOSED = "retroposed"


@dataclass
class ParalogPair:
    gene_a: str
    gene_b: str
    species: str
    sequence_identity: float
    query_coverage: float
    subject_coverage: float
    synonymous_distance: float | None
    paralog_type: ParalogType | None


@dataclass
class ParalogCluster:
    cluster_id: str
    genes: list[str]
    species: str
    n_members: int
    mean_identity: float
    paralog_type: ParalogType | None
    expansion_ratio: float


@dataclass
class PipelineResult:
    clusters: list[ParalogCluster]
    pairs: list[ParalogPair]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
