"""Data models for Off-Target Cluster Detector."""

from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class OffTargetSite:
    chromosome: str
    position: int
    strand: str
    cfd_score: float
    mismatches: int


@dataclass
class OffTargetCluster:
    cluster_id: str
    chromosome: str
    start: int
    end: int
    n_sites: int
    mean_cfd_score: float
    max_cfd_score: float
    sites: list[OffTargetSite] = field(default_factory=list)
    annotation: str = ""


@dataclass
class ClusterResult:
    clusters: list[OffTargetCluster]
    total_sites: int
    clustered_fraction: float
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
