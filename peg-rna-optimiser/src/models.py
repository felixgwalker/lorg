"""Data models for pegRNA Optimiser."""

from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class PegRNACandidate:
    pbs_length: int
    rt_length: int
    pbs_sequence: str
    rt_sequence: str
    spacer: str
    efficiency_score: float
    synthesis_complexity: float
    feature_scores: dict[str, float] = field(default_factory=dict)
    is_pareto_optimal: bool = False


@dataclass
class OptimisationResult:
    candidates: list[PegRNACandidate]
    pareto_front: list[PegRNACandidate]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
