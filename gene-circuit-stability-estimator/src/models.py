"""Data models for Gene Circuit Stability Estimator."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class CircuitBehaviour(Enum):
    STABLE_ON = "stable_on"
    STABLE_OFF = "stable_off"
    OSCILLATING = "oscillating"
    BISTABLE = "bistable"
    UNSTABLE = "unstable"


@dataclass
class CircuitNode:
    node_id: str
    protein: str
    basal_expression: float
    degradation_rate: float
    hill_coefficient: float


@dataclass
class CircuitEdge:
    source: str
    target: str
    interaction_type: str
    k_half: float
    max_effect: float


@dataclass
class StabilityAnalysis:
    circuit_id: str
    steady_states: list[dict[str, float]]
    behaviour: CircuitBehaviour
    lyapunov_exponent: float | None
    period_hours: float | None
    basin_of_attraction_size: float | None
    robustness_score: float


@dataclass
class PipelineResult:
    analysis: StabilityAnalysis | None
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
