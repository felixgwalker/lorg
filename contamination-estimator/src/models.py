"""Data models for Contamination Estimator."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class ContaminationMethod(Enum):
    MT_CONSENSUS = "mt_consensus"
    NUCLEAR_X = "nuclear_x"
    ANGSD = "ANGSD"
    SCHMUTZI = "schmutzi"


@dataclass
class ContaminationEstimate:
    method: ContaminationMethod
    contamination_rate: float
    ci_lower: float
    ci_upper: float
    n_sites: int
    confidence: str


@dataclass
class PipelineResult:
    estimates: list[ContaminationEstimate]
    combined_estimate: float | None
    passes_threshold: bool
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
