"""Data models for Edit Outcome Simulator."""

from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class IndelPrediction:
    indel_type: str
    indel_size: int
    indel_sequence: str | None
    frequency: float
    frameshift: bool


@dataclass
class EditOutcomeDistribution:
    target_id: str
    guide: str
    cas_variant: str
    indel_outcomes: list[IndelPrediction]
    frameshift_rate: float
    most_frequent_outcome: IndelPrediction | None
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
