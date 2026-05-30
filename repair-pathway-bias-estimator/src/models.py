"""Data models for Repair Pathway Bias Estimator."""

from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class MicrohomologyMatch:
    sequence: str
    length: int
    left_position: int
    right_position: int
    deletion_size: int
    gc_fraction: float
    mh_score: float
    predicted_frequency: float


@dataclass
class RepairPathwayResult:
    cell_type: str
    nhej_probability: float
    mmej_probability: float
    hdr_probability: float
    mmej_products: list[MicrohomologyMatch]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
