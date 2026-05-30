"""Data models for Microhomology Repair Predictor."""

from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class MHMatch:
    sequence: str
    length: int
    left_start: int
    right_start: int
    deletion_size: int
    gc_fraction: float
    mh_score: float
    predicted_frequency: float


@dataclass
class MicrohomologyResult:
    cut_position: int
    mh_products: list[MHMatch]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
