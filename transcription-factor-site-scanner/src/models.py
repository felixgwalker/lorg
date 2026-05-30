"""Data models for Transcription Factor Site Scanner."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class StrandSense(Enum):
    FORWARD = "+"
    REVERSE = "-"
    BOTH = "."


@dataclass
class PWMModel:
    tf_name: str
    matrix_id: str
    ic_content: float
    consensus: str
    width: int


@dataclass
class TFBSHit:
    chrom: str
    start: int
    end: int
    strand: StrandSense
    tf_name: str
    matrix_id: str
    score: float
    p_value: float
    matched_sequence: str


@dataclass
class PipelineResult:
    hits: list[TFBSHit]
    n_sequences_scanned: int
    n_pwms_used: int
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
