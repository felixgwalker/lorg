"""Data models for Repeat Element Classifier."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class RepeatClass(Enum):
    LINE = "LINE"
    SINE = "SINE"
    LTR = "LTR"
    DNA_TRANSPOSON = "DNA_transposon"
    SATELLITE = "satellite"
    SIMPLE_REPEAT = "simple_repeat"
    LOW_COMPLEXITY = "low_complexity"
    UNKNOWN = "unknown"


@dataclass
class RepeatElement:
    chrom: str
    start: int
    end: int
    element_id: str
    repeat_class: RepeatClass
    repeat_family: str | None
    divergence_percent: float
    length_bp: int
    strand: str


@dataclass
class RepeatLandscape:
    repeat_class: RepeatClass
    n_elements: int
    total_length_bp: int
    mean_divergence: float
    fraction_of_genome: float


@dataclass
class PipelineResult:
    elements: list[RepeatElement]
    landscape: list[RepeatLandscape]
    total_repeat_fraction: float
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
