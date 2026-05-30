"""Data models for Chromatin Accessibility Scorer."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class AccessibilityClass(Enum):
    OPEN = "open"
    INTERMEDIATE = "intermediate"
    CLOSED = "closed"


@dataclass
class ChromatinRegion:
    chrom: str
    start: int
    end: int
    peak_id: str
    cell_type: str | None


@dataclass
class AccessibilityScore:
    region: ChromatinRegion
    raw_signal: float
    normalised_score: float
    accessibility_class: AccessibilityClass
    enrichment_over_background: float
    q_value: float | None


@dataclass
class PipelineResult:
    scores: list[AccessibilityScore]
    n_open_regions: int
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
