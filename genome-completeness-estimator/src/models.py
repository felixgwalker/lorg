"""Data models for Genome Completeness Estimator."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class BUSCOStatus(Enum):
    COMPLETE_SINGLE = "complete_single"
    COMPLETE_DUPLICATE = "complete_duplicate"
    FRAGMENTED = "fragmented"
    MISSING = "missing"


@dataclass
class BUSCOResult:
    busco_id: str
    status: BUSCOStatus
    contig: str | None
    start: int | None
    end: int | None
    score: float | None
    length: int | None


@dataclass
class CompletenessStats:
    lineage: str
    n_buscos_assessed: int
    complete_single: int
    complete_duplicate: int
    fragmented: int
    missing: int
    complete_fraction: float
    fragmented_fraction: float
    missing_fraction: float


@dataclass
class PipelineResult:
    stats: CompletenessStats | None
    busco_results: list[BUSCOResult]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
