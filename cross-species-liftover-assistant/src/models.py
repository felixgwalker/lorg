"""Data models for Cross-Species Liftover Assistant."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class LiftoverStatus(Enum):
    SUCCESS = "success"
    PARTIAL = "partial"
    FAILED = "failed"
    MULTI_MAPPING = "multi_mapping"


@dataclass
class LiftoverRecord:
    source_chrom: str
    source_start: int
    source_end: int
    source_species: str
    target_chrom: str | None
    target_start: int | None
    target_end: int | None
    target_species: str
    status: LiftoverStatus
    identity: float | None
    strand: str | None
    name: str | None = None


@dataclass
class PipelineResult:
    records: list[LiftoverRecord]
    n_input: int
    n_success: int
    n_failed: int
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
