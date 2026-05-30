"""Data models for Alternative Splicing Detector."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class SplicingEventType(Enum):
    EXON_SKIPPING = "exon_skipping"
    INTRON_RETENTION = "intron_retention"
    ALT_5_PRIME = "alt_5_prime"
    ALT_3_PRIME = "alt_3_prime"
    MUTUALLY_EXCLUSIVE = "mutually_exclusive"
    ALT_FIRST_EXON = "alt_first_exon"
    ALT_LAST_EXON = "alt_last_exon"


@dataclass
class SplicingEvent:
    event_id: str
    gene_id: str
    event_type: SplicingEventType
    chrom: str
    start: int
    end: int
    psi_condition_a: float
    psi_condition_b: float
    delta_psi: float
    p_value: float
    fdr: float
    significant: bool


@dataclass
class PipelineResult:
    events: list[SplicingEvent]
    n_genes_with_events: int
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
