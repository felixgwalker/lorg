"""Data models for Gene Model Validator."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class ValidationError(Enum):
    NO_START_CODON = "no_start_codon"
    NO_STOP_CODON = "no_stop_codon"
    INTERNAL_STOP = "internal_stop"
    NON_CANONICAL_SPLICE = "non_canonical_splice"
    FRAME_DISRUPTION = "frame_disruption"
    OVERLAPPING_EXONS = "overlapping_exons"
    SHORT_INTRON = "short_intron"


@dataclass
class GeneModelValidation:
    gene_id: str
    transcript_id: str
    chrom: str
    n_exons: int
    cds_length_bp: int
    errors: list[ValidationError]
    warnings: list[str]
    is_valid: bool
    protein_sequence: str | None


@dataclass
class PipelineResult:
    validations: list[GeneModelValidation]
    n_valid: int
    n_invalid: int
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
