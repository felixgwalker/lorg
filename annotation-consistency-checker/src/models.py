"""Data models for Annotation Consistency Checker."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class InconsistencyType(Enum):
    DUPLICATE_GENE_ID = "duplicate_gene_id"
    OVERLAPPING_FEATURES = "overlapping_features"
    MISSING_CDS = "missing_cds"
    EXON_OUTSIDE_TRANSCRIPT = "exon_outside_transcript"
    TRANSCRIPT_OUTSIDE_GENE = "transcript_outside_gene"
    STRAND_MISMATCH = "strand_mismatch"
    NEGATIVE_COORDINATES = "negative_coordinates"
    CHROMOSOME_NOT_IN_FASTA = "chromosome_not_in_fasta"


@dataclass
class AnnotationIssue:
    feature_id: str
    feature_type: str
    chrom: str
    start: int
    end: int
    issue_type: InconsistencyType
    description: str


@dataclass
class PipelineResult:
    issues: list[AnnotationIssue]
    n_genes_checked: int
    n_transcripts_checked: int
    passes_qc: bool
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
