"""Data models for Conserved Synteny Detector."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class SyntenyOrientation(Enum):
    SAME = "same"
    INVERTED = "inverted"
    MIXED = "mixed"


@dataclass
class SyntenyBlock:
    block_id: str
    species_a: str
    chrom_a: str
    start_a: int
    end_a: int
    species_b: str
    chrom_b: str
    start_b: int
    end_b: int
    orientation: SyntenyOrientation
    n_anchors: int
    coverage_a: float
    coverage_b: float


@dataclass
class SyntenyAnchor:
    gene_a: str
    gene_b: str
    position_a: int
    position_b: int
    block_id: str


@dataclass
class PipelineResult:
    synteny_blocks: list[SyntenyBlock]
    anchors: list[SyntenyAnchor]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
