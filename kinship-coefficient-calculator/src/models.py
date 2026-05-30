"""Data models for Kinship Coefficient Calculator."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class KinshipMethod(Enum):
    KING = "KING"
    GENOMIC_RELATEDNESS = "genomic_relatedness"
    IBD = "IBD"


class RelationshipClass(Enum):
    IDENTICAL = "identical"
    FIRST_DEGREE = "first_degree"
    SECOND_DEGREE = "second_degree"
    THIRD_DEGREE = "third_degree"
    UNRELATED = "unrelated"


@dataclass
class KinshipPair:
    sample_a: str
    sample_b: str
    kinship_coefficient: float
    ibd0: float | None
    ibd1: float | None
    ibd2: float | None
    relationship_class: RelationshipClass
    method: KinshipMethod


@dataclass
class PipelineResult:
    kinship_pairs: list[KinshipPair]
    n_samples: int
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
