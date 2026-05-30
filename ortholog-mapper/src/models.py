"""Data models for Ortholog Mapper."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class OrthologType(Enum):
    ONE_TO_ONE = "one_to_one"
    ONE_TO_MANY = "one_to_many"
    MANY_TO_ONE = "many_to_one"
    MANY_TO_MANY = "many_to_many"


class OrthologMethod(Enum):
    RECIPROCAL_BEST_HITS = "reciprocal_best_hits"
    OMA = "OMA"
    INPARANOID = "inparanoid"


@dataclass
class Ortholog:
    query_gene: str
    query_species: str
    target_gene: str
    target_species: str
    ortholog_type: OrthologType
    method: OrthologMethod
    sequence_identity: float | None
    evalue: float | None
    synteny_support: bool


@dataclass
class OrthologGroup:
    group_id: str
    members: list[Ortholog]
    n_species: int
    is_universal: bool


@dataclass
class PipelineResult:
    ortholog_groups: list[OrthologGroup]
    n_query_genes: int
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
