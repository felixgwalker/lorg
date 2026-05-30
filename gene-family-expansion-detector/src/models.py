"""Data models for Gene Family Expansion Detector."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class ExpansionSignificance(Enum):
    SIGNIFICANT = "significant"
    MARGINAL = "marginal"
    NOT_SIGNIFICANT = "not_significant"


@dataclass
class FamilySizeProfile:
    family_id: str
    size_per_species: dict[str, int]
    mean_size: float
    sd_size: float
    max_species: str
    max_size: int


@dataclass
class FamilyExpansion:
    family: FamilySizeProfile
    expanded_species: str
    fold_expansion: float
    p_value: float
    significance: ExpansionSignificance
    branch_of_expansion: str | None
    functional_annotation: str | None


@dataclass
class PipelineResult:
    expansions: list[FamilyExpansion]
    n_families_assessed: int
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
