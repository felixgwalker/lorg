"""Data models for Synthetic Promoter Designer."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class PromoterType(Enum):
    CONSTITUTIVE = "constitutive"
    INDUCIBLE = "inducible"
    TISSUE_SPECIFIC = "tissue_specific"
    CELL_CYCLE_REGULATED = "cell_cycle_regulated"


@dataclass
class TFBSInsert:
    tf_name: str
    motif: str
    position: int
    strand: str
    spacing_to_next: int | None


@dataclass
class SyntheticPromoter:
    promoter_id: str
    sequence: str
    promoter_type: PromoterType
    predicted_strength: float
    gc_content: float
    tfbs_inserts: list[TFBSInsert]
    core_elements: list[str]
    n_variants: int


@dataclass
class PipelineResult:
    designs: list[SyntheticPromoter]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
