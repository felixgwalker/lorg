"""Data models for Prime Edit Design Assistant."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class EditType(Enum):
    SUBSTITUTION = "substitution"
    INSERTION = "insertion"
    DELETION = "deletion"


@dataclass
class EditSpec:
    position: int
    ref_allele: str
    alt_allele: str
    edit_type: EditType = EditType.SUBSTITUTION


@dataclass
class PegRNADesign:
    design_id: str
    spacer: str
    pbs: str
    rt_template: str
    pam: str
    pam_position: int
    strand: str
    pbs_gc: float
    rt_length: int
    on_target_score: float
    nick_guide: str | None = None


@dataclass
class PipelineResult:
    designs: list[PegRNADesign]
    edit_spec: EditSpec | None
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
