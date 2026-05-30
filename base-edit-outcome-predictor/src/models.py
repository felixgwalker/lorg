"""Data models for Base Edit Outcome Predictor."""

from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class BaseEditSite:
    position_in_window: int
    base: str
    trinucleotide_context: str
    edit_probability: float
    is_bystander: bool


@dataclass
class EditingWindowResult:
    target_id: str
    spacer: str
    editor: str
    window_start: int
    window_end: int
    editable_sites: list[BaseEditSite]
    primary_product_probability: float
    bystander_probability: float
    indel_frequency: float


@dataclass
class PipelineResult:
    outcomes: list[EditingWindowResult]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
