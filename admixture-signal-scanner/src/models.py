"""Data models for Admixture Signal Scanner."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class AdmixtureModel(Enum):
    UNSUPERVISED = "unsupervised"
    SUPERVISED = "supervised"


@dataclass
class AncestryComponent:
    component_id: int
    label: str | None
    reference_population: str | None


@dataclass
class SampleAncestry:
    sample_id: str
    population: str | None
    proportions: dict[int, float]
    dominant_component: int
    admixed: bool


@dataclass
class AdmixtureRun:
    k: int
    cross_validation_error: float | None
    log_likelihood: float | None
    components: list[AncestryComponent]


@dataclass
class PipelineResult:
    best_k: int
    runs: list[AdmixtureRun]
    sample_ancestries: list[SampleAncestry]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
