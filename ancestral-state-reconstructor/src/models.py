"""Data models for Ancestral State Reconstructor."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class ReconstructionMethod(Enum):
    PARSIMONY = "parsimony"
    ML = "maximum_likelihood"
    BAYESIAN = "bayesian"


@dataclass
class AncestralSiteState:
    site: int
    node_id: str
    ancestral_base: str
    posterior_probability: float
    method: ReconstructionMethod


@dataclass
class AncestralSequence:
    node_id: str
    node_label: str | None
    sequence: str
    n_uncertain_sites: int
    method: ReconstructionMethod


@dataclass
class PipelineResult:
    ancestral_sequences: list[AncestralSequence]
    site_states: list[AncestralSiteState]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
