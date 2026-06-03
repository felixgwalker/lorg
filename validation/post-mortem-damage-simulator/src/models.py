"""Data models for Post-Mortem Damage Simulator."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class DamageModel(Enum):
    BRIGGS = "briggs"
    UNIFORM = "uniform"
    DOUBLE_STRANDED = "double_stranded"
    SINGLE_STRANDED = "single_stranded"


@dataclass
class DamageParameters:
    model: DamageModel
    nick_freq_5prime: float
    nick_freq_3prime: float
    overhang_length: float
    deamination_rate_ds: float
    deamination_rate_ss: float
    mean_fragment_length: int
    fragmentation_lambda: float


@dataclass
class SimulatedRead:
    read_id: str
    sequence: str
    original_sequence: str
    chrom: str
    start: int
    end: int
    strand: str
    n_ct_substitutions: int
    n_ga_substitutions: int


@dataclass
class PipelineResult:
    simulated_reads: list[SimulatedRead]
    damage_parameters: DamageParameters
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
