"""Data models for Prime Edit Efficiency Predictor."""

from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class PegRNASpec:
    design_id: str
    spacer: str
    pbs: str
    rt_template: str
    nick_guide: str | None = None


@dataclass
class EfficiencyFeatures:
    pbs_gc: float
    rt_length: int
    rt_gc: float
    nick_distance: int | None
    spacer_tm: float
    rt_mfe_approx: float


@dataclass
class EfficiencyPrediction:
    design_id: str
    efficiency_score: float
    efficiency_band: str
    features: EfficiencyFeatures


@dataclass
class PipelineResult:
    predictions: list[EfficiencyPrediction]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
