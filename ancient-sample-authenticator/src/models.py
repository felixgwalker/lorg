"""Data models for Ancient Sample Authenticator."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class AuthenticationVerdict(Enum):
    AUTHENTIC = "authentic"
    LIKELY_AUTHENTIC = "likely_authentic"
    UNCERTAIN = "uncertain"
    MODERN_CONTAMINATION = "modern_contamination"
    FAILED = "failed"


@dataclass
class AuthenticationCriterion:
    name: str
    value: float | bool | str
    passed: bool
    threshold: float | str | None
    weight: float


@dataclass
class AuthenticationResult:
    sample_id: str
    verdict: AuthenticationVerdict
    composite_score: float
    criteria: list[AuthenticationCriterion]
    contamination_estimate: float | None
    mean_fragment_length: float | None
    ct_rate_5prime: float | None


@dataclass
class PipelineResult:
    result: AuthenticationResult | None
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
