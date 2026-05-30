"""Data models for PAM Flexibility Predictor."""

from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class PAMSite:
    position: int
    strand: str
    pam_sequence: str
    score: float


@dataclass
class CasVariantScore:
    variant_name: str
    pam_motif: str
    site_count: int
    density_per_kb: float
    sites: list[PAMSite] = field(default_factory=list)


@dataclass
class PAMFlexibilityResult:
    variant_scores: list[CasVariantScore]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
