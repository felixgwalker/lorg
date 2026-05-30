"""Data models for Safe Harbour Integration Finder."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class HarbourTier(Enum):
    VALIDATED = "validated"
    CANDIDATE = "candidate"
    SPECULATIVE = "speculative"


@dataclass
class IntegrationSite:
    chromosome: str
    start: int
    end: int
    strand: str = "."
    distance_to_oncogene: int | None = None
    distance_to_regulatory: int | None = None
    repeat_density: float | None = None
    tier: HarbourTier = HarbourTier.CANDIDATE
    notes: list[str] = field(default_factory=list)


@dataclass
class SafeHarbourCandidate:
    site: IntegrationSite
    composite_score: float
    known_harbour_match: str | None = None


@dataclass
class PipelineResult:
    candidates: list[SafeHarbourCandidate]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
