"""Data models for Biosafety Risk Assessor."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class BiosafetyLevel(Enum):
    BSL1 = "BSL-1"
    BSL2 = "BSL-2"
    BSL3 = "BSL-3"
    BSL4 = "BSL-4"
    UNKNOWN = "unknown"


class RiskCategory(Enum):
    HORIZONTAL_GENE_TRANSFER = "horizontal_gene_transfer"
    ANTIBIOTIC_RESISTANCE = "antibiotic_resistance"
    VIRULENCE_FACTOR = "virulence_factor"
    TOXIN = "toxin"
    PATHOGEN_HOMOLOGY = "pathogen_homology"
    REPLICATION_ENHANCER = "replication_enhancer"


@dataclass
class BiosafetyFlag:
    category: RiskCategory
    description: str
    severity: str
    evidence: str
    gene_hit: str | None


@dataclass
class BiosafetyAssessment:
    sequence_id: str
    recommended_bsl: BiosafetyLevel
    composite_risk_score: float
    flags: list[BiosafetyFlag]
    containment_recommendations: list[str]
    passes_screening: bool


@dataclass
class PipelineResult:
    assessments: list[BiosafetyAssessment]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
