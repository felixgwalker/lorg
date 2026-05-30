"""Data models for Allele Frequency Comparator."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class Population(Enum):
    AFR = "afr"
    AMR = "amr"
    ASJ = "asj"
    EAS = "eas"
    FIN = "fin"
    NFE = "nfe"
    SAS = "sas"
    OTH = "oth"
    ALL = "all"


@dataclass
class PopulationAF:
    population: Population
    ac: int
    an: int
    af: float
    homozygote_count: int | None


@dataclass
class FrequencyComparison:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str | None
    population_afs: dict[str, PopulationAF]
    max_af: float
    min_af: float
    fold_change: float
    most_common_pop: str
    rarest_pop: str
    fst_estimate: float | None


@dataclass
class DifferentialVariant:
    comparison: FrequencyComparison
    is_population_specific: bool
    specificity_population: str | None
    p_value: float | None
    significant: bool


@dataclass
class PipelineResult:
    comparisons: list[DifferentialVariant]
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
