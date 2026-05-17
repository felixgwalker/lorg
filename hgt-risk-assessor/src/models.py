"""Core data structures for the HGT risk assessment pipeline."""

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

# ---------------------------------------------------------------------------
# Three-layer model types (added in v0.2)
# ---------------------------------------------------------------------------


class ScoreBand(Enum):
    """Risk bands for the three-layer HGT Risk Index (configurable in config.py)."""
    LOW       = "low"
    MODERATE  = "moderate"
    HIGH      = "high"
    VERY_HIGH = "very_high"


@dataclass
class FeatureResult:
    """
    Result of a single feature extraction.

    score=None / available=False indicates the feature could not be computed
    (database absent, optional input not provided, placeholder not yet integrated).
    The layer aggregator re-normalises over available features only.
    """
    feature_name: str
    layer: str              # "transfer_opportunity" | "establishment" | "consequence"
    score: Optional[float]  # 0.0–1.0; None if unavailable
    weight: float           # within-layer weight from config.LAYER_FEATURE_WEIGHTS
    available: bool
    evidence: dict
    interpretation: str = ""
    source: str = ""        # "signal_reuse" | "computed" | "placeholder"


@dataclass
class LayerResult:
    """Aggregated score for one of the three risk layers."""
    layer_name: str
    layer_score: float          # 0.0–1.0, re-normalised over available features
    feature_results: list[FeatureResult]
    active_weight_sum: float    # sum of weights of available features
    completeness: float         # available_features / total_features


@dataclass
class ThreeLayerResult:
    """
    Combined output of the three-layer HGT risk model.

    hgt_risk_index is a weighted combination of the three layer scores using
    the selected weight profile.  All component scores and feature contributions
    are preserved for full interpretability.
    """
    transfer_layer: LayerResult
    establishment_layer: LayerResult
    consequence_layer: LayerResult
    hgt_risk_index: float
    score_band: ScoreBand
    weight_profile_name: str
    weight_profile: dict        # {layer_name: weight} for the selected profile
    overall_completeness: float
    explanation: str
    top_contributors: list[str]
    risk_reducers: list[str]
    missing_important_features: list[str]


class RiskLevel(Enum):
    LOW = "Low"
    MEDIUM = "Medium"
    HIGH = "High"
    CRITICAL = "Critical"

    @property
    def css_class(self) -> str:
        return f"risk-{self.value.lower()}"


class InputFormat(Enum):
    FASTA = "fasta"
    VCF = "vcf"


@dataclass
class HostProfile:
    identifier: str
    gc_content: float           # 0.0–1.0
    organism_name: str = ""
    source: str = ""            # "ncbi_entrez" | "user_supplied" | "estimated"


@dataclass
class QuerySequence:
    sequence: str
    identifier: str
    length: int
    gc_content: float           # 0.0–1.0
    source_format: InputFormat


@dataclass
class BlastHit:
    query_id: str
    subject_id: str
    pct_identity: float
    alignment_length: int
    query_coverage: float
    evalue: float
    bit_score: float
    subject_description: str = ""


@dataclass
class AttCSite:
    start: int
    end: int
    sequence: str               # first 20 chars of match
    spacer_length: int
    strand: str                 # "+" | "-"


@dataclass
class PhasterRegion:
    region_id: str
    completeness: str           # "intact" | "questionable" | "incomplete"
    start: int
    end: int
    gc_content: float           # 0.0–1.0
    num_cds: int


@dataclass
class PhasterResult:
    status: str                 # "complete" | "error"
    regions: list[PhasterRegion]
    raw_response: dict


@dataclass
class SignalResult:
    signal_name: str
    score: Optional[float]      # None = skipped
    weight: float
    evidence: dict
    warning: str = ""
    skipped: bool = False

    @property
    def weighted_contribution(self) -> Optional[float]:
        if self.score is None:
            return None
        return self.score * self.weight


@dataclass
class AggregationResult:
    signal_results: list[SignalResult]
    risk_index: float
    risk_level: RiskLevel
    active_weight_sum: float    # sum of weights for non-skipped signals
    skipped_signals: list[str]


@dataclass
class PipelineResult:
    query: QuerySequence
    host: HostProfile
    aggregation: AggregationResult
    run_timestamp: str
    pipeline_version: str
    three_layer: Optional["ThreeLayerResult"] = None   # populated when three-layer model runs
