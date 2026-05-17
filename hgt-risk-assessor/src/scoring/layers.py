"""
Three-layer HGT risk score computation.

Takes the five existing SignalResult objects (from the flat pipeline) plus a
QuerySequence and HostProfile, runs the new feature extractors, and assembles
the Transfer Opportunity / Establishment / Consequence layer scores and the
combined HGT Risk Index.

Missing features are handled by re-normalising weights over available
features only — identical in principle to how the flat aggregator handles
missing BLAST databases.
"""

from __future__ import annotations

import logging
from typing import Optional

from src.config import (
    LAYER_FEATURE_WEIGHTS,
    SCORE_BANDS,
    SCORE_BAND_TOP,
    WEIGHT_PROFILES,
    DEFAULT_WEIGHT_PROFILE,
)
from src.models import (
    FeatureResult,
    HostProfile,
    LayerResult,
    QuerySequence,
    ScoreBand,
    SignalResult,
    ThreeLayerResult,
)
from src.scoring import features as F
from src.scoring.explanation import generate_explanation

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Score-band classifier
# ---------------------------------------------------------------------------

def classify_band(index: float) -> ScoreBand:
    """Map a risk index in [0, 1] to a ScoreBand using configurable thresholds."""
    for upper, label in SCORE_BANDS:
        if index < upper:
            return ScoreBand(label)
    return ScoreBand(SCORE_BAND_TOP)


# ---------------------------------------------------------------------------
# Layer aggregation
# ---------------------------------------------------------------------------

def _aggregate_layer(features: list[FeatureResult]) -> LayerResult:
    """
    Compute a layer score from a list of FeatureResult objects.

    Weights of unavailable features are excluded and the remaining weights
    are re-normalised so they sum to 1.0.  The layer score is the weighted
    average of available feature scores.
    """
    available = [f for f in features if f.available and f.score is not None]
    total     = len(features)

    if not available:
        layer_name = features[0].layer if features else "unknown"
        return LayerResult(
            layer_name=layer_name,
            layer_score=0.0,
            feature_results=features,
            active_weight_sum=0.0,
            completeness=0.0,
        )

    active_weight_sum = sum(f.weight for f in available)
    if active_weight_sum == 0:
        layer_score = 0.0
    else:
        layer_score = sum(
            f.score * (f.weight / active_weight_sum)
            for f in available
        )
    layer_score = max(0.0, min(1.0, layer_score))

    return LayerResult(
        layer_name=available[0].layer,
        layer_score=layer_score,
        feature_results=features,
        active_weight_sum=active_weight_sum,
        completeness=len(available) / total if total > 0 else 0.0,
    )


# ---------------------------------------------------------------------------
# Signal score lookup helper
# ---------------------------------------------------------------------------

def _sig_score(signals: list[SignalResult], name: str) -> Optional[float]:
    """Return the score for a named signal, or None if skipped/unavailable."""
    for s in signals:
        if s.signal_name == name:
            return s.score if not s.skipped else None
    return None


def _sig_evidence(signals: list[SignalResult], name: str) -> dict:
    for s in signals:
        if s.signal_name == name:
            return s.evidence
    return {}


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def compute_three_layer(
    signal_results: list[SignalResult],
    query: QuerySequence,
    host: HostProfile,
    weight_profile_name: str = DEFAULT_WEIGHT_PROFILE,
    donor_taxon: Optional[str] = None,
    recipient_taxon: Optional[str] = None,
    flanking_sequence: Optional[str] = None,
    **kwargs,
) -> ThreeLayerResult:
    """
    Compute the three-layer HGT Risk Index from existing signals and new features.

    Parameters
    ----------
    signal_results     : output of the flat pipeline's five signal runners
    query              : parsed query sequence
    host               : resolved host profile
    weight_profile_name: one of the keys in config.WEIGHT_PROFILES
    donor_taxon        : optional donor organism name (for taxonomic distance)
    recipient_taxon    : optional recipient organism name (defaults to host.identifier)
    flanking_sequence  : optional flanking sequence string (for repeat density)
    """
    if weight_profile_name not in WEIGHT_PROFILES:
        logger.warning(
            f"Unknown weight profile {weight_profile_name!r}; "
            f"using '{DEFAULT_WEIGHT_PROFILE}'."
        )
        weight_profile_name = DEFAULT_WEIGHT_PROFILE
    profile = WEIGHT_PROFILES[weight_profile_name]

    # Use host as recipient if recipient_taxon not explicitly provided
    effective_recipient = recipient_taxon or host.identifier

    # ------------------------------------------------------------------
    # Layer A: Transfer Opportunity
    # ------------------------------------------------------------------
    transfer_features: list[FeatureResult] = [
        F.compute_is_element_match(signal_score=_sig_score(signal_results, "is_proximity")),
        F.compute_integron_association(signal_score=_sig_score(signal_results, "integron")),
        F.compute_conjugative_element(signal_score=_sig_score(signal_results, "conjugative")),
        F.compute_plasmid_context(),
        F.compute_transposase_proximity(),
        F.compute_repeat_density(query=query, flanking_sequence=flanking_sequence),
    ]
    transfer_layer = _aggregate_layer(transfer_features)

    # ------------------------------------------------------------------
    # Layer B: Establishment
    # ------------------------------------------------------------------
    gc_signal_score = _sig_score(signal_results, "gc_content")
    gc_evidence     = _sig_evidence(signal_results, "gc_content")
    establishment_features: list[FeatureResult] = [
        F.compute_gc_deviation(signal_score=gc_signal_score, evidence=gc_evidence),
        F.compute_codon_usage_distance(query=query, host=host),
        F.compute_taxonomic_distance(
            donor_taxon=donor_taxon,
            recipient_taxon=effective_recipient,
        ),
        F.compute_promoter_plausibility(query=query),
        F.compute_sequence_complexity(query=query),
    ]
    establishment_layer = _aggregate_layer(establishment_features)

    # ------------------------------------------------------------------
    # Layer C: Functional Consequence
    # ------------------------------------------------------------------
    consequence_features: list[FeatureResult] = [
        F.compute_prophage_context(signal_score=_sig_score(signal_results, "prophage")),
        F.compute_amr_content(),
        F.compute_virulence_flags(),
        F.compute_gene_completeness(query=query),
        F.compute_payload_count(query=query),
    ]
    consequence_layer = _aggregate_layer(consequence_features)

    # ------------------------------------------------------------------
    # Combined HGT Risk Index
    # ------------------------------------------------------------------
    # Only include layers with at least one active feature; re-normalise.
    layers = [
        ("transfer_opportunity", transfer_layer),
        ("establishment",        establishment_layer),
        ("consequence",          consequence_layer),
    ]
    active_layers = [(k, lr) for k, lr in layers if lr.active_weight_sum > 0]
    active_profile_sum = sum(profile[k] for k, _ in active_layers)

    if active_profile_sum == 0:
        hgt_risk_index = 0.0
    else:
        hgt_risk_index = sum(
            lr.layer_score * (profile[k] / active_profile_sum)
            for k, lr in active_layers
        )
    hgt_risk_index = max(0.0, min(1.0, hgt_risk_index))

    all_features = transfer_features + establishment_features + consequence_features
    total_features   = len(all_features)
    active_features  = sum(1 for f in all_features if f.available)
    overall_completeness = active_features / total_features if total_features > 0 else 0.0

    score_band = classify_band(hgt_risk_index)

    explanation, top_contributors, risk_reducers, missing_important = generate_explanation(
        transfer_layer=transfer_layer,
        establishment_layer=establishment_layer,
        consequence_layer=consequence_layer,
        hgt_risk_index=hgt_risk_index,
        score_band=score_band,
        overall_completeness=overall_completeness,
    )

    logger.info(
        f"Three-layer HGT Risk Index: {hgt_risk_index:.3f} ({score_band.value})  "
        f"completeness={overall_completeness:.0%}  profile={weight_profile_name}"
    )

    return ThreeLayerResult(
        transfer_layer=transfer_layer,
        establishment_layer=establishment_layer,
        consequence_layer=consequence_layer,
        hgt_risk_index=hgt_risk_index,
        score_band=score_band,
        weight_profile_name=weight_profile_name,
        weight_profile=profile,
        overall_completeness=overall_completeness,
        explanation=explanation,
        top_contributors=top_contributors,
        risk_reducers=risk_reducers,
        missing_important_features=missing_important,
    )
