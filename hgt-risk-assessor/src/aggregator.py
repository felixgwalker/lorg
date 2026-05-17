"""Weighted aggregation of signal scores into a single risk index."""

from src.config import FLAT_THRESHOLDS, FLAT_TOP
from src.models import AggregationResult, RiskLevel, SignalResult

SIGNAL_WEIGHTS: dict[str, float] = {
    "is_proximity":  0.25,
    "conjugative":   0.25,
    "integron":      0.20,
    "gc_content":    0.15,
    "prophage":      0.15,
}


def classify(risk_index: float) -> RiskLevel:
    """Classify a flat risk index using configurable thresholds (see config.FLAT_THRESHOLDS)."""
    for upper, label in FLAT_THRESHOLDS:
        if risk_index < upper:
            return RiskLevel(label)
    return RiskLevel(FLAT_TOP)


def aggregate(signal_results: list[SignalResult]) -> AggregationResult:
    """
    Compute risk index from signal results.

    Skipped signals are excluded; the weights of the remaining active signals
    are re-normalised so they sum to 1.0.  This avoids artificially deflating
    the score when databases are unavailable, but the report must clearly note
    how many signals contributed.
    """
    active = [s for s in signal_results if not s.skipped and s.score is not None]
    skipped = [s.signal_name for s in signal_results if s.skipped or s.score is None]

    if not active:
        return AggregationResult(
            signal_results=signal_results,
            risk_index=0.0,
            risk_level=RiskLevel.LOW,
            active_weight_sum=0.0,
            skipped_signals=skipped,
        )

    active_weight_sum = sum(s.weight for s in active)
    risk_index = sum(
        s.score * (s.weight / active_weight_sum)
        for s in active
    )
    risk_index = max(0.0, min(1.0, risk_index))     # clamp for floating-point safety

    return AggregationResult(
        signal_results=signal_results,
        risk_index=risk_index,
        risk_level=classify(risk_index),
        active_weight_sum=active_weight_sum,
        skipped_signals=skipped,
    )
