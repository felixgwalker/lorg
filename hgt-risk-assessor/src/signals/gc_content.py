"""Signal 2 — GC content deviation from host genome."""

from src.aggregator import SIGNAL_WEIGHTS
from src.models import HostProfile, QuerySequence, SignalResult

WEIGHT = SIGNAL_WEIGHTS["gc_content"]

# A deviation of ≥25 percentage points maps to a score of 1.0.
# This threshold reflects the empirical observation that horizontally acquired
# regions typically deviate by ≥10–15 pp from host GC; 25 pp is a conservative
# upper bound for maximal concern.
MAX_DEVIATION = 0.25


def run(query: QuerySequence, host: HostProfile, **kwargs) -> SignalResult:
    deviation = abs(query.gc_content - host.gc_content)
    score = min(deviation / MAX_DEVIATION, 1.0)

    return SignalResult(
        signal_name="gc_content",
        score=score,
        weight=WEIGHT,
        evidence={
            "query_gc":        round(query.gc_content * 100, 2),
            "host_gc":         round(host.gc_content * 100, 2),
            "deviation_pct":   round(deviation * 100, 2),
            "max_deviation_pct": MAX_DEVIATION * 100,
        },
    )
