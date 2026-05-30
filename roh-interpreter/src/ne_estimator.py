"""Ne estimation from ROH length class proportions."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

from src.froh_calculator import FROHResult, AUTOSOMAL_GENOME_LENGTH

if TYPE_CHECKING:
    from src.roh_detector import ROHSegment


@dataclass
class NeEstimate:
    """Per-individual Ne trajectory derived from ROH length class proportions.

    Fields
    ------
    individual_id:
        Sample identifier.
    ne_recent:
        Ne estimated from *long* ROH (>1 Mb), reflecting the last ~2–5
        generations.
    ne_moderate:
        Ne estimated from *medium* ROH (100 kb–1 Mb), reflecting ~10–50
        generations ago.
    ne_ancient:
        Ne estimated from *short* ROH (<100 kb), reflecting ~50–200
        generations ago.
    generation_time_years:
        Generation time in years used for calendar-year conversions.
    ci_low:
        Lower 95 % confidence bound on ne_recent (Poisson approximation).
    ci_high:
        Upper 95 % confidence bound on ne_recent (Poisson approximation).
    """
    individual_id: str
    ne_recent: float
    ne_moderate: float
    ne_ancient: float
    generation_time_years: float
    ci_low: float
    ci_high: float


# Representative generation numbers for the Ne ≈ 100/(2·g·FROH_class) formula
# long  ROH → g_recent    (midpoint of 2–5 range)
# medium ROH → g_moderate  (midpoint of 10–50 range)
# short  ROH → g_ancient   (midpoint of 50–200 range)
_G_RECENT = 3.5      # long ROH: 2–5 generations
_G_MODERATE = 25.0   # medium ROH: 10–50 generations
_G_ANCIENT = 100.0   # short ROH: 50–200 generations

# Cap on Ne to avoid infinity for zero-FROH classes
_NE_MAX = 1_000_000.0


def estimate_ne_from_roh(
    roh_segments: list[ROHSegment],
    generation_time_years: float = 25.0,
    genome_length: int = AUTOSOMAL_GENOME_LENGTH,
    individual_id: str | None = None,
) -> NeEstimate:
    """Estimate Ne trajectory from ROH segments for a single individual.

    Uses the formula::

        Ne(g) ≈ 100 / (2 × g × FROH_class)

    where *g* is a representative generation depth for each ROH class and
    *FROH_class* is the fraction of the genome covered by that ROH class.

    Parameters
    ----------
    roh_segments:
        All ROH segments for a single individual (any chromosomes).  May be
        empty — in that case all Ne values are set to the maximum cap.
    generation_time_years:
        Generation time in years (used only for annotation; does not alter Ne
        arithmetic).
    genome_length:
        Autosomal genome length in bp (denominator for FROH calculation).
    individual_id:
        Explicit sample identifier.  If ``None`` the ID is read from the
        first element of *roh_segments*.  Required when *roh_segments* is
        empty.

    Returns
    -------
    NeEstimate
        Dataclass with ne_recent, ne_moderate, ne_ancient, generation_time_years,
        ci_low, and ci_high.
    """
    from src.froh_calculator import calculate_froh

    if individual_id is None:
        individual_id = roh_segments[0].individual_id if roh_segments else "unknown"
    froh_result = calculate_froh(roh_segments, total_genome_bp=genome_length)

    ne_recent = _ne_from_froh(froh_result.froh_long, _G_RECENT)
    ne_moderate = _ne_from_froh(froh_result.froh_medium, _G_MODERATE)
    ne_ancient = _ne_from_froh(froh_result.froh_short, _G_ANCIENT)

    # 95 % CI on ne_recent using Poisson approximation on the ROH count.
    # Under a Poisson count model: SE(FROH) ≈ sqrt(n_long) / L → propagate to Ne.
    ci_low, ci_high = _ne_recent_ci(
        froh_result.froh_long, froh_result.n_long, genome_length, _G_RECENT
    )

    return NeEstimate(
        individual_id=individual_id,
        ne_recent=round(ne_recent, 1),
        ne_moderate=round(ne_moderate, 1),
        ne_ancient=round(ne_ancient, 1),
        generation_time_years=generation_time_years,
        ci_low=round(ci_low, 1),
        ci_high=round(ci_high, 1),
    )


def _ne_from_froh(froh_class: float, g: float) -> float:
    """Apply Ne ≈ 100 / (2 × g × FROH_class), capped at _NE_MAX."""
    if froh_class <= 0.0:
        return _NE_MAX
    ne = 100.0 / (2.0 * g * froh_class)
    return min(ne, _NE_MAX)


def _ne_recent_ci(
    froh_long: float,
    n_long: int,
    genome_length: int,
    g: float,
) -> tuple[float, float]:
    """Approximate 95 % CI for ne_recent using Poisson count propagation.

    Assumes ROH count follows a Poisson distribution.  The standard error of
    FROH_long is estimated as sqrt(n_long) * mean_roh_length / genome_length,
    where mean_roh_length is back-computed from FROH and n.

    When n_long == 0 the CI is (Ne_max, Ne_max) — effectively unresolved.
    """
    if n_long <= 0 or froh_long <= 0.0:
        return (_NE_MAX, _NE_MAX)

    # mean ROH length (bp), back-computed from FROH and count
    mean_roh_len = (froh_long * genome_length) / n_long
    # SE(total_bp) = sqrt(n_long) * mean_roh_len  (Poisson count model)
    # SE(FROH_long) = SE(total_bp) / genome_length = mean_roh_len * sqrt(n_long) / L
    se_froh = mean_roh_len * math.sqrt(n_long) / genome_length

    # 1.96-sigma band on FROH_long
    froh_lo = max(froh_long - 1.96 * se_froh, 1e-9)
    froh_hi = froh_long + 1.96 * se_froh

    ne_lo = min(_NE_MAX, 100.0 / (2.0 * g * froh_hi))
    ne_hi = min(_NE_MAX, 100.0 / (2.0 * g * froh_lo))
    return (ne_lo, ne_hi)


# ---------------------------------------------------------------------------
# Legacy interface kept for pipeline.py compatibility
# ---------------------------------------------------------------------------

@dataclass
class _NePoint:
    """Legacy Ne trajectory point (internal use)."""
    individual_id: str
    time_class: str
    generations_ago: float
    ne: float


def estimate_ne(
    froh_results: list[FROHResult],
    generation_time_years: float = 6.0,
    genome_length: int = AUTOSOMAL_GENOME_LENGTH,
) -> list[_NePoint]:
    """Estimate Ne from per-individual FROHResult objects (pipeline helper).

    Produces three Ne trajectory points per individual (long/medium/short),
    using the same Ne ≈ 100/(2·g·FROH_class) formula.
    """
    estimates: list[_NePoint] = []
    for fr in froh_results:
        for time_class, froh_class, gen_ago in [
            ("long",   fr.froh_long,   _G_RECENT),
            ("medium", fr.froh_medium, _G_MODERATE),
            ("short",  fr.froh_short,  _G_ANCIENT),
        ]:
            ne = _ne_from_froh(froh_class, gen_ago)
            estimates.append(_NePoint(
                individual_id=fr.individual_id,
                time_class=time_class,
                generations_ago=gen_ago,
                ne=round(ne, 1),
            ))
    return estimates
