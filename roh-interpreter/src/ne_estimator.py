"""Ne estimation from ROH length class proportions."""

from __future__ import annotations

from dataclasses import dataclass

from src.froh_calculator import FROHResult, AUTOSOMAL_GENOME_LENGTH


@dataclass
class NeEstimate:
    """Ne trajectory point for one individual."""
    individual_id: str
    time_class: str
    generations_ago: float
    ne: float


def estimate_ne(
    froh_results: list[FROHResult],
    generation_time_years: float = 6.0,
    genome_length: int = AUTOSOMAL_GENOME_LENGTH,
) -> list[NeEstimate]:
    """Estimate Ne from ROH length classes using Ne ≈ 1/(2*g*f).

    Three time windows correspond to short/medium/long ROH classes:
      - LONG  (>10 Mb)  → recent bottleneck (~50 gen)
      - MEDIUM (1-10 Mb) → intermediate (~200 gen)
      - SHORT (100kb-1Mb) → ancient (~1000 gen)
    """
    estimates: list[NeEstimate] = []
    for fr in froh_results:
        for time_class, bp, gen_ago in [
            ("LONG", fr.bp_long, 50.0),
            ("MEDIUM", fr.bp_medium, 200.0),
            ("SHORT", fr.bp_short, 1000.0),
        ]:
            f = bp / genome_length if genome_length > 0 else 0.0
            if f > 0:
                ne = 1.0 / (2.0 * gen_ago * f)
            else:
                ne = float("inf")
            ne = min(ne, 1_000_000.0)
            estimates.append(NeEstimate(
                individual_id=fr.individual_id,
                time_class=time_class,
                generations_ago=gen_ago,
                ne=round(ne, 1),
            ))
    return estimates
