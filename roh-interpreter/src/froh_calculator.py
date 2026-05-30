"""FROH computation per individual."""

from __future__ import annotations

from dataclasses import dataclass

from src.roh_detector import ROHSegment


AUTOSOMAL_GENOME_LENGTH = 2_700_000_000


@dataclass
class FROHResult:
    """FROH and length-class summary for one individual."""
    individual_id: str
    total_roh_bp: int
    n_roh: int
    froh: float
    froh_short: float
    froh_medium: float
    froh_long: float
    n_short: int
    n_medium: int
    n_long: int
    bp_short: int
    bp_medium: int
    bp_long: int


def calculate_froh(
    roh_segments: list[ROHSegment],
    total_genome_bp: int = AUTOSOMAL_GENOME_LENGTH,
) -> FROHResult:
    """Compute genome-wide FROH and per-length-class FROH.

    FROH = sum(ROH lengths) / total_genome_bp

    Also returns froh_short, froh_medium, froh_long — the contribution of
    each length class to the overall inbreeding coefficient.

    Parameters
    ----------
    roh_segments:
        All ROH segments for a *single* individual (any mix of chromosomes).
    total_genome_bp:
        Autosomal genome length used as denominator.  Defaults to the human
        autosomal length (~2.7 Gb).
    """
    segs = roh_segments
    total_bp = sum(s.length_bp for s in segs)
    bp_short = sum(s.length_bp for s in segs if s.length_class == "short")
    bp_medium = sum(s.length_bp for s in segs if s.length_class == "medium")
    bp_long = sum(s.length_bp for s in segs if s.length_class == "long")

    denom = total_genome_bp if total_genome_bp > 0 else 1
    froh = total_bp / denom
    froh_short = bp_short / denom
    froh_medium = bp_medium / denom
    froh_long = bp_long / denom

    individual_id = segs[0].individual_id if segs else "unknown"

    return FROHResult(
        individual_id=individual_id,
        total_roh_bp=total_bp,
        n_roh=len(segs),
        froh=round(froh, 6),
        froh_short=round(froh_short, 6),
        froh_medium=round(froh_medium, 6),
        froh_long=round(froh_long, 6),
        n_short=sum(1 for s in segs if s.length_class == "short"),
        n_medium=sum(1 for s in segs if s.length_class == "medium"),
        n_long=sum(1 for s in segs if s.length_class == "long"),
        bp_short=bp_short,
        bp_medium=bp_medium,
        bp_long=bp_long,
    )


def compute_froh(
    individual_id: str,
    segments: list[ROHSegment],
    genome_length: int = AUTOSOMAL_GENOME_LENGTH,
) -> FROHResult:
    """Compute FROH for *individual_id* from a mixed-individual segment list.

    Filters *segments* to those belonging to *individual_id*, then delegates
    to :func:`calculate_froh`.
    """
    ind_segs = [s for s in segments if s.individual_id == individual_id]
    result = calculate_froh(ind_segs, total_genome_bp=genome_length)
    # Ensure the individual_id is set correctly even if ind_segs is empty
    return FROHResult(
        individual_id=individual_id,
        total_roh_bp=result.total_roh_bp,
        n_roh=result.n_roh,
        froh=result.froh,
        froh_short=result.froh_short,
        froh_medium=result.froh_medium,
        froh_long=result.froh_long,
        n_short=result.n_short,
        n_medium=result.n_medium,
        n_long=result.n_long,
        bp_short=result.bp_short,
        bp_medium=result.bp_medium,
        bp_long=result.bp_long,
    )
