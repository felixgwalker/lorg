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
    n_short: int
    n_medium: int
    n_long: int
    bp_short: int
    bp_medium: int
    bp_long: int


def compute_froh(
    individual_id: str,
    segments: list[ROHSegment],
    genome_length: int = AUTOSOMAL_GENOME_LENGTH,
) -> FROHResult:
    """Compute FROH = sum(ROH lengths) / autosomal genome length."""
    ind_segs = [s for s in segments if s.individual_id == individual_id]
    total_bp = sum(s.length_bp for s in ind_segs)
    n_short = sum(1 for s in ind_segs if s.length_class == "SHORT")
    n_medium = sum(1 for s in ind_segs if s.length_class == "MEDIUM")
    n_long = sum(1 for s in ind_segs if s.length_class == "LONG")
    bp_short = sum(s.length_bp for s in ind_segs if s.length_class == "SHORT")
    bp_medium = sum(s.length_bp for s in ind_segs if s.length_class == "MEDIUM")
    bp_long = sum(s.length_bp for s in ind_segs if s.length_class == "LONG")
    froh = total_bp / genome_length if genome_length > 0 else 0.0
    return FROHResult(
        individual_id=individual_id,
        total_roh_bp=total_bp,
        n_roh=len(ind_segs),
        froh=round(froh, 6),
        n_short=n_short,
        n_medium=n_medium,
        n_long=n_long,
        bp_short=bp_short,
        bp_medium=bp_medium,
        bp_long=bp_long,
    )
