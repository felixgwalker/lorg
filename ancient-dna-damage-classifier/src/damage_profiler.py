"""
Damage profiler: computes position-specific C→T (5' end) and G→A (3' end)
substitution frequencies from a BAM file.

For each read that passes quality filters, the first `n_terminal` aligned
positions are examined for C→T mismatches (5' cytosine deamination) and the
last `n_terminal` positions for G→A mismatches (3' cytosine deamination on the
complementary strand).  The method uses pysam's get_aligned_pairs(with_seq=True)
which reconstructs reference bases from the MD tag.

Key design decisions:
- Forward-strand reads contribute to 5' C→T arrays.
- Reverse-strand reads contribute to 3' G→A arrays (representing the 5' end of
  the complementary strand as sequenced in reverse complement).
- Positional arrays are indexed by enumeration order of aligned pairs, NOT by
  raw query_pos, so soft-clipped bases do not distort the terminal indices.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pysam

from src.config import MIN_MD_TAG_FRACTION

logger = logging.getLogger(__name__)


@dataclass
class ReadFeatures:
    """Per-read features extracted during BAM traversal."""

    read_id: str
    read_length: int
    template_length: int            # abs(tlen) for paired-end, query_length for SE
    ct_terminal_count: int          # C→T mismatches in first n_terminal aligned positions
    ga_terminal_count: int          # G→A mismatches in last n_terminal aligned positions
    ct_terminal_opportunities: int  # reference C bases in the 5' terminal window
    ga_terminal_opportunities: int  # reference G bases in the 3' terminal window
    is_paired: bool
    is_reverse: bool


@dataclass
class DamageProfile:
    """Aggregated damage statistics for the whole library."""

    n_terminal: int
    n_reads_total: int
    n_reads_passed: int
    n_reads_no_md: int              # reads skipped due to missing MD tags

    # Shape: (n_terminal,) — per-position counts
    ct_count: np.ndarray            # C→T substitutions at 5' position i
    c_count: np.ndarray             # reference C bases at 5' position i (denominator)
    ga_count: np.ndarray            # G→A substitutions at 3' position i (0 = most terminal)
    g_count: np.ndarray             # reference G bases at 3' position i

    read_features: list[ReadFeatures] = field(default_factory=list)

    @property
    def ct_rate(self) -> np.ndarray:
        """C→T substitution rate per 5' position. Returns 0.0 where c_count == 0."""
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(self.c_count > 0, self.ct_count / self.c_count, 0.0)

    @property
    def ga_rate(self) -> np.ndarray:
        """G→A substitution rate per 3' position (index 0 = most terminal)."""
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(self.g_count > 0, self.ga_count / self.g_count, 0.0)


def profile_damage(
    bam_path: Path,
    min_mapq: int,
    min_length: int,
    n_terminal: int,
) -> DamageProfile:
    """
    Traverse a BAM file and compute per-position damage frequencies.

    Args:
        bam_path:    Path to a coordinate-sorted, indexed BAM file.
        min_mapq:    Minimum mapping quality; reads below this are skipped.
        min_length:  Minimum read length in bp; shorter reads are skipped.
        n_terminal:  Number of terminal positions to profile at each end.

    Returns:
        A DamageProfile containing positional counts and per-read features.

    Raises:
        FileNotFoundError: If bam_path does not exist.
        ValueError:        If the BAM has no index or too few MD-tagged reads.
    """
    bam_path = Path(bam_path)
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    ct_count = np.zeros(n_terminal, dtype=np.float64)
    c_count  = np.zeros(n_terminal, dtype=np.float64)
    ga_count = np.zeros(n_terminal, dtype=np.float64)
    g_count  = np.zeros(n_terminal, dtype=np.float64)

    read_features: list[ReadFeatures] = []
    n_reads_total = 0
    n_reads_passed = 0
    n_reads_no_md = 0
    n_reads_no_seq = 0

    try:
        bam = pysam.AlignmentFile(str(bam_path), "rb")
    except ValueError as exc:
        raise ValueError(f"Cannot open BAM file: {exc}") from exc

    with bam:
        # Verify index is present
        try:
            bam.check_index()
        except (ValueError, AttributeError):
            raise ValueError(
                f"BAM index not found for {bam_path}. "
                "Run 'samtools index <bam>' before classifying."
            )

        for read in bam.fetch():
            n_reads_total += 1

            # ── primary quality filters ──────────────────────────────────────
            if read.is_unmapped:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if read.query_sequence is None:
                n_reads_no_seq += 1
                continue
            if read.query_length < min_length:
                continue

            rf = _process_read(
                read, n_terminal, ct_count, c_count, ga_count, g_count
            )
            if rf is None:
                n_reads_no_md += 1
                continue

            n_reads_passed += 1
            read_features.append(rf)

    if n_reads_no_seq > 0:
        logger.warning(
            "%d reads had no stored query sequence and were skipped.", n_reads_no_seq
        )

    # ── MD tag coverage check ────────────────────────────────────────────────
    total_attempted = n_reads_passed + n_reads_no_md
    if total_attempted > 0:
        md_fraction = n_reads_passed / total_attempted
        if md_fraction < MIN_MD_TAG_FRACTION:
            raise ValueError(
                f"Only {md_fraction:.1%} of reads have valid MD tags "
                f"(minimum required: {MIN_MD_TAG_FRACTION:.0%}). "
                "Re-generate the BAM with 'samtools calmd' to add MD tags."
            )
        if n_reads_no_md > 0:
            logger.warning(
                "%d reads lacked MD tags and were skipped (%.1f%% of attempted reads).",
                n_reads_no_md,
                (1 - md_fraction) * 100,
            )

    if n_reads_passed == 0:
        logger.warning(
            "No reads passed all filters (total=%d). "
            "Consider relaxing --min-mapq or --min-length.",
            n_reads_total,
        )

    return DamageProfile(
        n_terminal=n_terminal,
        n_reads_total=n_reads_total,
        n_reads_passed=n_reads_passed,
        n_reads_no_md=n_reads_no_md,
        ct_count=ct_count,
        c_count=c_count,
        ga_count=ga_count,
        g_count=g_count,
        read_features=read_features,
    )


def _process_read(
    read: pysam.AlignedSegment,
    n_terminal: int,
    ct_count: np.ndarray,
    c_count: np.ndarray,
    ga_count: np.ndarray,
    g_count: np.ndarray,
) -> ReadFeatures | None:
    """
    Extract damage information from a single aligned read and accumulate counts.

    Uses get_aligned_pairs(with_seq=True) to obtain (query_pos, ref_pos, ref_base)
    triples.  ref_base is lowercased by pysam; we call .upper() before comparison.
    Returns None if the read has no usable aligned pairs with reference sequence
    (indicates missing MD tags).

    Positional arrays are accumulated by enumeration index of aligned pairs,
    NOT by raw query_pos, so soft-clipped bases do not shift terminal indices.

    For forward-strand reads: examine first n_terminal pairs for C→T (5' damage).
    For reverse-strand reads: examine last n_terminal pairs reversed for G→A
        (3' end of molecule = 5' end of complementary strand).
    """
    query_seq = read.query_sequence

    # Collect all aligned pairs with ref sequence, skipping None entries
    # (None query_pos = deletion; None ref_pos = insertion; None ref_base = no MD)
    try:
        raw_pairs = read.get_aligned_pairs(with_seq=True)
    except Exception:
        return None

    aligned_pairs = [
        (qp, rp, rb)
        for qp, rp, rb in raw_pairs
        if qp is not None and rp is not None and rb is not None
    ]

    if not aligned_pairs:
        return None

    # Check if MD tags were present (ref_base would all be None if absent)
    # Already filtered above: if aligned_pairs is non-empty, we have ref bases.

    ct_terminal_count = 0
    ct_terminal_opportunities = 0
    ga_terminal_count = 0
    ga_terminal_opportunities = 0

    if not read.is_reverse:
        # ── forward strand: 5' C→T damage ───────────────────────────────────
        five_prime_pairs = aligned_pairs[:n_terminal]
        for i, (qpos, _rpos, rbase) in enumerate(five_prime_pairs):
            rb = rbase.upper()
            qb = query_seq[qpos].upper()
            if rb == "C":
                c_count[i] += 1
                ct_terminal_opportunities += 1
                if qb == "T":
                    ct_count[i] += 1
                    ct_terminal_count += 1
    else:
        # ── reverse strand: 3' G→A damage ───────────────────────────────────
        # The 3' end of the molecule corresponds to the end of aligned pairs
        # for reverse-strand reads.  Reverse so index 0 = most terminal.
        three_prime_pairs = list(reversed(aligned_pairs[-n_terminal:]))
        for i, (qpos, _rpos, rbase) in enumerate(three_prime_pairs):
            rb = rbase.upper()
            qb = query_seq[qpos].upper()
            if rb == "G":
                g_count[i] += 1
                ga_terminal_opportunities += 1
                if qb == "A":
                    ga_count[i] += 1
                    ga_terminal_count += 1

    # ── template length (fragment size indicator) ────────────────────────────
    tlen = abs(read.template_length)
    template_length = tlen if tlen > 0 else read.query_length

    return ReadFeatures(
        read_id=read.query_name or "",
        read_length=read.query_length,
        template_length=template_length,
        ct_terminal_count=ct_terminal_count,
        ga_terminal_count=ga_terminal_count,
        ct_terminal_opportunities=ct_terminal_opportunities,
        ga_terminal_opportunities=ga_terminal_opportunities,
        is_paired=read.is_paired,
        is_reverse=read.is_reverse,
    )
