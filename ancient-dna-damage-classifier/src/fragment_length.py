"""
Fragment length distribution analysis.

For paired-end BAM data, template_length (TLEN field, absolute value) is used
as the fragment size estimate.  For single-end or unpaired reads, query_length
is used instead (template_length is 0 in those cases).

These values were extracted during BAM traversal by damage_profiler and stored
in ReadFeatures.template_length, so this module works purely from in-memory
data with no additional file I/O.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

from src.damage_profiler import ReadFeatures

logger = logging.getLogger(__name__)


@dataclass
class FragmentLengthStats:
    """Summary statistics of fragment/read length distribution."""

    mean: float
    median: float
    std: float
    min_len: int
    max_len: int
    n_reads: int
    is_paired: bool         # True if majority of reads had TLEN > 0
    hist_counts: np.ndarray
    hist_edges: np.ndarray  # length n_bins + 1


def compute_fragment_lengths(read_features: list[ReadFeatures]) -> FragmentLengthStats:
    """
    Compute summary statistics of the fragment/template length distribution.

    Args:
        read_features: List of per-read feature dicts from damage_profiler.

    Returns:
        A FragmentLengthStats dataclass.

    Raises:
        ValueError: If read_features is empty.
    """
    if not read_features:
        raise ValueError(
            "No read features available for fragment length analysis. "
            "The BAM file may contain no reads passing quality filters."
        )

    lengths = np.array([rf.template_length for rf in read_features], dtype=np.int64)

    # Determine whether the majority of reads are from paired-end sequencing
    n_paired = sum(1 for rf in read_features if rf.is_paired)
    is_paired = n_paired > len(read_features) / 2

    # Compute histogram with sensible bin width
    min_len = int(lengths.min())
    max_len = int(lengths.max())

    if max_len > min_len:
        # ~50 bins or 1 bp per bin for very short ranges
        n_bins = min(50, max_len - min_len)
        hist_counts, hist_edges = np.histogram(lengths, bins=n_bins)
    else:
        # All reads the same length — degenerate histogram
        hist_counts = np.array([len(lengths)], dtype=np.float64)
        hist_edges = np.array([min_len - 0.5, max_len + 0.5])

    return FragmentLengthStats(
        mean=float(np.mean(lengths)),
        median=float(np.median(lengths)),
        std=float(np.std(lengths)),
        min_len=min_len,
        max_len=max_len,
        n_reads=len(lengths),
        is_paired=is_paired,
        hist_counts=hist_counts.astype(np.float64),
        hist_edges=hist_edges,
    )
