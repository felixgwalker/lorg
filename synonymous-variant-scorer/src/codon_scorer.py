"""Codon usage bias scoring against human codon frequency table."""

from __future__ import annotations

import math

from src.codon_tables import HUMAN_CODON_FREQ
from src.vcf_parser import SynonymousVariant


def score_codon_usage(variant: SynonymousVariant) -> float:
    """Return codon usage bias impact score in [0, 1].

    Higher score = alt codon is rarer relative to ref codon.
    """
    ref_freq = HUMAN_CODON_FREQ.get(variant.ref_codon.upper(), 0.5)
    alt_freq = HUMAN_CODON_FREQ.get(variant.alt_codon.upper(), 0.5)

    if ref_freq <= 0 or alt_freq <= 0:
        return 0.0

    delta = ref_freq - alt_freq

    if delta <= 0:
        return 0.0

    log_ratio = math.log2(ref_freq / alt_freq)
    score = min(1.0, log_ratio / 4.0)
    return max(0.0, score)
