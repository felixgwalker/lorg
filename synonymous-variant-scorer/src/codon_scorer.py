"""Codon usage bias scoring using organism-specific tAI proxy tables.

The primary public interface is ``score_codon_usage``, which accepts either:
  - a ``SynonymousVariant`` object (used by the pipeline), or
  - explicit ``codon_ref`` / ``codon_alt`` / ``organism`` strings (used directly).

Returned score in [0, 1]: higher = alt codon is rarer relative to ref codon
(greater potential for translational slowing / yield reduction).

The ``delta_tAI`` value (alt_freq - ref_freq, i.e. negative = worse codon) is
also available via ``score_codon_usage_delta``.
"""

from __future__ import annotations

import math
from typing import Union

from src.codon_tables import (
    HUMAN_CODON_FREQ,
    get_codon_usage_human,
    get_codon_usage_ecoli,
)
from src.vcf_parser import SynonymousVariant

# Supported organism identifiers
_ORGANISM_TABLES: dict[str, dict[str, float]] = {}  # populated lazily


def _get_table(organism: str) -> dict[str, float]:
    """Return the codon frequency table for *organism* (case-insensitive)."""
    key = organism.lower().strip()
    if key not in _ORGANISM_TABLES:
        if key in ("human", "homo_sapiens", "h_sapiens", "hs"):
            _ORGANISM_TABLES[key] = get_codon_usage_human()
        elif key in ("ecoli", "e_coli", "escherichia_coli", "ec"):
            _ORGANISM_TABLES[key] = get_codon_usage_ecoli()
        else:
            raise ValueError(
                f"Unknown organism '{organism}'. "
                "Supported values: 'human', 'ecoli'."
            )
    return _ORGANISM_TABLES[key]


def score_codon_usage_delta(
    codon_ref: str,
    codon_alt: str,
    organism: str = "human",
) -> float:
    """Return delta_tAI = tAI(alt) - tAI(ref).

    A negative value means the alt codon is rarer than the reference codon
    (potentially harmful for translational efficiency).

    Parameters
    ----------
    codon_ref:
        Reference codon (uppercase 3-letter DNA).
    codon_alt:
        Alternate codon (uppercase 3-letter DNA).
    organism:
        ``'human'`` (default) or ``'ecoli'``.

    Returns
    -------
    float
        delta_tAI in [-1, 1].
    """
    table = _get_table(organism)
    ref_freq = table.get(codon_ref.upper(), 0.5)
    alt_freq = table.get(codon_alt.upper(), 0.5)
    return alt_freq - ref_freq


def score_codon_usage(
    variant_or_ref: Union[SynonymousVariant, str],
    codon_alt: str | None = None,
    organism: str = "human",
) -> float:
    """Return codon usage bias impact score in [0, 1].

    Can be called in two ways:

    **Pipeline mode** (called by ``pipeline.py``)::

        score_codon_usage(variant)   # SynonymousVariant, human table

    **Direct mode** (called externally with explicit codons)::

        score_codon_usage(codon_ref, codon_alt, organism)

    In both cases a higher score indicates the alt codon is rarer than the
    ref codon — i.e. greater potential for ribosomal slowing / yield drop.
    Score of 0.0 means alt codon is equally or more frequent than ref.

    Parameters
    ----------
    variant_or_ref:
        Either a ``SynonymousVariant`` (pipeline mode) or a reference codon
        string (direct mode).
    codon_alt:
        Required in direct mode; the alternate codon.
    organism:
        ``'human'`` (default) or ``'ecoli'``.  Ignored in pipeline mode
        (always uses human table).

    Returns
    -------
    float
        Impact score in [0, 1].
    """
    if isinstance(variant_or_ref, SynonymousVariant):
        variant = variant_or_ref
        ref_codon = variant.ref_codon.upper()
        alt_codon = variant.alt_codon.upper()
        table = HUMAN_CODON_FREQ
    else:
        if codon_alt is None:
            raise TypeError(
                "score_codon_usage() requires codon_alt when called with a string ref codon."
            )
        ref_codon = variant_or_ref.upper()
        alt_codon = codon_alt.upper()
        table = _get_table(organism)

    ref_freq = table.get(ref_codon, 0.5)
    alt_freq = table.get(alt_codon, 0.5)

    if ref_freq <= 0 or alt_freq <= 0:
        return 0.0

    # Only score as impactful when the alt codon is rarer than the ref codon.
    delta = ref_freq - alt_freq
    if delta <= 0:
        return 0.0

    # Log2 ratio, capped at 4 bits (= 16-fold difference) → score in [0, 1].
    log_ratio = math.log2(ref_freq / alt_freq)
    score = min(1.0, log_ratio / 4.0)
    return max(0.0, score)
