"""Cotranslational folding scorer: CAI-proxy for co-translational folding risk.

Synonymous variants that introduce rarer codons can disrupt ribosomal pausing
patterns that are thought to be conserved to allow correct co-translational
folding of nascent peptides.  This scorer estimates the risk using:

  1. **CAI proxy** — the Codon Adaptation Index for the alt codon relative to
     the optimal codon for the same amino acid in human cells.  A low CAI
     means the alt codon is suboptimal, predicting ribosomal slowing.
  2. **Ramp region weighting** — variants in the first ~50 codons of the CDS
     (the 5' translation ramp) carry higher risk because this region is
     thought to be optimised for slow decoding that prevents ribosome queuing.
  3. **Rare-codon introduction** — introduction of a codon below the rare-
     codon frequency threshold adds an additional risk component.

Public API
----------
``score_folding(variant)``
    Pipeline mode; accepts ``SynonymousVariant``.  Returns [0, 1].

``score_protein_folding(codon_ref, codon_alt)``
    Direct mode; returns fold_risk_score in [0, 1].
"""

from __future__ import annotations

from typing import Union

from src.codon_tables import (
    HUMAN_CODON_FREQ,
    RARE_CODON_THRESHOLD,
    RARE_CODONS,
    STANDARD_GENETIC_CODE,
)
from src.vcf_parser import SynonymousVariant

# ---------------------------------------------------------------------------
# Optimal codons per amino acid (human): the codon with the highest RSCU in
# HUMAN_CODON_FREQ for each amino acid.
# ---------------------------------------------------------------------------
def _build_optimal_codons() -> dict[str, str]:
    """Return mapping amino_acid -> most frequent codon (human)."""
    aa_best: dict[str, tuple[str, float]] = {}
    for codon, aa in STANDARD_GENETIC_CODE.items():
        freq = HUMAN_CODON_FREQ.get(codon, 0.0)
        if aa not in aa_best or freq > aa_best[aa][1]:
            aa_best[aa] = (codon, freq)
    return {aa: codon for aa, (codon, _) in aa_best.items()}


_OPTIMAL_CODON: dict[str, str] = _build_optimal_codons()

# Ramp window: first N codons of the CDS
_RAMP_WINDOW = 50


def _cai_codon(codon: str) -> float:
    """Return the CAI value for a single codon.

    CAI(codon) = freq(codon) / max_freq(synonymous_family).
    For single-codon amino acids (Met, Trp) or stops, CAI = 1.0.

    Parameters
    ----------
    codon:
        Uppercase 3-letter DNA codon.

    Returns
    -------
    float
        CAI value in [0, 1].
    """
    codon = codon.upper()
    aa = STANDARD_GENETIC_CODE.get(codon)
    if aa is None:
        return 0.5  # unknown codon

    optimal = _OPTIMAL_CODON.get(aa, codon)
    max_freq = HUMAN_CODON_FREQ.get(optimal, 1.0)
    codon_freq = HUMAN_CODON_FREQ.get(codon, 0.0)

    if max_freq <= 0:
        return 1.0

    return min(1.0, codon_freq / max_freq)


def score_protein_folding(codon_ref: str, codon_alt: str) -> float:
    """Return fold_risk_score in [0, 1] based on CAI proxy for co-translational folding.

    A higher score indicates greater risk of aberrant co-translational folding
    due to the synonymous substitution.

    The score combines:
      - CAI drop: (CAI_ref - CAI_alt) normalised to [0, 1]
      - Rare codon introduction bonus (if alt codon crosses rare threshold)

    Parameters
    ----------
    codon_ref:
        Reference codon (uppercase 3-letter DNA).
    codon_alt:
        Alternate codon (uppercase 3-letter DNA).

    Returns
    -------
    float
        fold_risk_score in [0, 1].
    """
    cai_ref = _cai_codon(codon_ref)
    cai_alt = _cai_codon(codon_alt)

    # CAI drop: only penalise when alt is worse than ref
    cai_drop = max(0.0, cai_ref - cai_alt)

    # Rare codon component: extra risk if alt codon is below rare threshold
    ref_rare = codon_ref.upper() in RARE_CODONS
    alt_rare = codon_alt.upper() in RARE_CODONS
    rare_bonus = 0.0
    if alt_rare and not ref_rare:
        alt_freq = HUMAN_CODON_FREQ.get(codon_alt.upper(), 0.0)
        # How far below threshold is the alt codon?
        rare_bonus = min(0.5, (RARE_CODON_THRESHOLD - alt_freq) / RARE_CODON_THRESHOLD)

    fold_risk = min(1.0, cai_drop + rare_bonus)
    return max(0.0, fold_risk)


def _ramp_disruption(variant: SynonymousVariant) -> float:
    """Assess disruption of the 5' translation ramp region."""
    codon_position = variant.exon_pos // 3
    if codon_position > _RAMP_WINDOW:
        return 0.0
    cai_ref = _cai_codon(variant.ref_codon)
    cai_alt = _cai_codon(variant.alt_codon)
    cai_drop = max(0.0, cai_ref - cai_alt)
    ramp_weight = 1.0 - codon_position / _RAMP_WINDOW
    return min(1.0, cai_drop * ramp_weight * 2.0)


def score_folding(
    variant_or_ref: Union[SynonymousVariant, str],
    codon_alt: str | None = None,
) -> float:
    """Return cotranslational folding disruption score in [0, 1].

    **Pipeline mode** (called by ``pipeline.py``)::

        score_folding(variant)  →  float in [0, 1]

    **Direct mode**::

        score_folding(codon_ref, codon_alt)  →  fold_risk_score in [0, 1]

    In pipeline mode the score combines:
      - CAI-proxy risk (score_protein_folding)
      - Ramp region weighting (5' CDS position sensitivity)

    Parameters
    ----------
    variant_or_ref:
        ``SynonymousVariant`` (pipeline mode) or reference codon string.
    codon_alt:
        Required in direct mode.

    Returns
    -------
    float
        Folding risk score in [0, 1].
    """
    if isinstance(variant_or_ref, SynonymousVariant):
        variant = variant_or_ref
        cai_score = score_protein_folding(variant.ref_codon, variant.alt_codon)
        ramp_score = _ramp_disruption(variant)
        return min(1.0, max(0.0, 0.6 * cai_score + 0.4 * ramp_score))
    else:
        if codon_alt is None:
            raise TypeError(
                "score_folding() requires codon_alt when called with a string ref codon."
            )
        return score_protein_folding(variant_or_ref, codon_alt)
