"""mRNA stability scoring using a simplified nearest-neighbour RNA folding model.

The scorer estimates the change in minimum free energy (delta_MFE) for a
40-nucleotide window centred on the variant site, using standard nearest-
neighbour thermodynamic parameters for RNA (SantaLucia & Turner, 2004).

No external dependencies are required — the nearest-neighbour stacking
parameters are hardcoded from published tables.

Public API
----------
``score_mrna_stability(variant)``
    Called by the pipeline; accepts a ``SynonymousVariant``.  Returns a
    normalised impact score in [0, 1].

``score_mrna_stability(codon_ref, codon_alt, context_seq)``
    Called directly with explicit codon strings.  Returns delta_MFE in kcal/mol.
    A negative delta_MFE means the alt codon makes the local structure more
    stable (lower energy); positive means less stable (destabilised).
"""

from __future__ import annotations

import math
from typing import Union

from src.vcf_parser import SynonymousVariant

# ---------------------------------------------------------------------------
# Nearest-neighbour RNA stacking free energies (ΔG°37, kcal/mol)
# Propagation parameters from Xia et al. (1998) / Mathews et al. (1999).
# Keys are 5'-XY-3' dinucleotide pairs written as "XY" for the 5'→3' strand.
# Values are ΔG°37 in kcal/mol.
# ---------------------------------------------------------------------------
_NN_DG: dict[str, float] = {
    "AA": -0.93, "AU": -1.10, "AC": -2.24, "AG": -2.08,
    "UA": -1.33, "UU": -0.93, "UC": -2.11, "UG": -2.35,
    "CA": -2.11, "CU": -1.36, "CC": -3.26, "CG": -3.42,
    "GA": -2.35, "GU": -2.24, "GC": -3.42, "GG": -3.26,
}

# Initiation penalty for an AU terminal pair (kcal/mol)
_INIT_AU = 0.45
# Initiation penalty for a GC terminal pair (kcal/mol)
_INIT_GC = 0.98

# Hairpin loop closing penalty (simplified; kcal/mol per unpaired base, capped)
_HAIRPIN_MIN_SIZE = 3
_HAIRPIN_PENALTY = 5.4  # kcal/mol for a 3-base loop

# Window size (nt) centred on the mutated codon
_WINDOW = 40


def _dna_to_rna(seq: str) -> str:
    """Convert DNA sequence to RNA (T -> U, uppercase)."""
    return seq.upper().replace("T", "U")


def _nn_energy(seq: str) -> float:
    """Compute simplified nearest-neighbour stacking energy for an RNA duplex.

    This is a linear approximation: sum the stacking parameters for all
    consecutive dinucleotide pairs in the sequence.  This does not fold the
    molecule; it estimates the energetic propensity for stable secondary
    structure formation (lower = more structured / stable).

    Parameters
    ----------
    seq:
        RNA sequence (U-based, uppercase).

    Returns
    -------
    float
        Summed nearest-neighbour ΔG in kcal/mol.
    """
    rna = _dna_to_rna(seq)
    energy = 0.0
    for i in range(len(rna) - 1):
        dinuc = rna[i:i + 2]
        energy += _NN_DG.get(dinuc, -1.5)  # fall back to average for non-ACGU

    # Add initiation term based on terminal nucleotides
    if rna:
        energy += _INIT_AU if rna[0] in "AU" else _INIT_GC
        energy += _INIT_AU if rna[-1] in "AU" else _INIT_GC

    return energy


def _mfe_estimate(seq: str) -> float:
    """Estimate MFE for a short RNA sequence.

    Uses the nearest-neighbour stacking energy divided by a normalisation
    factor to approximate the minimum free energy of the most stable hairpin
    or stem-loop that can form.  The estimate is qualitative but monotonically
    related to real MFE for short windows.

    Parameters
    ----------
    seq:
        Nucleotide sequence (DNA or RNA, will be converted).

    Returns
    -------
    float
        Estimated MFE in kcal/mol (negative = more stable).
    """
    rna = _dna_to_rna(seq)
    n = len(rna)
    if n < _HAIRPIN_MIN_SIZE + 2:
        return 0.0

    # Score the best candidate stem: try all hairpin lengths from 4 to n//2
    # For each stem length, compute the stacking energy of the stem portion.
    best_mfe = 0.0
    for stem_len in range(2, n // 2 + 1):
        loop_len = n - 2 * stem_len
        if loop_len < _HAIRPIN_MIN_SIZE:
            continue
        # 5' stem arm
        arm5 = rna[:stem_len]
        # 3' stem arm (reverse complement check skipped; use stacking of arm5)
        arm3 = rna[n - stem_len:]

        # Stacking energy of the 5' arm
        stem_energy = 0.0
        for i in range(stem_len - 1):
            dinuc5 = arm5[i:i + 2]
            stem_energy += _NN_DG.get(dinuc5, -1.5)

        # Penalise for loop entropy (simplified: +1.0 kcal/mol per base > 3)
        loop_penalty = _HAIRPIN_PENALTY + max(0, loop_len - 3) * 1.0

        candidate_mfe = stem_energy + loop_penalty
        if candidate_mfe < best_mfe:
            best_mfe = candidate_mfe

    return best_mfe


def _build_window(codon_ref: str, codon_alt: str, context_seq: str | None) -> tuple[str, str]:
    """Construct 40-nt window sequences for reference and alt.

    If *context_seq* is provided, embed the codon at the centre; otherwise
    pad with a neutral GC-balanced flanking sequence.

    Returns
    -------
    tuple[str, str]
        (ref_window, alt_window) each of length *_WINDOW*.
    """
    half = _WINDOW // 2
    flank_len = half - len(codon_ref) // 2

    if context_seq and len(context_seq) >= _WINDOW:
        # Find approximate codon centre in context
        centre = len(context_seq) // 2
        start = max(0, centre - half)
        end = start + _WINDOW
        if end > len(context_seq):
            end = len(context_seq)
            start = max(0, end - _WINDOW)
        ref_window = context_seq[start:end]
        # Replace codon in window
        codon_start = centre - start - len(codon_ref) // 2
        codon_start = max(0, min(codon_start, len(ref_window) - len(codon_ref)))
        alt_window = (
            ref_window[:codon_start]
            + codon_alt
            + ref_window[codon_start + len(codon_ref):]
        )
        # Trim/pad to _WINDOW
        ref_window = (ref_window + "G" * _WINDOW)[:_WINDOW]
        alt_window = (alt_window + "G" * _WINDOW)[:_WINDOW]
    else:
        # Neutral flanking: alternating GC to avoid extreme bias
        flank = ("GC" * (flank_len // 2 + 2))[:flank_len]
        ref_window = (flank + codon_ref + flank)[:_WINDOW]
        alt_window = (flank + codon_alt + flank)[:_WINDOW]
        # Pad to exactly _WINDOW
        ref_window = (ref_window + "G" * _WINDOW)[:_WINDOW]
        alt_window = (alt_window + "G" * _WINDOW)[:_WINDOW]

    return ref_window, alt_window


def score_mrna_stability_delta(
    codon_ref: str,
    codon_alt: str,
    context_seq: str | None = None,
) -> float:
    """Return delta_MFE = MFE(alt) - MFE(ref) in kcal/mol.

    A negative delta_MFE indicates the alt codon stabilises local mRNA
    secondary structure (more structured → typically slower translation /
    lower ribosome processivity in the 5' UTR or near splice sites).
    A positive value indicates destabilisation.

    Parameters
    ----------
    codon_ref:
        Reference codon (3-letter DNA/RNA string).
    codon_alt:
        Alternate codon (3-letter DNA/RNA string).
    context_seq:
        Optional 40+ nt genomic sequence centred on the variant.  If omitted,
        a synthetic GC-balanced flanking sequence is used.

    Returns
    -------
    float
        delta_MFE in kcal/mol.
    """
    ref_window, alt_window = _build_window(codon_ref, codon_alt, context_seq)
    mfe_ref = _mfe_estimate(ref_window)
    mfe_alt = _mfe_estimate(alt_window)
    return mfe_alt - mfe_ref


def score_mrna_stability(
    variant_or_ref: Union[SynonymousVariant, str],
    codon_alt: str | None = None,
    context_seq: str | None = None,
) -> float:
    """Return mRNA stability impact score.

    **Pipeline mode** (called by ``pipeline.py``)::

        score_mrna_stability(variant)  →  float in [0, 1]

    **Direct mode** (called externally)::

        score_mrna_stability(codon_ref, codon_alt, context_seq)  →  delta_MFE (kcal/mol)

    In pipeline mode the function returns a normalised [0, 1] impact score
    where higher = greater predicted destabilisation of local mRNA structure.

    Parameters
    ----------
    variant_or_ref:
        ``SynonymousVariant`` (pipeline) or reference codon string (direct).
    codon_alt:
        Required in direct mode.
    context_seq:
        Optional 40-nt window sequence (direct mode only; ignored in pipeline
        mode since the variant carries no genomic context string).

    Returns
    -------
    float
        Impact score in [0, 1] (pipeline mode) or delta_MFE in kcal/mol
        (direct mode).
    """
    if isinstance(variant_or_ref, SynonymousVariant):
        variant = variant_or_ref
        delta = score_mrna_stability_delta(variant.ref_codon, variant.alt_codon)
        # Normalise: ±5 kcal/mol range → [0, 1].
        # Destabilisation (positive delta) raises risk; stabilisation lowers it.
        # We treat extreme stabilisation at 5' end as also potentially disruptive
        # (pausing), so we take |delta| for the normalised score.
        normalised = min(1.0, abs(delta) / 5.0)
        return max(0.0, normalised)
    else:
        if codon_alt is None:
            raise TypeError(
                "score_mrna_stability() requires codon_alt when called with "
                "a string ref codon."
            )
        return score_mrna_stability_delta(variant_or_ref, codon_alt, context_seq)
