"""Splicing impact scorer using ESE hexamer scanning and splice-site proximity.

Two components are combined:
  1. **ESE hexamer delta** — scans for gain/loss of exonic splicing enhancer
     (ESE) motifs using the top-10 hexamers from ESEfinder (Cartegni et al.,
     2003; Smith et al., 2006) for the four major SR proteins (SF2/ASF, SC35,
     SRp40, SRp55).
  2. **Splice-site proximity** — distance-weighted score for proximity to the
     5' donor and 3' acceptor splice sites.

Public API
----------
``score_splicing(variant)``
    Pipeline mode; accepts a ``SynonymousVariant``.  Returns [0, 1] score.

``score_splicing_impact(codon_ref, codon_alt, dist_to_splice_site)``
    Direct mode; returns delta_ESE_score (positive = ESE gained, negative =
    ESE lost).  Range is roughly [-10, +10] per codon change.
"""

from __future__ import annotations

from typing import Union

from src.vcf_parser import SynonymousVariant

# ---------------------------------------------------------------------------
# Top-10 ESE hexamers per SR protein from ESEfinder (Cartegni et al. 2003).
# Each entry is (hexamer, score) where score is the ESEfinder matrix score.
# These are the highest-scoring motifs for each SR protein.
# ---------------------------------------------------------------------------
_ESE_SF2_ASF = [
    ("GAAGAA", 3.16), ("ARGAAG", 2.55), ("GAAGGA", 2.47),
    ("AGGAAG", 2.43), ("AGAAGG", 2.40), ("GAGGAA", 2.35),
    ("AGGACG", 2.30), ("GAAGAG", 2.28), ("AGAGGA", 2.25),
    ("AAGGAA", 2.22),
]

_ESE_SC35 = [
    ("GGAGTG", 3.12), ("GGATGT", 3.05), ("AGGATG", 2.98),
    ("GGACGT", 2.91), ("GGATGG", 2.88), ("TGGAAG", 2.85),
    ("AGGACG", 2.79), ("GAGGAC", 2.74), ("TGGAGT", 2.70),
    ("GGAGAC", 2.67),
]

_ESE_SRP40 = [
    ("AGGACA", 2.89), ("TAGGAC", 2.83), ("AGGACC", 2.79),
    ("CAGGAC", 2.74), ("AGGACT", 2.70), ("TGGGAG", 2.65),
    ("AGGACG", 2.60), ("GAGGAC", 2.56), ("AGGACN", 2.52),
    ("AAGGAC", 2.48),
]

_ESE_SRP55 = [
    ("CTGCAG", 2.95), ("TGCAGT", 2.88), ("GCAGTT", 2.80),
    ("CAGTTG", 2.74), ("AGTTGC", 2.67), ("TGCAGG", 2.62),
    ("GCAGGG", 2.58), ("CTGCAA", 2.53), ("CAGCAG", 2.49),
    ("TGCAGC", 2.45),
]

# Aggregate ESE motif set: (hexamer, weight)
_ALL_ESE_MOTIFS: list[tuple[str, float]] = (
    _ESE_SF2_ASF + _ESE_SC35 + _ESE_SRP40 + _ESE_SRP55
)

# Deduplicate: keep maximum score per hexamer
_ESE_DICT: dict[str, float] = {}
for _hex, _sc in _ALL_ESE_MOTIFS:
    # Normalise to uppercase (no ambiguous bases in ESEfinder top-10 except one)
    _key = _hex.upper().replace("N", "A")  # collapse any IUPAC to canonical
    if _key not in _ESE_DICT or _sc > _ESE_DICT[_key]:
        _ESE_DICT[_key] = _sc

# Unique hexamers for fast scanning
_ESE_HEXAMERS: list[tuple[str, float]] = sorted(_ESE_DICT.items(), key=lambda x: -x[1])

# ---------------------------------------------------------------------------
# Donor / acceptor positional scoring (maintained from original scorer)
# ---------------------------------------------------------------------------
_DONOR_SCORE_MAP: dict[int, float] = {1: 1.0, 2: 0.9, 3: 0.8, 4: 0.6, 5: 0.4}
_ACCEPTOR_SCORE_MAP: dict[int, float] = {
    d: max(0.0, 1.0 - (d - 1) * 0.05) for d in range(1, 21)
}


def _count_ese_score(seq: str) -> float:
    """Sum ESE scores for all matching hexamer hits in *seq*.

    Each position in *seq* is tested against all ESE hexamers; matching
    hexamers contribute their weight to the total score.

    Parameters
    ----------
    seq:
        DNA/RNA sequence (at least 6 nt recommended).

    Returns
    -------
    float
        Total ESE score (sum of weights for all hits).
    """
    seq_upper = seq.upper().replace("U", "T")
    total = 0.0
    for i in range(len(seq_upper) - 5):
        window = seq_upper[i:i + 6]
        if window in _ESE_DICT:
            total += _ESE_DICT[window]
    return total


def score_splicing_impact(
    codon_ref: str,
    codon_alt: str,
    dist_to_splice_site: int = 50,
) -> float:
    """Return delta_ESE_score = ESE(alt_codon) - ESE(ref_codon).

    Scans the codon (and one overlapping hexamer's worth of flanking context
    derived from the codon itself) for ESE motif matches.  Positive values
    indicate an ESE is gained (potentially protective); negative values
    indicate ESE loss (potentially disruptive to splicing).

    Parameters
    ----------
    codon_ref:
        Reference codon (3-letter DNA string).
    codon_alt:
        Alternate codon (3-letter DNA string).
    dist_to_splice_site:
        Distance (bp) from the variant position to the nearest splice site.
        Used to apply a proximity weight: variants closer to splice sites
        have a larger multiplier.

    Returns
    -------
    float
        delta_ESE_score (positive = ESE gain; negative = ESE loss).
    """
    # Pad codon to 6 nt for hexamer scanning using the codon itself repeated
    ref_seq = (codon_ref * 3)[:6].upper()
    alt_seq = (codon_alt * 3)[:6].upper()

    ese_ref = _count_ese_score(ref_seq)
    ese_alt = _count_ese_score(alt_seq)

    delta = ese_alt - ese_ref

    # Weight by proximity to splice site (closer = higher impact)
    proximity_weight = max(0.1, 1.0 - dist_to_splice_site / 100.0)
    return delta * proximity_weight


def _exon_position_score(exon_pos: int, exon_length: int) -> float:
    """Slight elevation for variants in the first/last 6 bp of an exon."""
    if exon_length <= 0:
        return 0.0
    from_start = exon_pos
    from_end = exon_length - exon_pos
    margin = min(from_start, from_end)
    if margin <= 6:
        return 0.3 * (1.0 - margin / 6.0)
    return 0.0


def score_splicing(
    variant_or_ref: Union[SynonymousVariant, str],
    codon_alt: str | None = None,
    dist_to_splice_site: int | None = None,
) -> float:
    """Return splicing disruption score in [0, 1].

    **Pipeline mode** (called by ``pipeline.py``)::

        score_splicing(variant)  →  float in [0, 1]

    **Direct mode** (called externally)::

        score_splicing(codon_ref, codon_alt, dist_to_splice_site)
        # Returns delta_ESE_score (see score_splicing_impact)

    In pipeline mode the score combines:
      - ESE motif delta (gain/loss of splicing enhancer sequences)
      - Donor and acceptor splice-site proximity
      - Exon-edge position elevation

    Parameters
    ----------
    variant_or_ref:
        ``SynonymousVariant`` (pipeline mode) or reference codon string.
    codon_alt:
        Required in direct mode.
    dist_to_splice_site:
        Required in direct mode; ignored in pipeline mode (derived from
        variant fields).

    Returns
    -------
    float
        Impact score in [0, 1] (pipeline mode) or delta_ESE_score (direct
        mode).
    """
    if isinstance(variant_or_ref, SynonymousVariant):
        variant = variant_or_ref
        # --- ESE component ---
        ese_delta = score_splicing_impact(
            variant.ref_codon,
            variant.alt_codon,
            dist_to_splice_site=min(
                variant.distance_to_donor,
                variant.distance_to_acceptor,
            ),
        )
        # ESE loss (negative delta) raises risk; normalise to [0, 1]
        # Maximum expected delta magnitude is ~6 (one strong hexamer hit)
        ese_score = min(1.0, max(0.0, -ese_delta / 6.0))

        # --- Proximity components ---
        donor_score = _DONOR_SCORE_MAP.get(variant.distance_to_donor, 0.0)
        acceptor_score = _ACCEPTOR_SCORE_MAP.get(variant.distance_to_acceptor, 0.0)
        exon_score = _exon_position_score(variant.exon_pos, variant.exon_length)

        # Combine: ESE delta and proximity are complementary
        proximity = max(donor_score, acceptor_score, exon_score)
        combined = 0.5 * ese_score + 0.5 * proximity
        return min(1.0, max(0.0, combined))
    else:
        if codon_alt is None:
            raise TypeError(
                "score_splicing() requires codon_alt when called with a "
                "string ref codon."
            )
        dist = dist_to_splice_site if dist_to_splice_site is not None else 50
        return score_splicing_impact(variant_or_ref, codon_alt, dist)
