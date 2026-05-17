"""Cotranslational folding scorer: codon ramp disruption and rare codon clustering."""

from __future__ import annotations

from src.codon_tables import HUMAN_CODON_FREQ, RARE_CODONS
from src.vcf_parser import SynonymousVariant


_RAMP_WINDOW = 50


def score_folding(variant: SynonymousVariant) -> float:
    """Return cotranslational folding disruption score in [0, 1].

    Combines ramp region disruption (5' end of CDS) and rare codon
    introduction weighted by local clustering potential.
    """
    ramp_score = _ramp_disruption(variant)
    cluster_score = _rare_codon_introduced(variant)
    return min(1.0, 0.5 * ramp_score + 0.5 * cluster_score)


def _ramp_disruption(variant: SynonymousVariant) -> float:
    codon_position = variant.exon_pos // 3
    if codon_position > _RAMP_WINDOW:
        return 0.0
    ref_freq = HUMAN_CODON_FREQ.get(variant.ref_codon.upper(), 0.5)
    alt_freq = HUMAN_CODON_FREQ.get(variant.alt_codon.upper(), 0.5)
    freq_drop = max(0.0, ref_freq - alt_freq)
    ramp_weight = 1.0 - codon_position / _RAMP_WINDOW
    return min(1.0, freq_drop * ramp_weight * 3.0)


def _rare_codon_introduced(variant: SynonymousVariant) -> float:
    ref_rare = variant.ref_codon.upper() in RARE_CODONS
    alt_rare = variant.alt_codon.upper() in RARE_CODONS
    if alt_rare and not ref_rare:
        alt_freq = HUMAN_CODON_FREQ.get(variant.alt_codon.upper(), 0.1)
        return min(1.0, (RARE_CODON_THRESHOLD - alt_freq) / RARE_CODON_THRESHOLD + 0.5)
    return 0.0


RARE_CODON_THRESHOLD = 0.15
