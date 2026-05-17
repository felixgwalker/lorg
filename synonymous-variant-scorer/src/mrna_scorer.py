"""mRNA stability scoring: GC content change + CpG motif changes."""

from __future__ import annotations

from src.vcf_parser import SynonymousVariant


def score_mrna_stability(variant: SynonymousVariant) -> float:
    """Return mRNA stability impact score in [0, 1].

    Combines GC-content delta and CpG motif gain/loss contributions.
    """
    gc_score = _gc_delta_score(variant.ref_codon, variant.alt_codon)
    cpg_score = _cpg_score(variant.ref_codon, variant.alt_codon)
    return min(1.0, 0.6 * gc_score + 0.4 * cpg_score)


def _gc_content(codon: str) -> float:
    if not codon:
        return 0.0
    gc = sum(1 for b in codon.upper() if b in "GC")
    return gc / len(codon)


def _gc_delta_score(ref: str, alt: str) -> float:
    delta = abs(_gc_content(ref) - _gc_content(alt))
    return min(1.0, delta / (1.0 / 3))


def _cpg_score(ref: str, alt: str) -> float:
    ref_cpg = _count_cpg(ref)
    alt_cpg = _count_cpg(alt)
    delta = abs(ref_cpg - alt_cpg)
    return min(1.0, delta * 0.5)


def _count_cpg(seq: str) -> int:
    s = seq.upper()
    count = 0
    for i in range(len(s) - 1):
        if s[i] == "C" and s[i + 1] == "G":
            count += 1
    return count
