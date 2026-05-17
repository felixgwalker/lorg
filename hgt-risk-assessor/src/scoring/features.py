"""
Feature extraction for the three-layer HGT risk model.

Each public function returns a FeatureResult.  Features backed by external
databases or optional inputs degrade gracefully: they return
FeatureResult(available=False, score=None) rather than raising or fabricating
data.  The caller (layers.py) re-normalises weights over available features
only.

Feature names are the canonical keys from config.LAYER_FEATURE_WEIGHTS.

Scientific caveats are documented inline.  Where an implementation is a first-
pass approximation, it is clearly labelled with a NOTE or TODO comment.
"""

from __future__ import annotations

import logging
import re
from typing import Optional

from src.config import LAYER_FEATURE_WEIGHTS
from src.models import FeatureResult, HostProfile, QuerySequence

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def _rc(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _w(layer: str, name: str) -> float:
    """Look up canonical within-layer weight."""
    return LAYER_FEATURE_WEIGHTS[layer][name]


# ---------------------------------------------------------------------------
# Codon usage reference tables
# ---------------------------------------------------------------------------
# Approximate per-1000-codon frequencies for common laboratory organisms.
# Source: Kazusa Codon Usage Database / Nakamura et al. 2000, Nucleic Acids Res.
# These are rounded approximations; replace with CoCoPUTs data for higher accuracy.
# Stop codons (TAA, TAG, TGA) are excluded from the comparison.
#
# NOTE: These values encode known strong biases (e.g., AGA/AGG rare in E. coli,
# CTG >> other Leu codons) and are scientifically reasonable for screening
# purposes.  They should NOT be treated as exact calibration data.

_CODON_TABLES: dict[str, dict[str, float]] = {
    "ecoli_k12": {
        # Phe
        "TTT": 23.0, "TTC": 17.0,
        # Leu
        "TTA": 14.0, "TTG": 13.0, "CTT": 11.0, "CTC": 10.0,
        "CTA":  4.0, "CTG": 52.0,
        # Ile
        "ATT": 28.0, "ATC": 27.0, "ATA":  5.0,
        # Met
        "ATG": 27.0,
        # Val
        "GTT": 19.0, "GTC": 16.0, "GTA": 11.0, "GTG": 27.0,
        # Ser
        "TCT":  8.0, "TCC":  8.0, "TCA":  7.0, "TCG":  9.0,
        "AGT":  9.0, "AGC": 16.0,
        # Pro
        "CCT":  7.0, "CCC":  5.0, "CCA":  9.0, "CCG": 25.0,
        # Thr
        "ACT": 10.0, "ACC": 23.0, "ACA":  7.0, "ACG": 13.0,
        # Ala
        "GCT": 15.0, "GCC": 22.0, "GCA": 17.0, "GCG": 33.0,
        # Tyr
        "TAT": 16.0, "TAC": 12.0,
        # Cys
        "TGT":  5.0, "TGC":  7.0,
        # Trp
        "TGG": 15.0,
        # His
        "CAT": 12.0, "CAC": 10.0,
        # Gln
        "CAA": 15.0, "CAG": 30.0,
        # Asn
        "AAT": 17.0, "AAC": 22.0,
        # Lys
        "AAA": 34.0, "AAG": 12.0,
        # Asp
        "GAT": 29.0, "GAC": 23.0,
        # Glu
        "GAA": 39.0, "GAG": 18.0,
        # Arg  — AGA/AGG are rare codons in E. coli
        "CGT": 21.0, "CGC": 20.0, "CGA":  4.0, "CGG":  5.0,
        "AGA":  2.0, "AGG":  2.0,
        # Gly
        "GGT": 24.0, "GGC": 29.0, "GGA":  9.0, "GGG":  7.0,
    },
    "bsubtilis_168": {
        # Phe
        "TTT": 22.0, "TTC": 16.0,
        # Leu
        "TTA": 19.0, "TTG": 17.0, "CTT": 14.0, "CTC":  9.0,
        "CTA":  7.0, "CTG": 21.0,
        # Ile
        "ATT": 30.0, "ATC": 18.0, "ATA": 10.0,
        # Met
        "ATG": 25.0,
        # Val
        "GTT": 22.0, "GTC": 11.0, "GTA": 13.0, "GTG": 21.0,
        # Ser
        "TCT": 11.0, "TCC":  9.0, "TCA": 11.0, "TCG":  7.0,
        "AGT": 12.0, "AGC": 11.0,
        # Pro
        "CCT": 11.0, "CCC":  7.0, "CCA": 14.0, "CCG": 11.0,
        # Thr
        "ACT": 15.0, "ACC": 16.0, "ACA": 14.0, "ACG": 11.0,
        # Ala
        "GCT": 22.0, "GCC": 16.0, "GCA": 21.0, "GCG": 16.0,
        # Tyr
        "TAT": 18.0, "TAC": 11.0,
        # Cys
        "TGT":  6.0, "TGC":  5.0,
        # Trp
        "TGG": 12.0,
        # His
        "CAT": 14.0, "CAC":  8.0,
        # Gln
        "CAA": 20.0, "CAG": 17.0,
        # Asn
        "AAT": 24.0, "AAC": 16.0,
        # Lys
        "AAA": 36.0, "AAG": 13.0,
        # Asp
        "GAT": 33.0, "GAC": 17.0,
        # Glu
        "GAA": 40.0, "GAG": 15.0,
        # Arg
        "CGT": 14.0, "CGC":  8.0, "CGA":  8.0, "CGG":  6.0,
        "AGA": 14.0, "AGG":  8.0,
        # Gly
        "GGT": 21.0, "GGC": 16.0, "GGA": 15.0, "GGG": 11.0,
    },
}

# Synonymous groups (stop codons excluded)
_SYN_GROUPS: list[list[str]] = [
    ["TTT", "TTC"], ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    ["ATT", "ATC", "ATA"], ["ATG"],
    ["GTT", "GTC", "GTA", "GTG"],
    ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    ["CCT", "CCC", "CCA", "CCG"],
    ["ACT", "ACC", "ACA", "ACG"],
    ["GCT", "GCC", "GCA", "GCG"],
    ["TAT", "TAC"], ["TGT", "TGC"], ["TGG"],
    ["CAT", "CAC"], ["CAA", "CAG"],
    ["AAT", "AAC"], ["AAA", "AAG"],
    ["GAT", "GAC"], ["GAA", "GAG"],
    ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    ["GGT", "GGC", "GGA", "GGG"],
]
_ALL_SENSE_CODONS: set[str] = {c for grp in _SYN_GROUPS for c in grp}


def _normalise_codon_table(raw: dict[str, float]) -> dict[str, float]:
    """Normalise a raw codon frequency table to sum to 1.0 over sense codons."""
    total = sum(v for k, v in raw.items() if k in _ALL_SENSE_CODONS)
    if total == 0:
        return {c: 1.0 / 61 for c in _ALL_SENSE_CODONS}
    return {c: raw.get(c, 0.0) / total for c in _ALL_SENSE_CODONS}


# Pre-normalised reference tables
_REFS: dict[str, dict[str, float]] = {
    name: _normalise_codon_table(table)
    for name, table in _CODON_TABLES.items()
}

# Equal-usage reference: uniform across all synonymous groups
_equal_table: dict[str, float] = {}
for _grp in _SYN_GROUPS:
    _p = 1.0 / len(_grp) / len(_SYN_GROUPS)
    for _c in _grp:
        _equal_table[_c] = _p
_REFS["equal"] = _normalise_codon_table(_equal_table)


def _select_reference(host_id: str) -> tuple[str, dict[str, float]]:
    """Choose the nearest reference codon table from the host identifier."""
    hid = host_id.lower().replace(" ", "").replace(".", "").replace("_", "")
    if any(x in hid for x in ("ecoli", "escherichiacoli", "coli")):
        return "ecoli_k12", _REFS["ecoli_k12"]
    if any(x in hid for x in ("bacillus", "bsubtilis", "subtilis")):
        return "bsubtilis_168", _REFS["bsubtilis_168"]
    return "equal (no host-specific reference available)", _REFS["equal"]


def _query_codon_freqs(sequence: str) -> dict[str, float]:
    """Count codons in all three forward reading frames, return normalised frequencies."""
    counts: dict[str, float] = {c: 0.0 for c in _ALL_SENSE_CODONS}
    seq = sequence.upper()
    total = 0
    for frame in range(3):
        for i in range(frame, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            if codon in counts:
                counts[codon] += 1
                total += 1
    if total == 0:
        return counts
    return {c: v / total for c, v in counts.items()}


# ---------------------------------------------------------------------------
# ORF utilities (shared by gene_completeness and payload_count)
# ---------------------------------------------------------------------------

_STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})


def _scan_orfs(sequence: str, min_nt: int = 90) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """
    Scan all 6 reading frames for ORFs.

    Returns (complete, partial) where each element is (start, end) in
    forward-strand coordinates (0-based, end exclusive).

    complete = ATG + in-frame stop, length >= min_nt
    partial  = ATG with no in-frame stop before sequence end (truncated),
               length >= min_nt
    """
    complete: list[tuple[int, int]] = []
    partial:  list[tuple[int, int]] = []
    seq_len = len(sequence)

    for strand_seq, forward in ((sequence.upper(), True), (_rc(sequence).upper(), False)):
        for frame in range(3):
            in_orf = False
            orf_start_in_strand = 0
            i = frame
            while i + 3 <= len(strand_seq):
                codon = strand_seq[i:i + 3]
                if not in_orf:
                    if codon == "ATG":
                        in_orf = True
                        orf_start_in_strand = i
                else:
                    if codon in _STOP_CODONS:
                        orf_end = i + 3
                        orf_len = orf_end - orf_start_in_strand
                        if orf_len >= min_nt:
                            if forward:
                                complete.append((orf_start_in_strand, orf_end))
                            else:
                                # Convert reverse-complement coordinates to forward
                                fwd_start = seq_len - orf_end
                                fwd_end   = seq_len - orf_start_in_strand
                                complete.append((fwd_start, fwd_end))
                        in_orf = False
                i += 3
            # Truncated ORF at end of sequence (no stop codon found)
            if in_orf:
                orf_len = i - orf_start_in_strand
                if orf_len >= min_nt:
                    if forward:
                        partial.append((orf_start_in_strand, i))
                    else:
                        fwd_start = seq_len - i
                        fwd_end   = seq_len - orf_start_in_strand
                        partial.append((fwd_start, fwd_end))

    return complete, partial


# ---------------------------------------------------------------------------
# Promoter motif utilities
# ---------------------------------------------------------------------------

def _hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def _count_promoter_like(sequence: str, max_mm_10: int = 2, max_mm_35: int = 2) -> int:
    """
    Count sigma-70-like prokaryotic promoter elements in the sequence.

    Looks for the canonical -10 (TATAAT) and -35 (TTGACA) boxes with up to
    max_mm mismatches each, with 15–19 bp spacer between them.

    NOTE: This is a heuristic proxy, NOT a definitive promoter predictor.
    Many false positives are expected in random sequence.  The feature is
    most useful as a relative comparator, not an absolute measure.
    """
    seq = sequence.upper()
    n = len(seq)
    matches = 0
    BOX_35 = "TTGACA"
    BOX_10 = "TATAAT"
    for i in range(n - 6):
        if _hamming(seq[i:i + 6], BOX_35) <= max_mm_35:
            # -35 box found at position i; look for -10 box at i+21 to i+25
            for spacer in range(15, 20):
                j = i + 6 + spacer
                if j + 6 <= n:
                    if _hamming(seq[j:j + 6], BOX_10) <= max_mm_10:
                        matches += 1
    return matches


# ---------------------------------------------------------------------------
# Repeat density utilities
# ---------------------------------------------------------------------------

def _repeat_density(sequence: str, min_len: int = 8) -> float:
    """
    Estimate the density of direct repeats and palindromes in the sequence.

    Uses a sliding-window hash approach: count distinct k-mers (k=min_len)
    that appear more than once.  Returns (redundant bp) / (total bp).

    NOTE: This is a rough approximation.  A proper implementation would use
    suffix arrays or RepeatMasker.
    """
    seq = sequence.upper()
    n = len(seq)
    if n < min_len * 2:
        return 0.0
    kmer_counts: dict[str, int] = {}
    for i in range(n - min_len + 1):
        km = seq[i:i + min_len]
        kmer_counts[km] = kmer_counts.get(km, 0) + 1
    # Count base positions covered by repeated kmers
    covered = 0
    for i in range(n - min_len + 1):
        if kmer_counts.get(seq[i:i + min_len], 1) > 1:
            covered += 1
    return covered / n


# ---------------------------------------------------------------------------
# A. MOBILITY / TRANSFER FEATURES
# ---------------------------------------------------------------------------

def compute_is_element_match(signal_score: Optional[float], **kwargs) -> FeatureResult:
    """
    IS element match — reuses the is_proximity signal score directly.

    The is_proximity signal (ISfinder BLAST) is the authoritative IS element
    assessment.  This feature acts as a pass-through into the transfer layer.
    """
    layer = "transfer_opportunity"
    name  = "is_element_match"
    if signal_score is None:
        return FeatureResult(
            feature_name=name, layer=layer,
            score=None, weight=_w(layer, name), available=False,
            evidence={}, source="signal_reuse",
            interpretation="ISfinder BLAST DB unavailable; signal not computed.",
        )
    return FeatureResult(
        feature_name=name, layer=layer,
        score=signal_score, weight=_w(layer, name), available=True,
        evidence={"source_signal": "is_proximity"},
        source="signal_reuse",
        interpretation=f"Score {signal_score:.3f} from ISfinder BLAST (reused).",
    )


def compute_integron_association(signal_score: Optional[float], **kwargs) -> FeatureResult:
    """Integron association — reuses the integron signal score."""
    layer = "transfer_opportunity"
    name  = "integron_association"
    if signal_score is None:
        return FeatureResult(
            feature_name=name, layer=layer,
            score=None, weight=_w(layer, name), available=False,
            evidence={}, source="signal_reuse",
            interpretation="Integron signal unavailable.",
        )
    return FeatureResult(
        feature_name=name, layer=layer,
        score=signal_score, weight=_w(layer, name), available=True,
        evidence={"source_signal": "integron"},
        source="signal_reuse",
        interpretation=f"Score {signal_score:.3f} from attC regex + INTEGRALL BLAST (reused).",
    )


def compute_conjugative_element(signal_score: Optional[float], **kwargs) -> FeatureResult:
    """Conjugative element homology — reuses the conjugative signal score."""
    layer = "transfer_opportunity"
    name  = "conjugative_element"
    if signal_score is None:
        return FeatureResult(
            feature_name=name, layer=layer,
            score=None, weight=_w(layer, name), available=False,
            evidence={}, source="signal_reuse",
            interpretation="Conjugative element protein DB unavailable.",
        )
    return FeatureResult(
        feature_name=name, layer=layer,
        score=signal_score, weight=_w(layer, name), available=True,
        evidence={"source_signal": "conjugative"},
        source="signal_reuse",
        interpretation=f"Score {signal_score:.3f} from conjugative element BLASTX (reused).",
    )


def compute_plasmid_context(**kwargs) -> FeatureResult:
    """
    Plasmid-origin probability (placeholder).

    Future implementation: classify the sequence as plasmid-borne or
    chromosomal using PlasClass, mlplasmids, or BLAST against RefSeq
    plasmid collection.  Current output is always unavailable.

    TODO: integrate PlasClass (Pellow et al. 2020, PLOS Comput Biol) or
          a BLAST-against-NCBI-plasmid-division approach.
    """
    layer = "transfer_opportunity"
    name  = "plasmid_context"
    return FeatureResult(
        feature_name=name, layer=layer,
        score=None, weight=_w(layer, name), available=False,
        evidence={"status": "placeholder"},
        source="placeholder",
        interpretation=(
            "Plasmid classifier not yet integrated.  "
            "Future: PlasClass or RefSeq plasmid BLAST."
        ),
    )


def compute_transposase_proximity(annotation_path: Optional[str] = None, **kwargs) -> FeatureResult:
    """
    Distance to nearest transposase/integrase annotation (placeholder).

    Future implementation: parse GFF3/GenBank annotation to identify
    transposase and integrase annotations, then compute bp distance from
    the query region to the nearest annotated mobile enzyme.  Requires
    genomic context input (GFF or GenBank).

    TODO: accept --annotation GFF argument in CLI and parse here.
    """
    layer = "transfer_opportunity"
    name  = "transposase_proximity"
    return FeatureResult(
        feature_name=name, layer=layer,
        score=None, weight=_w(layer, name), available=False,
        evidence={"status": "placeholder", "annotation_provided": annotation_path is not None},
        source="placeholder",
        interpretation=(
            "Genomic annotation not provided.  "
            "Future: compute distance to nearest transposase from GFF3/GenBank input."
        ),
    )


def compute_repeat_density(
    query: QuerySequence,
    flanking_sequence: Optional[str] = None,
    **kwargs,
) -> FeatureResult:
    """
    Repeat density in flanking sequence (or query if flanking unavailable).

    High repeat density in flanking regions is associated with mobile element
    insertion sites and recombination hotspots.  When flanking sequence is not
    provided, the query sequence itself is scanned as a proxy — this is weaker
    evidence and is noted in the interpretation.

    Scoring: repeat density (fraction of bases covered by repeated 8-mers)
    mapped to 0–1 using a 0.3 saturation threshold (30% repeat density → 1.0).
    """
    layer = "transfer_opportunity"
    name  = "repeat_density"
    seq_to_scan = flanking_sequence if flanking_sequence else query.sequence
    context_label = "flanking sequence" if flanking_sequence else "query sequence (no flanking provided)"

    density = _repeat_density(seq_to_scan)
    # 30% repeat density saturates the score
    score = min(density / 0.30, 1.0)

    return FeatureResult(
        feature_name=name, layer=layer,
        score=score, weight=_w(layer, name), available=True,
        evidence={
            "density": round(density, 4),
            "scanned": context_label,
            "min_kmer_len": 8,
        },
        source="computed",
        interpretation=(
            f"Repeat density {density:.1%} in {context_label}. "
            f"{'High' if density > 0.15 else 'Low'} density of repeated 8-mers."
            + ("  (Proxy: flanking sequence not provided.)" if not flanking_sequence else "")
        ),
    )


# ---------------------------------------------------------------------------
# B. ESTABLISHMENT / HOST COMPATIBILITY FEATURES
# ---------------------------------------------------------------------------

def compute_gc_deviation(signal_score: Optional[float], evidence: Optional[dict] = None, **kwargs) -> FeatureResult:
    """
    GC content deviation — reuses the gc_content signal score.

    A large deviation from host GC content is both a classical marker of
    foreign origin and a potential barrier to host replication fidelity.
    """
    layer = "establishment"
    name  = "gc_deviation"
    if signal_score is None:
        return FeatureResult(
            feature_name=name, layer=layer,
            score=None, weight=_w(layer, name), available=False,
            evidence={}, source="signal_reuse",
            interpretation="GC content signal unavailable.",
        )
    ev = evidence or {}
    return FeatureResult(
        feature_name=name, layer=layer,
        score=signal_score, weight=_w(layer, name), available=True,
        evidence={
            "source_signal": "gc_content",
            "query_gc_pct": ev.get("query_gc"),
            "host_gc_pct":  ev.get("host_gc"),
            "deviation_pct": ev.get("deviation_pct"),
        },
        source="signal_reuse",
        interpretation=(
            f"GC deviation {ev.get('deviation_pct', '?')}% "
            f"(query {ev.get('query_gc', '?')}%, "
            f"host {ev.get('host_gc', '?')}%)."
        ),
    )


def compute_codon_usage_distance(query: QuerySequence, host: HostProfile, **kwargs) -> FeatureResult:
    """
    Codon usage distance between query element and host.

    Computes codon frequencies across all three forward reading frames of the
    query, then calculates the Total Variation Distance (TVD = L1/2) against
    the nearest host-specific reference codon table.

    High TVD: element uses significantly different codons than the host →
    may indicate foreign origin and/or reduced expression fitness.

    Reference tables available: E. coli K-12, B. subtilis 168.
    Falls back to an equal-usage null reference if host is unrecognised.

    NOTE: Frequencies are computed from raw sequence (all frames), not from
    predicted CDS only.  This is a first-pass approximation.  A full
    implementation would use annotated CDS coordinates.

    Scoring: TVD mapped directly to score (TVD ∈ [0, 1] by definition).
    """
    layer = "establishment"
    name  = "codon_usage_distance"

    query_freqs = _query_codon_freqs(query.sequence)
    ref_name, ref_freqs = _select_reference(host.identifier)

    l1 = sum(abs(query_freqs[c] - ref_freqs[c]) for c in _ALL_SENSE_CODONS)
    tvd = l1 / 2.0          # Total Variation Distance ∈ [0, 1]

    # High TVD → high concern (element looks foreign relative to host codon profile)
    score = min(tvd * 2.0, 1.0)  # scale: TVD of 0.5 → max concern

    # Identify the most-diverged codons for the evidence dict
    diffs = sorted(
        ((c, round(query_freqs[c] - ref_freqs[c], 4)) for c in _ALL_SENSE_CODONS),
        key=lambda x: abs(x[1]), reverse=True
    )[:5]

    return FeatureResult(
        feature_name=name, layer=layer,
        score=score, weight=_w(layer, name), available=True,
        evidence={
            "tvd": round(tvd, 4),
            "reference_used": ref_name,
            "top_diverged_codons": [{"codon": c, "delta": d} for c, d in diffs],
            "note": (
                "Approximate: codon frequencies computed across all 3 forward frames. "
                "Replace reference tables with CoCoPUTs data for higher accuracy."
            ),
        },
        source="computed",
        interpretation=(
            f"Codon usage TVD = {tvd:.3f} vs {ref_name}. "
            f"{'Large divergence suggests different evolutionary origin.' if tvd > 0.25 else 'Modest divergence from host codon profile.'}"
        ),
    )


def compute_taxonomic_distance(
    donor_taxon: Optional[str] = None,
    recipient_taxon: Optional[str] = None,
    **kwargs,
) -> FeatureResult:
    """
    Categorical taxonomic distance between donor and recipient organisms.

    Compares lineage strings from a small curated lookup of common laboratory
    and environmental organisms.  Returns a distance in [0, 1] based on the
    deepest shared taxonomic rank.

    This is a categorical proxy, NOT an evolutionary distance.  It is most
    useful for flagging very distant transfers (e.g., plant → bacterium) and
    should not be over-interpreted for closely related organisms.

    If either taxon is absent or unrecognised, the feature degrades gracefully.
    """
    layer = "establishment"
    name  = "taxonomic_distance"

    if not donor_taxon or not recipient_taxon:
        return FeatureResult(
            feature_name=name, layer=layer,
            score=None, weight=_w(layer, name), available=False,
            evidence={"donor": donor_taxon, "recipient": recipient_taxon},
            source="computed",
            interpretation="Donor and/or recipient taxon not provided; feature unavailable.",
        )

    score, detail = _categorical_tax_distance(donor_taxon, recipient_taxon)
    return FeatureResult(
        feature_name=name, layer=layer,
        score=score, weight=_w(layer, name), available=True,
        evidence=detail,
        source="computed",
        interpretation=(
            f"Categorical taxonomic distance {score:.2f}: "
            f"{detail.get('shared_rank', 'unknown')} is the deepest shared rank."
        ),
    )


# Curated lineage strings for common lab / environmental organisms.
# Format: semicolon-separated ranks from domain to genus.
_TAX_DB: dict[str, str] = {
    "ecoli":        "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia",
    "salmonella":   "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Salmonella",
    "klebsiella":   "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Klebsiella",
    "pseudomonas":  "Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas",
    "acinetobacter":"Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Moraxellaceae;Acinetobacter",
    "bacillus":     "Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus",
    "staphylococcus":"Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus",
    "streptococcus":"Bacteria;Firmicutes;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus",
    "clostridium":  "Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridium",
    "mycobacterium":"Bacteria;Actinobacteria;Actinomycetia;Corynebacteriales;Mycobacteriaceae;Mycobacterium",
    "streptomyces": "Bacteria;Actinobacteria;Actinomycetia;Streptomycetales;Streptomycetaceae;Streptomyces",
    "saccharomyces":"Eukaryota;Fungi;Ascomycota;Saccharomycetes;Saccharomycetales;Saccharomycetaceae;Saccharomyces",
    "human":        "Eukaryota;Metazoa;Chordata;Mammalia;Primates;Hominidae;Homo",
    "arabidopsis":  "Eukaryota;Plantae;Tracheophyta;Magnoliopsida;Brassicales;Brassicaceae;Arabidopsis",
}

# Rank depth → categorical distance score
_RANK_DISTANCE: dict[str, float] = {
    "genus":    0.10,
    "family":   0.30,
    "order":    0.50,
    "class":    0.65,
    "phylum":   0.80,
    "domain":   0.95,
    "none":     1.00,   # no shared ranks at all (e.g., prokaryote vs eukaryote)
}


def _normalise_taxon(name: str) -> str:
    return name.lower().replace(" ", "").replace("_", "").replace(".", "")


def _categorical_tax_distance(donor: str, recipient: str) -> tuple[float, dict]:
    d_key = _normalise_taxon(donor)
    r_key = _normalise_taxon(recipient)

    # Allow partial matches (e.g., "e.coli" → "ecoli")
    d_lineage = next((v for k, v in _TAX_DB.items() if k in d_key or d_key in k), None)
    r_lineage = next((v for k, v in _TAX_DB.items() if k in r_key or r_key in k), None)

    if d_lineage is None or r_lineage is None:
        return (0.5, {
            "donor_recognised":     d_lineage is not None,
            "recipient_recognised": r_lineage is not None,
            "note": "One or both taxa not in curated lookup; distance set to 0.5 (uncertain).",
            "shared_rank": "unknown",
        })

    d_ranks = d_lineage.split(";")
    r_ranks = r_lineage.split(";")
    rank_labels = ["domain", "phylum", "class", "order", "family", "genus"]
    shared_depth = 0
    for d_r, r_r in zip(d_ranks, r_ranks):
        if d_r == r_r:
            shared_depth += 1
        else:
            break

    # Map shared_depth to rank label (1=domain, 2=phylum, ...)
    if shared_depth == 0:
        shared_rank = "none"
    elif shared_depth >= len(rank_labels):
        shared_rank = "genus"           # same genus or deeper
    else:
        shared_rank = rank_labels[shared_depth - 1]

    distance = _RANK_DISTANCE.get(shared_rank, 0.5)
    return (distance, {
        "donor": donor,
        "recipient": recipient,
        "donor_lineage": d_lineage,
        "recipient_lineage": r_lineage,
        "shared_rank": shared_rank,
        "note": "Categorical proxy only; not an evolutionary distance metric.",
    })


def compute_promoter_plausibility(query: QuerySequence, **kwargs) -> FeatureResult:
    """
    Prokaryotic sigma-70 promoter motif plausibility (heuristic proxy).

    Scans for sigma-70-like -10 (TATAAT) and -35 (TTGACA) box pairs with
    15–19 bp spacer spacing.  A higher count suggests the element contains
    sequences that may function as promoters in a gram-negative host.

    IMPORTANT CAVEATS:
    - This is a simple pattern match, not a mechanistic expression predictor.
    - False positives are common in random sequence.
    - A positive result indicates potential for expression; it does NOT
      confirm that transcription will occur.
    - Eukaryotic promoter architectures are not assessed here.

    Scoring: min(count / 5, 1.0)  (5+ promoter-like motifs → maximum score).
    """
    layer = "establishment"
    name  = "promoter_plausibility"
    count = _count_promoter_like(query.sequence)
    score = min(count / 5.0, 1.0)

    return FeatureResult(
        feature_name=name, layer=layer,
        score=score, weight=_w(layer, name), available=True,
        evidence={
            "sigma70_like_count": count,
            "method": "-10/-35 consensus scan with ≤2 mismatches each",
            "caveat": "Heuristic proxy only. False positives expected.",
        },
        source="computed",
        interpretation=(
            f"{count} putative sigma-70-like promoter element(s) detected. "
            + ("May support expression in gram-negative host." if count > 0
               else "No strong promoter-like sequences found.")
            + "  (Heuristic proxy — not a definitive expression predictor.)"
        ),
    )


def compute_sequence_complexity(query: QuerySequence, **kwargs) -> FeatureResult:
    """
    Sequence length and estimated functional gene count as establishment burden proxies.

    Longer, gene-dense elements impose a higher metabolic burden on the recipient
    and require more compatibility factors to establish stably.  This feature
    scores the element as more 'complex' (and thus carrying a higher potential
    functional payload) as length and gene count increase.

    Scoring: combined length score (saturation at 20 kbp) and ORF count score
    (saturation at 20 complete ORFs ≥ 300 bp), averaged.
    """
    layer = "establishment"
    name  = "sequence_complexity"

    complete_orfs, _ = _scan_orfs(query.sequence, min_nt=300)
    orf_count = len(complete_orfs)

    length_score = min(query.length / 20_000, 1.0)
    orf_score    = min(orf_count / 20.0, 1.0)
    score        = (length_score + orf_score) / 2.0

    return FeatureResult(
        feature_name=name, layer=layer,
        score=score, weight=_w(layer, name), available=True,
        evidence={
            "sequence_length_bp": query.length,
            "complete_orfs_ge300bp": orf_count,
            "length_score": round(length_score, 3),
            "orf_score": round(orf_score, 3),
        },
        source="computed",
        interpretation=(
            f"Sequence {query.length:,} bp with {orf_count} complete ORF(s) ≥300 bp. "
            f"{'High' if score > 0.5 else 'Moderate' if score > 0.2 else 'Low'} complexity."
        ),
    )


# ---------------------------------------------------------------------------
# C. FUNCTIONAL CONSEQUENCE FEATURES
# ---------------------------------------------------------------------------

def compute_prophage_context(signal_score: Optional[float], **kwargs) -> FeatureResult:
    """Prophage context — reuses the prophage signal score."""
    layer = "consequence"
    name  = "prophage_context"
    if signal_score is None:
        return FeatureResult(
            feature_name=name, layer=layer,
            score=None, weight=_w(layer, name), available=False,
            evidence={}, source="signal_reuse",
            interpretation="Prophage signal unavailable (PHASTER API and local DB absent).",
        )
    return FeatureResult(
        feature_name=name, layer=layer,
        score=signal_score, weight=_w(layer, name), available=True,
        evidence={"source_signal": "prophage"},
        source="signal_reuse",
        interpretation=f"Score {signal_score:.3f} from PHASTER prophage analysis (reused).",
    )


def compute_amr_content(data_dir=None, **kwargs) -> FeatureResult:
    """
    Antimicrobial resistance (AMR) gene content (placeholder).

    Future implementation: BLAST against CARD (Comprehensive Antibiotic
    Resistance Database) protein homolog models, or use RGI (Resistance
    Gene Identifier) directly.

    TODO: integrate CARD/RGI.  Key markers: beta-lactamases (blaTEM, blaCTX),
          aminoglycoside modifying enzymes, efflux pump regulators (mexAB).

    This feature has the highest weight in the consequence layer because
    AMR-carrying elements are the primary regulatory concern for engineered
    organisms in clinical and environmental contexts.
    """
    layer = "consequence"
    name  = "amr_content"
    return FeatureResult(
        feature_name=name, layer=layer,
        score=None, weight=_w(layer, name), available=False,
        evidence={"status": "placeholder"},
        source="placeholder",
        interpretation=(
            "CARD/RGI integration not yet available.  "
            "Future: BLAST against CARD protein homolog models."
        ),
    )


def compute_virulence_flags(data_dir=None, **kwargs) -> FeatureResult:
    """
    Virulence, toxin, and persistence gene flagging (placeholder).

    Future implementation: BLAST against VFDB (Virulence Factor Database)
    or PATRIC virulence factor collection; search for known toxin-antitoxin
    (TA) system signatures; check for adhesin and invasion-associated genes.

    TODO: integrate VFDB BLAST or annotate via Prokka + VF annotation.
    """
    layer = "consequence"
    name  = "virulence_flags"
    return FeatureResult(
        feature_name=name, layer=layer,
        score=None, weight=_w(layer, name), available=False,
        evidence={"status": "placeholder"},
        source="placeholder",
        interpretation=(
            "VFDB/virulence annotation not yet integrated.  "
            "Future: BLAST against VFDB or Prokka + VF annotation."
        ),
    )


def compute_gene_completeness(query: QuerySequence, **kwargs) -> FeatureResult:
    """
    ORF integrity assessment.

    Distinguishes complete ORFs (ATG + in-frame stop, ≥90 bp) from partial
    ORFs truncated at the sequence boundary.  High completeness indicates that
    functional units appear intact and capable of producing full-length
    proteins in a recipient — this increases functional consequence risk.

    Scoring: fraction of (complete / (complete + partial)).
    When no ORFs are detected, returns score 0.0 (no detectable payload).
    """
    layer = "consequence"
    name  = "gene_completeness"

    complete, partial = _scan_orfs(query.sequence, min_nt=90)
    total = len(complete) + len(partial)

    if total == 0:
        score = 0.0
        interp = "No ORFs ≥90 bp detected; no apparent functional payload."
    else:
        score = len(complete) / total
        interp = (
            f"{len(complete)} complete ORF(s) and {len(partial)} partial ORF(s) ≥90 bp. "
            f"Completeness ratio {score:.2f}."
        )

    return FeatureResult(
        feature_name=name, layer=layer,
        score=score, weight=_w(layer, name), available=True,
        evidence={
            "complete_orfs_ge90bp": len(complete),
            "partial_orfs_ge90bp":  len(partial),
            "total_orfs":           total,
            "completeness_ratio":   round(score, 3),
        },
        source="computed",
        interpretation=interp,
    )


def compute_payload_count(query: QuerySequence, **kwargs) -> FeatureResult:
    """
    Number of candidate functional genes (multi-gene payload burden).

    Counts complete ORFs ≥300 bp (≈100 aa, a conservative minimum for a
    functional protein domain).  A higher count indicates a more complex
    functional payload with greater potential for diverse effects in a
    recipient organism.

    Scoring: min(count / 10, 1.0)  (10+ complete ORFs ≥300 bp → maximum score).
    """
    layer = "consequence"
    name  = "payload_count"

    complete, _ = _scan_orfs(query.sequence, min_nt=300)
    count = len(complete)
    score = min(count / 10.0, 1.0)

    return FeatureResult(
        feature_name=name, layer=layer,
        score=score, weight=_w(layer, name), available=True,
        evidence={
            "complete_orfs_ge300bp": count,
            "saturation_at": 10,
        },
        source="computed",
        interpretation=(
            f"{count} complete ORF(s) ≥300 bp detected. "
            f"{'High' if count >= 5 else 'Moderate' if count >= 2 else 'Low'} multi-gene payload burden."
        ),
    )
