import math
from typing import Any, Dict, List

from .sequence_analyzer import (
    extract_window, compute_repeat_score, compute_complexity_score,
    SequenceFeatures,
)
from .feature_annotator import gene_density_score, regulatory_proximity_score
from .config import (
    SCORE_WEIGHTS, SIZE_PENALTY_THRESHOLD_BP,
    TIER_HIGH, TIER_MODERATE,
    GC_OPEN_MIN, GC_OPEN_MAX,
)


def _chromatin_score(features: SequenceFeatures) -> float:
    """Proxy for open-chromatin permissiveness.

    GC content in the moderate range (GC_OPEN_MIN–GC_OPEN_MAX) and low repeat
    density are associated with gene-desert / open chromatin that tolerates
    large inserts.  Returns a value in [0, 1].
    """
    gc = features.gc_content
    # Triangular reward: peak at midpoint of the open range
    mid = (GC_OPEN_MIN + GC_OPEN_MAX) / 2.0
    half = (GC_OPEN_MAX - GC_OPEN_MIN) / 2.0
    gc_score = max(0.0, 1.0 - abs(gc - mid) / half) if half > 0 else 0.0

    # Penalise high repeat density (>0.3 is poor)
    repeat_penalty = min(features.repeat_density / 0.3, 1.0)
    repeat_score = 1.0 - repeat_penalty

    return (gc_score + repeat_score) / 2.0


def _sequence_complexity_score(features: SequenceFeatures) -> float:
    """Penalise structural motifs that stall polymerases or cause artefacts.

    G-quadruplexes, homopolymer runs, and dense inverted repeats are penalised.
    Returns a value in [0, 1] where 1 = no problematic motifs.
    """
    seq_len = max(features.sequence_length, 1)

    # Normalise counts relative to sequence length (per kb)
    per_kb = seq_len / 1000.0

    g4_rate = features.g_quadruplex_count / per_kb
    hp_rate = features.homopolymer_count / per_kb
    ir_rate = features.inverted_repeat_count / per_kb

    # Penalty caps: >2 G4/kb, >5 homopolymers/kb, >50 inverted repeats/kb = worst
    # (inverted repeat cap is high because even random sequence contains many
    #  coincidental 12-bp palindromes; only dense clusters are penalised)
    g4_penalty = min(g4_rate / 2.0, 1.0)
    hp_penalty = min(hp_rate / 5.0, 1.0)
    ir_penalty = min(ir_rate / 50.0, 1.0)

    penalty = (g4_penalty * 0.4 + hp_penalty * 0.4 + ir_penalty * 0.2)
    return 1.0 - penalty


def _gene_density_score(genomic_features: List[Dict[str, Any]]) -> float:
    """Penalise loci where the insert would disrupt annotated genes or exons.

    Returns a value in [0, 1].
    """
    if not genomic_features:
        return 1.0  # no annotation = assume safe

    gene_hits = sum(1 for f in genomic_features if f["feature_type"] == "gene")
    exon_hits = sum(1 for f in genomic_features if f["feature_type"] == "exon")
    reg_hits = sum(1 for f in genomic_features
                   if f["feature_type"] == "regulatory")

    # Each overlapping gene is a major penalty; exons are very bad;
    # regulatory elements carry moderate risk
    penalty = min(gene_hits * 0.25 + exon_hits * 0.40 + reg_hits * 0.10, 1.0)
    return 1.0 - penalty


def _size_penalty(insert_size_bp: int) -> float:
    """Log-linear penalty for insert size beyond SIZE_PENALTY_THRESHOLD_BP.

    Returns a value in [0, 1] where 1 = no penalty (small insert).
    """
    threshold = SIZE_PENALTY_THRESHOLD_BP
    if insert_size_bp <= threshold:
        return 1.0
    # Penalty grows log-linearly from threshold to 100 kb (full penalty)
    max_size = 100_000
    excess = insert_size_bp - threshold
    max_excess = max_size - threshold
    # log scale: log(1 + excess) / log(1 + max_excess)
    penalty = math.log1p(excess) / math.log1p(max_excess)
    return max(0.0, 1.0 - penalty)


def score_tolerance(
    seq_features: SequenceFeatures,
    genomic_features: List[Dict[str, Any]],
    insert_size_bp: int,
) -> Dict[str, Any]:
    """Compute composite insertion tolerance for a locus.

    Parameters
    ----------
    seq_features:
        SequenceFeatures returned by analyze_sequence().
    genomic_features:
        List of GenomicFeature dicts from annotate_locus().
    insert_size_bp:
        Size of the intended insert in base pairs.

    Returns
    -------
    dict with keys:
        chromatin_score, sequence_complexity_score, gene_density_score,
        size_penalty, composite_tolerance (all 0–1),
        tolerance_tier ("high" | "moderate" | "low")
    """
    cs = _chromatin_score(seq_features)
    sc = _sequence_complexity_score(seq_features)
    gd = _gene_density_score(genomic_features)
    sp = _size_penalty(insert_size_bp)

    w = SCORE_WEIGHTS
    composite = (
        w["chromatin_score"] * cs
        + w["sequence_complexity_score"] * sc
        + w["gene_density_score"] * gd
        + w["size_penalty"] * sp
    )
    composite = max(0.0, min(1.0, composite))

    if composite > TIER_HIGH:
        tier = "high"
    elif composite >= TIER_MODERATE:
        tier = "moderate"
    else:
        tier = "low"

    return {
        "chromatin_score": round(cs, 4),
        "sequence_complexity_score": round(sc, 4),
        "gene_density_score": round(gd, 4),
        "size_penalty": round(sp, 4),
        "composite_tolerance": round(composite, 4),
        "tolerance_tier": tier,
    }


def score_window(sequences, genes, chrom, win_start, win_end):
    seq = extract_window(sequences, chrom, win_start, win_end)

    gd_score = gene_density_score(genes, chrom, win_start, win_end)
    reg_score = regulatory_proximity_score(sequences, chrom, win_start, win_end)
    rep_score = compute_repeat_score(seq)
    cplx_score = compute_complexity_score(seq)
    total = gd_score + reg_score + rep_score + cplx_score

    return {
        "chrom": chrom,
        "start": win_start,
        "end": win_end,
        "gene_density_score": gd_score,
        "regulatory_score": reg_score,
        "repeat_score": rep_score,
        "complexity_score": cplx_score,
        "total_score": total,
    }


def scan_locus(sequences, genes, locus, window_size=1000, step=200, flank=50000):
    chrom = locus["chrom"]
    center = (locus["start"] + locus["end"]) // 2
    scan_start = max(0, center - flank)
    chrom_len = len(sequences.get(chrom, ""))
    scan_end = min(chrom_len if chrom_len else center + flank, center + flank)

    results = []
    pos = scan_start
    while pos + window_size <= scan_end:
        row = score_window(sequences, genes, chrom, pos, pos + window_size)
        row["locus_name"] = locus["name"]
        results.append(row)
        pos += step

    return results
