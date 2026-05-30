from dataclasses import dataclass, field
from typing import List


# ---------------------------------------------------------------------------
# VariantSet dataclass
# ---------------------------------------------------------------------------

@dataclass
class Variant:
    """Single variant record extracted from an alignment."""
    position: int          # 0-based position in seq1 (proxy) coordinate
    type: str              # "SNP", "SMALL_INS", "LARGE_INS", "SMALL_DEL", "LARGE_DEL"
    ref_len: int           # length in seq1 (0 for pure insertions)
    alt_len: int           # length in seq2 (0 for pure deletions)
    ref_allele: str        # actual bases from seq1 (empty string for insertions)
    alt_allele: str        # actual bases from seq2 (empty string for deletions)
    context_seq: str       # ±10 bp context from seq1
    chrom: str = "chr1"
    impact_category: str = ""
    impact_score: int = 0
    weight: int = 1
    variant_class: str = ""


@dataclass
class VariantSet:
    """Container for all variants parsed from a single alignment."""
    variants: List[Variant] = field(default_factory=list)
    snps: List[Variant] = field(default_factory=list)
    small_indels: List[Variant] = field(default_factory=list)   # 1–50 bp
    large_indels: List[Variant] = field(default_factory=list)   # >50 bp

    @property
    def total(self) -> int:
        return len(self.variants)


# ---------------------------------------------------------------------------
# Weights (must be defined before classify_variants)
# ---------------------------------------------------------------------------

WEIGHTS = {
    "SNV": 1,
    "SMALL_INS": 3,
    "SMALL_DEL": 3,
    "LARGE_INS": 10,
    "LARGE_DEL": 10,
    "SV_INS": 50,
    "SV_DEL": 50,
}

# ---------------------------------------------------------------------------
# classify_variants() — public API
# ---------------------------------------------------------------------------

_SMALL_INDEL_MAX = 50   # bp; <= this → small indel, > this → large indel


def classify_variants(alignment, chrom: str = "chr1") -> VariantSet:
    """Parse an Alignment namedtuple and return a VariantSet.

    Parameters
    ----------
    alignment : Alignment
        namedtuple from sequence_aligner.align_sequences with fields
        aligned_seq1, aligned_seq2, cigar.
    chrom : str
        Chromosome / sequence label to attach to each variant record.

    Returns
    -------
    VariantSet
        dataclass holding all variants plus convenience sub-lists.
    """
    al1 = alignment.aligned_seq1
    al2 = alignment.aligned_seq2
    seq1_ungapped = al1.replace("-", "")  # used for context extraction

    vs = VariantSet()

    # Walk the pairwise alignment column by column.
    seq1_pos = 0   # 0-based position in original (ungapped) seq1
    i = 0
    n = len(al1)

    while i < n:
        c1 = al1[i]
        c2 = al2[i]

        if c1 == "-":
            # Insertion in seq2 relative to seq1 — consume run.
            ins_bases: list[str] = []
            j = i
            while j < n and al1[j] == "-":
                ins_bases.append(al2[j])
                j += 1
            ins_seq = "".join(ins_bases)
            ins_len = len(ins_seq)
            vtype = "SMALL_INS" if ins_len <= _SMALL_INDEL_MAX else "LARGE_INS"
            ctx = _context(seq1_ungapped, seq1_pos, 10)
            v = Variant(
                position=seq1_pos,
                type=vtype,
                ref_len=0,
                alt_len=ins_len,
                ref_allele="",
                alt_allele=ins_seq,
                context_seq=ctx,
                chrom=chrom,
                variant_class=vtype,
                weight=WEIGHTS.get(vtype, 1),
            )
            vs.variants.append(v)
            if ins_len <= _SMALL_INDEL_MAX:
                vs.small_indels.append(v)
            else:
                vs.large_indels.append(v)
            i = j

        elif c2 == "-":
            # Deletion in seq2 relative to seq1 — consume run.
            del_bases: list[str] = []
            j = i
            while j < n and al2[j] == "-":
                del_bases.append(al1[j])
                j += 1
            del_seq = "".join(del_bases)
            del_len = len(del_seq)
            vtype = "SMALL_DEL" if del_len <= _SMALL_INDEL_MAX else "LARGE_DEL"
            ctx = _context(seq1_ungapped, seq1_pos, 10)
            v = Variant(
                position=seq1_pos,
                type=vtype,
                ref_len=del_len,
                alt_len=0,
                ref_allele=del_seq,
                alt_allele="",
                context_seq=ctx,
                chrom=chrom,
                variant_class=vtype,
                weight=WEIGHTS.get(vtype, 1),
            )
            vs.variants.append(v)
            if del_len <= _SMALL_INDEL_MAX:
                vs.small_indels.append(v)
            else:
                vs.large_indels.append(v)
            seq1_pos += del_len
            i = j

        else:
            # Both columns are bases.
            if c1 != c2:
                ctx = _context(seq1_ungapped, seq1_pos, 10)
                v = Variant(
                    position=seq1_pos,
                    type="SNP",
                    ref_len=1,
                    alt_len=1,
                    ref_allele=c1,
                    alt_allele=c2,
                    context_seq=ctx,
                    chrom=chrom,
                    variant_class="SNV",
                    weight=WEIGHTS.get("SNV", 1),
                )
                vs.variants.append(v)
                vs.snps.append(v)
            seq1_pos += 1
            i += 1

    return vs


def _context(seq: str, pos: int, radius: int = 10) -> str:
    """Return up to ±radius bases around *pos* in *seq*."""
    start = max(0, pos - radius)
    end = min(len(seq), pos + radius + 1)
    return seq[start:end]


def classify_variant(variant):
    vtype = variant["type"]
    length = variant.get("length", 1)

    if vtype == "SNV":
        return "SNV"
    elif vtype == "INS":
        if length <= 50:
            return "SMALL_INS"
        elif length <= 500:
            return "LARGE_INS"
        else:
            return "SV_INS"
    elif vtype == "DEL":
        if length <= 50:
            return "SMALL_DEL"
        elif length <= 500:
            return "LARGE_DEL"
        else:
            return "SV_DEL"
    return "SNV"


def classify_all_variants(variants):
    for v in variants:
        v["variant_class"] = classify_variant(v)
        v["weight"] = WEIGHTS.get(v["variant_class"], 1)
    return variants


def compute_burden(variants, total_genome_bp):
    class_counts = {}
    for v in variants:
        vc = v.get("variant_class", "SNV")
        class_counts[vc] = class_counts.get(vc, 0) + 1

    total_edits = sum(class_counts.values())
    weighted_burden = sum(
        count * WEIGHTS.get(vc, 1) for vc, count in class_counts.items()
    )

    total_mb = max(total_genome_bp / 1_000_000, 1e-6)
    normalized_burden = weighted_burden / total_mb

    return {
        "class_counts": class_counts,
        "total_edits": total_edits,
        "weighted_burden": weighted_burden,
        "normalized_burden_per_mb": round(normalized_burden, 2),
        "total_genome_bp": total_genome_bp,
    }
