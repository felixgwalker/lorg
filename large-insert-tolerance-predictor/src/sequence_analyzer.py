import random
import math
import re
from dataclasses import dataclass
from typing import Dict, Any


@dataclass
class SequenceFeatures:
    """Sequence-level features computed for a genomic interval."""
    gc_content: float = 0.0          # fraction 0–1
    repeat_density: float = 0.0      # fraction of bases in simple repeats
    g_quadruplex_count: int = 0      # G4 motif occurrences
    inverted_repeat_count: int = 0   # palindromes > 8 bp
    homopolymer_count: int = 0       # runs of ≥ 6 identical bases
    sequence_length: int = 0


def analyze_sequence(seq: str) -> SequenceFeatures:
    """Compute sequence features from a DNA string.

    Returns a SequenceFeatures dataclass with:
    - gc_content: fraction of G/C bases
    - repeat_density: fraction of bases covered by simple tandem repeats
      (periods 1–6 detected by a sliding window approach)
    - g_quadruplex_count: occurrences of G₃+N₁₋₇G₃+N₁₋₇G₃+N₁₋₇G₃+
    - inverted_repeat_count: palindromic sequences > 8 bp
    - homopolymer_count: runs of ≥ 6 identical bases
    """
    if not seq:
        return SequenceFeatures()

    seq_upper = seq.upper()
    n = len(seq_upper)

    # GC content
    gc = sum(1 for b in seq_upper if b in ("G", "C"))
    gc_content = gc / n

    # Simple repeat density: tandem repeats of period 1–6.
    # For each period p, scan with a sliding window: if the next p bases
    # match the current p bases, mark all those positions as repetitive.
    covered = bytearray(n)  # 0/1 per position
    for period in range(1, 7):
        i = 0
        while i + period * 2 <= n:
            unit = seq_upper[i:i + period]
            # Extend run as far as it goes
            j = i + period
            while j + period <= n and seq_upper[j:j + period] == unit:
                j += period
            run_len = j - i
            if run_len >= period * 2:  # at least 2 copies = tandem repeat
                for k in range(i, j):
                    covered[k] = 1
                i = j
            else:
                i += 1
    repeat_density = sum(covered) / n

    # G-quadruplex: G{3,}N{1,7}G{3,}N{1,7}G{3,}N{1,7}G{3,}
    g4_pattern = re.compile(r"G{3,}[ACGTN]{1,7}G{3,}[ACGTN]{1,7}G{3,}[ACGTN]{1,7}G{3,}",
                            re.IGNORECASE)
    g_quadruplex_count = sum(1 for _ in g4_pattern.finditer(seq_upper))

    # Inverted repeats (palindromes) > 8 bp.
    # We search for sequences of length 8..min(30, n//2) that appear
    # alongside their reverse complement within a local window.
    complement = str.maketrans("ACGTN", "TGCAN")
    inverted_repeat_count = 0
    min_arm = 6   # minimum arm length (palindrome total = 2 * arm = ≥ 12 bp)
    max_arm = 20
    # Slide along the sequence looking for even-length palindromic centres
    for arm in range(min_arm, max_arm + 1):
        for i in range(n - 2 * arm + 1):
            left = seq_upper[i:i + arm]
            right = seq_upper[i + arm:i + 2 * arm]
            rc_right = right.translate(complement)[::-1]
            if left == rc_right:
                inverted_repeat_count += 1

    # Homopolymer runs ≥ 6
    homopolymer_count = 0
    i = 0
    while i < n:
        j = i + 1
        while j < n and seq_upper[j] == seq_upper[i]:
            j += 1
        if j - i >= 6:
            homopolymer_count += 1
        i = j

    return SequenceFeatures(
        gc_content=gc_content,
        repeat_density=repeat_density,
        g_quadruplex_count=g_quadruplex_count,
        inverted_repeat_count=inverted_repeat_count,
        homopolymer_count=homopolymer_count,
        sequence_length=n,
    )


def extract_window(sequences, chrom, start, end):
    seq = sequences.get(chrom, "")
    if not seq:
        return ""
    start = max(0, start)
    end = min(len(seq), end)
    return seq[start:end]


def compute_repeat_score(seq):
    if not seq:
        return 10
    seq_upper = seq.upper()
    length = len(seq_upper)

    homopolymer_count = 0
    i = 0
    while i < length:
        j = i + 1
        while j < length and seq_upper[j] == seq_upper[i]:
            j += 1
        run_len = j - i
        if run_len >= 4:
            homopolymer_count += run_len
        i = j

    k = 6
    kmers = {}
    for idx in range(length - k + 1):
        km = seq_upper[idx:idx + k]
        kmers[km] = kmers.get(km, 0) + 1

    repeat_kmer_count = sum(v for v in kmers.values() if v > 1) * k
    repeat_fraction = (homopolymer_count + repeat_kmer_count) / max(length, 1)
    repeat_fraction = min(repeat_fraction, 1.0)

    if repeat_fraction < 0.10:
        return 20
    elif repeat_fraction <= 0.30:
        return 10
    else:
        return 0


def compute_complexity_score(seq):
    if not seq:
        return 10
    seq_upper = seq.upper()
    length = len(seq_upper)
    k = 4
    if length < k:
        return 10
    total_kmers = length - k + 1
    unique_kmers = len(set(seq_upper[i:i + k] for i in range(total_kmers)))
    complexity = unique_kmers / max(total_kmers, 1)

    if complexity > 0.7:
        return 20
    elif complexity >= 0.5:
        return 10
    else:
        return 0


def make_demo_sequences(genome_size=200000):
    random.seed(42)
    bases = "ACGT"
    chroms = {}

    def random_seq(n):
        return "".join(random.choices(bases, k=n))

    seq1 = list(random_seq(genome_size))
    for pos in range(30000, 30100):
        seq1[pos] = "A"
    chroms["chr1"] = "".join(seq1)

    seq2 = list(random_seq(genome_size))
    for pos in range(60000, 60060):
        seq2[pos] = seq2[pos % 6 + 60000 - 6] if pos >= 6 else seq2[0]
    chroms["chr2"] = "".join(seq2)

    return chroms
