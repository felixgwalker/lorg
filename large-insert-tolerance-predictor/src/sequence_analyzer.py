import random
import math


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
