"""Build k-mer index of 20-mers from genome FASTA."""

import random
from pathlib import Path

KMER_LEN = 20

# Ambiguity codes included: N stays N; standard complement table
_COMP_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_complement(seq: str) -> str:
    """Return reverse complement; N maps to N, all standard bases handled."""
    return seq.translate(_COMP_TABLE)[::-1]


def build_genome_index(genome: dict[str, str]) -> dict[str, list[tuple[str, int, str]]]:
    """Return dict mapping 20-mer -> list of (chrom, pos, strand).

    Indexes both forward and reverse-complement strands.  Kmers containing N
    are skipped so the index only holds unambiguous 20-mers.
    """
    index: dict[str, list[tuple[str, int, str]]] = {}
    for chrom, seq in genome.items():
        seq = seq.upper()
        rc = reverse_complement(seq)
        n = len(seq)
        for i in range(n - KMER_LEN + 1):
            fwd_kmer = seq[i:i + KMER_LEN]
            if "N" not in fwd_kmer:
                index.setdefault(fwd_kmer, []).append((chrom, i, "+"))
            rev_pos = n - i - KMER_LEN
            rev_kmer = rc[i:i + KMER_LEN]
            if "N" not in rev_kmer:
                index.setdefault(rev_kmer, []).append((chrom, rev_pos, "-"))
    return index


# Keep the old name as an alias so existing code continues to work
def build_index(genome: dict[str, str]) -> dict[str, list[tuple[str, int, str]]]:
    """Return dict mapping 20-mer -> list of (chrom, pos, strand)."""
    index: dict[str, list[tuple[str, int, str]]] = {}
    for chrom, seq in genome.items():
        seq = seq.upper()
        rc = reverse_complement(seq)
        n = len(seq)
        for i in range(n - KMER_LEN + 1):
            fwd_kmer = seq[i:i + KMER_LEN]
            if "N" not in fwd_kmer:
                index.setdefault(fwd_kmer, []).append((chrom, i, "+"))
            rev_pos = n - i - KMER_LEN
            rev_kmer = rc[i:i + KMER_LEN]
            if "N" not in rev_kmer:
                index.setdefault(rev_kmer, []).append((chrom, rev_pos, "-"))
    return index


def load_fasta_genome(path: Path) -> dict[str, str]:
    """Parse FASTA file into dict of chrom->sequence."""
    genome: dict[str, str] = {}
    chrom, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if chrom is not None:
                    genome[chrom] = "".join(parts).upper()
                chrom = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if chrom is not None:
        genome[chrom] = "".join(parts).upper()
    return genome


def make_demo_genome(seed: int = 42) -> dict[str, str]:
    """Generate a synthetic 5000bp genome with planted guide targets."""
    rng = random.Random(seed)
    bases = "ACGT"

    guides = [
        "ACGTACGTACGTACGTACGT",
        "TGCATGCATGCATGCATGCA",
        "GCTAGCTAGCTAGCTAGCTA",
    ]

    seq = "".join(rng.choice(bases) for _ in range(5000))
    seq = list(seq)

    # Plant on-target NGG hits for each guide.
    # PAM check reads seq[pos+20:pos+22]; NGG means those 2 bases must be "GG".
    positions = [200, 1500, 3000]
    for guide, pos in zip(guides, positions):
        for j, base in enumerate(guide):
            seq[pos + j] = base
        # PAM: positions +20, +21 must be "GG" (the NG is implied by N being any)
        seq[pos + 20] = "G"
        seq[pos + 21] = "G"
        seq[pos + 22] = "A"

    # Plant off-targets with 1-2 mismatches
    off_positions = [500, 800, 1200, 2000, 2500, 4000]
    for i, pos in enumerate(off_positions):
        guide = guides[i % len(guides)]
        seq_slice = list(guide)
        n_mm = (i % 2) + 1
        mm_positions = rng.sample(range(20), n_mm)
        for mp in mm_positions:
            alt = [b for b in bases if b != seq_slice[mp]]
            seq_slice[mp] = rng.choice(alt)
        for j, base in enumerate(seq_slice):
            seq[pos + j] = base
        seq[pos + 20] = "G"
        seq[pos + 21] = "G"
        seq[pos + 22] = "C"

    return {"chr1": "".join(seq)}
