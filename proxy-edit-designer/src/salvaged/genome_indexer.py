"""Build k-mer index of 20-mers from genome FASTA.

Salvaged from guide-rna-off-target-scorer (deleted stage1f).
"""

import random
from pathlib import Path

KMER_LEN = 20
_COMP_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP_TABLE)[::-1]


def build_genome_index(genome: dict[str, str]) -> dict[str, list[tuple[str, int, str]]]:
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
