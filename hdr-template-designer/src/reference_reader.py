"""Read FASTA reference and extract regions."""

import random
from pathlib import Path


def load_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    chrom, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if chrom is not None:
                    seqs[chrom] = "".join(parts).upper()
                chrom = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if chrom is not None:
        seqs[chrom] = "".join(parts).upper()
    return seqs


def extract_region(seq: str, start: int, end: int) -> str:
    """Extract subsequence [start, end) with boundary clamping."""
    start = max(0, start)
    end = min(len(seq), end)
    return seq[start:end]


def make_demo_reference(seed: int = 42) -> tuple[str, str]:
    """Return (chrom_name, 500bp synthetic sequence)."""
    rng = random.Random(seed)
    bases = "ACGT"
    seq = "".join(rng.choice(bases) for _ in range(500))
    seq = list(seq)
    # Plant a clear cut site at position 250 with NGG PAM
    # Guide target: 20nt at 230-250, PAM at 250-252
    guide = "AGTCAGTCAGTCAGTCAGTC"
    for i, b in enumerate(guide):
        seq[230 + i] = b
    seq[250] = "C"  # will be mutated to T (edit)
    seq[251] = "G"
    seq[252] = "G"
    return "chr1", "".join(seq).upper()
