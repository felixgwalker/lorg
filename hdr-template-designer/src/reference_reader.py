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


def read_locus(fasta_path: Path, chrom: str, start: int, end: int) -> str:
    """Load a reference sequence slice [start, end) from a FASTA file.

    Parameters
    ----------
    fasta_path : Path
        Path to the FASTA file.
    chrom : str
        Chromosome/contig name to extract from.
    start : int
        0-based start coordinate (inclusive).
    end : int
        0-based end coordinate (exclusive).

    Returns
    -------
    str
        Uppercase sequence slice.

    Raises
    ------
    KeyError
        If *chrom* is not present in the FASTA file.
    ValueError
        If *start* >= *end* or coordinates are negative.
    """
    if start < 0 or end < 0:
        raise ValueError(f"Coordinates must be non-negative (got start={start}, end={end})")
    if start >= end:
        raise ValueError(f"start ({start}) must be less than end ({end})")

    genome = load_fasta(fasta_path)
    if chrom not in genome:
        available = list(genome.keys())[:5]
        raise KeyError(
            f"Chromosome '{chrom}' not found in {fasta_path}. "
            f"Available (first 5): {available}"
        )

    return extract_region(genome[chrom], start, end)


def make_demo_locus(seed: int = 42) -> tuple[str, int]:
    """Return a synthetic 1000 bp reference sequence with a cut site at position 500.

    The sequence is constructed to have approximately 45 % GC content.  A 20 nt
    guide target is embedded immediately upstream of the cut site followed by an
    NGG PAM, so PAM-disruption logic has something to work with.

    Returns
    -------
    tuple[str, int]
        ``(sequence, cut_pos)`` where *sequence* is 1000 bp and *cut_pos* is 500.
    """
    rng = random.Random(seed)

    # Build a 1000 bp sequence biased toward ~45 % GC.
    # Weights: A=0.275, C=0.225, G=0.225, T=0.275 → GC = 0.45
    bases = "ACGT"
    weights = [0.275, 0.225, 0.225, 0.275]
    seq = list("".join(rng.choices(bases, weights=weights, k=1000)))

    cut_pos = 500

    # Embed a 20 nt guide target ending at cut_pos (positions 480-499).
    guide = "GCTACGATCGATCGATCGAT"  # ~45% GC guide sequence
    for i, base in enumerate(guide):
        seq[480 + i] = base

    # Place NGG PAM at positions 500-502 (immediately after cut).
    # The first base of the PAM region doubles as the edit target (C -> T).
    seq[cut_pos] = "C"      # ref allele; will be mutated to T in demo edit
    seq[cut_pos + 1] = "G"  # PAM G1
    seq[cut_pos + 2] = "G"  # PAM G2

    return "".join(seq).upper(), cut_pos


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
