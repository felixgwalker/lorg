"""Read FASTA reference and extract regions.

Salvaged from hdr-template-designer (deleted stage1f).
"""

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


def read_locus(fasta_path: Path, chrom: str, start: int, end: int) -> str:
    """Load a reference sequence slice [start, end) from a FASTA file."""
    if start < 0 or end < 0:
        raise ValueError(f"Coordinates must be non-negative (got start={start}, end={end})")
    if start >= end:
        raise ValueError(f"start ({start}) must be less than end ({end})")
    genome = load_fasta(fasta_path)
    if chrom not in genome:
        raise KeyError(f"Chromosome '{chrom}' not found in {fasta_path}")
    seq = genome[chrom]
    return seq[max(0, start):min(len(seq), end)]
