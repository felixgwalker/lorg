"""Parse FASTA/CLUSTAL alignments into a list of (name, sequence) pairs."""

from __future__ import annotations

import random


def read_alignment(file_path: str) -> list[tuple[str, str]]:
    """Read a FASTA or CLUSTAL alignment file; return list of (name, seq) tuples."""
    with open(file_path) as fh:
        content = fh.read()
    if content.lstrip().startswith(">"):
        return _parse_fasta(content)
    if content.lstrip().startswith("CLUSTAL"):
        return _parse_clustal(content)
    return _parse_fasta(content)


def generate_demo_alignment(
    n_species: int = 5,
    seq_length: int = 500,
    seed: int = 42,
) -> list[tuple[str, str]]:
    """Generate synthetic alignment with known divergences for demo mode."""
    rng = random.Random(seed)
    bases = "ACGT"
    species_names = ["Human", "Chimpanzee", "Gorilla", "Orangutan", "Gibbon"][:n_species]
    divergences = [0.0, 0.01, 0.03, 0.08, 0.20][:n_species]

    ancestor = "".join(rng.choices(bases, k=seq_length))
    result: list[tuple[str, str]] = []
    for name, div in zip(species_names, divergences):
        seq = list(ancestor)
        for i in range(seq_length):
            if rng.random() < div:
                current = seq[i]
                others = [b for b in bases if b != current]
                seq[i] = rng.choice(others)
        result.append((name, "".join(seq)))
    return result


def validate_alignment(sequences: list[tuple[str, str]]) -> None:
    """Raise ValueError if sequences differ in length or fewer than 2 are present."""
    if len(sequences) < 2:
        raise ValueError(f"Alignment must contain at least 2 sequences; got {len(sequences)}.")
    lengths = {len(s) for _, s in sequences}
    if len(lengths) > 1:
        raise ValueError(
            f"All sequences must be the same length. Found lengths: {sorted(lengths)}."
        )


def _parse_fasta(content: str) -> list[tuple[str, str]]:
    sequences: list[tuple[str, str]] = []
    name = ""
    parts: list[str] = []
    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name:
                sequences.append((name, "".join(parts).upper()))
            name = line[1:].split()[0]
            parts = []
        else:
            parts.append(line.replace(" ", "").replace("-", "-"))
    if name:
        sequences.append((name, "".join(parts).upper()))
    return sequences


def _parse_clustal(content: str) -> list[tuple[str, str]]:
    seqs: dict[str, list[str]] = {}
    order: list[str] = []
    for line in content.splitlines():
        if not line.strip() or line.startswith("CLUSTAL") or line.startswith(" ") or line.startswith("*"):
            continue
        parts = line.split()
        if len(parts) >= 2 and not parts[0].startswith("*"):
            name = parts[0]
            seq_part = parts[1]
            if name not in seqs:
                seqs[name] = []
                order.append(name)
            seqs[name].append(seq_part.upper())
    return [(n, "".join(seqs[n])) for n in order]
