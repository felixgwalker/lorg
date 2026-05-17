"""Parse guide RNA sequences from FASTA or CSV."""

import csv
from pathlib import Path


def parse_guides(path: Path) -> list[tuple[str, str]]:
    """Return list of (name, sequence) tuples from FASTA or CSV."""
    suffix = path.suffix.lower()
    if suffix in (".fa", ".fasta", ".fna"):
        return _parse_fasta(path)
    elif suffix == ".csv":
        return _parse_csv(path)
    else:
        raise ValueError(f"Unsupported guide file format: {suffix}")


def _parse_fasta(path: Path) -> list[tuple[str, str]]:
    guides = []
    name, seq_parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    guides.append((name, "".join(seq_parts).upper()))
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if name is not None:
        guides.append((name, "".join(seq_parts).upper()))
    return guides


def _parse_csv(path: Path) -> list[tuple[str, str]]:
    guides = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            name = row.get("name") or row.get("id") or row.get("Name") or row.get("ID", "guide")
            seq = (row.get("sequence") or row.get("seq") or row.get("Sequence") or "").strip().upper()
            if seq:
                guides.append((name, seq))
    return guides


def validate_guide(seq: str) -> str:
    """Return cleaned 20nt protospacer or raise."""
    seq = seq.upper().strip()
    if len(seq) != 20:
        raise ValueError(f"Guide must be 20 nt, got {len(seq)}: {seq}")
    valid = set("ACGT")
    invalid = set(seq) - valid
    if invalid:
        raise ValueError(f"Invalid bases in guide: {invalid}")
    return seq


def make_demo_guides() -> list[tuple[str, str]]:
    return [
        ("guide_1", "ACGTACGTACGTACGTACGT"),
        ("guide_2", "TGCATGCATGCATGCATGCA"),
        ("guide_3", "GCTAGCTAGCTAGCTAGCTA"),
    ]
