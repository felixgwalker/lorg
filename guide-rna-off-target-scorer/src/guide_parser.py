"""Parse guide RNA sequences from FASTA or CSV, or from a raw sequence string."""

import csv
import re
import warnings
from pathlib import Path


def parse_guides(input_str) -> list[tuple[str, str]]:
    """Return list of (name, sequence) tuples.

    Accepts:
    - A Path object pointing to a FASTA (.fa / .fasta / .fna) or CSV file.
    - A str that is a valid file path to one of the above formats.
    - A raw nucleotide sequence string (20-nt, ACGTN).  It is assigned the
      name "guide_1".
    - A multi-guide raw string where guides are separated by commas or
      whitespace, each optionally prefixed with a name and '=' (e.g.
      "guide_1=ACGTACGTACGTACGTACGT,guide_2=TGCATGCATGCATGCATGCA").

    Warnings are printed (not raised) for guides containing N.
    """
    # Coerce to string for initial inspection
    input_as_str = str(input_str) if not isinstance(input_str, str) else input_str

    # Try to treat as a file path first
    candidate_path = Path(input_as_str)
    if candidate_path.exists() and candidate_path.is_file():
        suffix = candidate_path.suffix.lower()
        if suffix in (".fa", ".fasta", ".fna"):
            return _parse_fasta(candidate_path)
        elif suffix == ".csv":
            return _parse_csv(candidate_path)
        else:
            raise ValueError(f"Unsupported guide file format: {suffix}")

    # If a Path object was passed but didn't resolve, raise early
    if isinstance(input_str, Path):
        raise FileNotFoundError(f"Guide file not found: {input_str}")

    # Otherwise treat as raw sequence input
    return _parse_raw_string(input_as_str)


def _parse_raw_string(raw: str) -> list[tuple[str, str]]:
    """Parse a raw guide string — single sequence or comma/whitespace separated."""
    import re as _re
    raw = raw.strip()
    guides = []

    # Split on commas or whitespace runs
    tokens = _re.split(r"[,\s]+", raw)
    tokens = [t for t in tokens if t]

    for idx, token in enumerate(tokens, 1):
        if "=" in token:
            name, seq = token.split("=", 1)
        else:
            name = f"guide_{idx}"
            seq = token
        seq = seq.upper().strip()
        if not seq:
            continue
        if "N" in seq:
            warnings.warn(f"Guide {name!r} contains N bases: {seq}", UserWarning, stacklevel=3)
        guides.append((name, seq))

    if not guides:
        raise ValueError(f"No valid guide sequences found in input: {raw!r}")
    return guides


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
    """Return cleaned 20-nt protospacer or raise ValueError.

    Raises ValueError for wrong length or non-ACGTN bases.
    Warns (does not raise) if the sequence contains N.
    """
    seq = seq.upper().strip()
    if len(seq) != 20:
        raise ValueError(f"Guide must be 20 nt, got {len(seq)}: {seq}")
    valid = set("ACGTN")
    invalid = set(seq) - valid
    if invalid:
        raise ValueError(f"Invalid bases in guide (only ACGTN allowed): {invalid}")
    if "N" in seq:
        warnings.warn(f"Guide sequence contains N bases: {seq}", UserWarning, stacklevel=2)
    return seq


def make_demo_guides() -> list[tuple[str, str]]:
    return [
        ("guide_1", "ACGTACGTACGTACGTACGT"),
        ("guide_2", "TGCATGCATGCATGCATGCA"),
        ("guide_3", "GCTAGCTAGCTAGCTAGCTA"),
    ]
