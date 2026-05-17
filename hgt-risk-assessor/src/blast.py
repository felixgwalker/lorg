"""BLAST+ subprocess wrapper and shared hit-scoring utilities."""

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from src.models import BlastHit

logger = logging.getLogger(__name__)

# Tabular output format — 8 columns
BLAST_OUTFMT = "6 qseqid sseqid pident length qcovs evalue bitscore stitle"


def db_exists(db_path: Path, dbtype: str = "nucl") -> bool:
    """Return True if a BLAST database index file is present."""
    ext = ".nhr" if dbtype == "nucl" else ".phr"
    return db_path.with_suffix(ext).exists()


def write_query_fasta(sequence: str, identifier: str, tmp_dir: Path) -> Path:
    fasta_path = tmp_dir / "query.fasta"
    fasta_path.write_text(f">{identifier}\n{sequence}\n", encoding="utf-8")
    return fasta_path


def run_blast(
    query_fasta: Path,
    db_path: Path,
    blast_type: str,            # "blastn" or "blastx"
    evalue: float = 1e-5,
    threads: int = 4,
    blast_bin_dir: Optional[Path] = None,
    query_length: Optional[int] = None,
) -> list[BlastHit]:
    """Run BLAST and return parsed hits. Returns [] on any failure."""
    binary = blast_type
    if blast_bin_dir:
        binary = str(blast_bin_dir / blast_type)

    cmd = [
        binary,
        "-query", query_fasta.as_posix(),
        "-db", db_path.as_posix(),
        "-outfmt", BLAST_OUTFMT,
        "-evalue", str(evalue),
        "-num_threads", str(threads),
    ]

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    except FileNotFoundError:
        logger.warning(
            f"{blast_type} not found. Install BLAST+ and ensure it is on PATH, "
            "or use --blast-bin-dir."
        )
        return []
    except subprocess.TimeoutExpired:
        logger.warning(f"{blast_type} timed out after 300 s.")
        return []

    if proc.returncode != 0:
        logger.warning(f"{blast_type} exited {proc.returncode}: {proc.stderr[:300]}")
        return []

    return _parse_tabular(proc.stdout)


def _parse_tabular(output: str) -> list[BlastHit]:
    hits: list[BlastHit] = []
    for line in output.strip().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t", maxsplit=7)
        if len(parts) < 7:
            continue
        try:
            hits.append(BlastHit(
                query_id=parts[0],
                subject_id=parts[1],
                pct_identity=float(parts[2]),
                alignment_length=int(parts[3]),
                query_coverage=float(parts[4]),
                evalue=float(parts[5]),
                bit_score=float(parts[6]),
                subject_description=parts[7] if len(parts) > 7 else "",
            ))
        except (ValueError, IndexError) as exc:
            logger.debug(f"Skipping malformed BLAST line: {exc}")
    return hits


def blast_score_from_hits(hits: list[BlastHit], cap: float = 3.0) -> float:
    """
    Aggregate BLAST hits into a 0.0–1.0 score.

    Each hit contributes (pct_identity/100) * (query_coverage/100).
    The sum is capped at `cap` then normalised.  A cap of 3 means three
    full-identity, full-coverage hits saturate the score at 1.0.
    """
    if not hits:
        return 0.0
    total = sum(
        (h.pct_identity / 100.0) * (h.query_coverage / 100.0)
        for h in hits
    )
    return min(total, cap) / cap
