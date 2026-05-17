"""Signal 1 — IS element proximity (ISfinder BLAST)."""

import logging
import tempfile
from pathlib import Path
from typing import Optional

from src.aggregator import SIGNAL_WEIGHTS
from src.blast import blast_score_from_hits, db_exists, run_blast, write_query_fasta
from src.models import HostProfile, QuerySequence, SignalResult

logger = logging.getLogger(__name__)

WEIGHT = SIGNAL_WEIGHTS["is_proximity"]


def run(
    query: QuerySequence,
    host: HostProfile,
    data_dir: Path,
    threads: int = 4,
    blast_bin_dir: Optional[Path] = None,
    **kwargs,
) -> SignalResult:
    db_path = data_dir / "isfinder" / "isfinder"

    if not db_exists(db_path, dbtype="nucl"):
        return SignalResult(
            signal_name="is_proximity",
            score=None,
            weight=WEIGHT,
            evidence={},
            warning=(
                "ISfinder BLAST DB not found. "
                "Run: python data/download_databases.py --isfinder"
            ),
            skipped=True,
        )

    with tempfile.TemporaryDirectory() as tmp:
        query_fasta = write_query_fasta(query.sequence, query.identifier, Path(tmp))
        hits = run_blast(
            query_fasta, db_path, "blastn",
            evalue=1e-5, threads=threads,
            blast_bin_dir=blast_bin_dir,
            query_length=query.length,
        )

    score = blast_score_from_hits(hits)

    return SignalResult(
        signal_name="is_proximity",
        score=score,
        weight=WEIGHT,
        evidence={
            "hit_count": len(hits),
            "top_hits": [
                {
                    "subject":       h.subject_id,
                    "description":   h.subject_description[:80],
                    "pct_identity":  h.pct_identity,
                    "query_coverage": h.query_coverage,
                    "evalue":        h.evalue,
                }
                for h in sorted(hits, key=lambda h: h.evalue)[:5]
            ],
        },
    )
