"""Signal 4 — Conjugative element homology (BLASTX against protein DB)."""

import logging
import tempfile
from pathlib import Path
from typing import Optional

from src.aggregator import SIGNAL_WEIGHTS
from src.blast import db_exists, run_blast, write_query_fasta
from src.models import BlastHit, HostProfile, QuerySequence, SignalResult

logger = logging.getLogger(__name__)

WEIGHT = SIGNAL_WEIGHTS["conjugative"]

# Key protein markers for conjugative mobility.  Hits against these markers
# receive a 1.5× score multiplier because they indicate the actual machinery
# of transfer, not just co-localisation.
_HIGH_CONCERN = {
    "relaxase", "mob", "virb4", "vird4", "t4ss",
    "trai", "traa", "nikb", "moba", "trag",
    "conjugat", "orit", "mobiliz",
}


def _marker_boost(hit: BlastHit) -> float:
    desc = hit.subject_description.lower()
    for marker in _HIGH_CONCERN:
        if marker in desc:
            return 1.5
    return 1.0


def _score(hits: list[BlastHit], cap: float = 3.0) -> float:
    if not hits:
        return 0.0
    total = sum(
        min((h.pct_identity / 100.0) * (h.query_coverage / 100.0) * _marker_boost(h), 1.0)
        for h in hits
    )
    return min(total, cap) / cap


def run(
    query: QuerySequence,
    host: HostProfile,
    data_dir: Path,
    threads: int = 4,
    blast_bin_dir: Optional[Path] = None,
    **kwargs,
) -> SignalResult:
    db_path = data_dir / "conjugative" / "conjugative_proteins"

    if not db_exists(db_path, dbtype="prot"):
        return SignalResult(
            signal_name="conjugative",
            score=None,
            weight=WEIGHT,
            evidence={},
            warning=(
                "Conjugative element protein DB not found. "
                "Run: python data/download_databases.py --conjugative --entrez-email you@example.com"
            ),
            skipped=True,
        )

    with tempfile.TemporaryDirectory() as tmp:
        query_fasta = write_query_fasta(query.sequence, query.identifier, Path(tmp))
        hits = run_blast(
            query_fasta, db_path, "blastx",
            evalue=1e-5, threads=threads,
            blast_bin_dir=blast_bin_dir,
            query_length=query.length,
        )

    score = _score(hits)
    high_concern_hits = [h for h in hits if _marker_boost(h) > 1.0]
    markers_found = sorted({
        m for h in high_concern_hits
        for m in _HIGH_CONCERN
        if m in h.subject_description.lower()
    })

    return SignalResult(
        signal_name="conjugative",
        score=score,
        weight=WEIGHT,
        evidence={
            "hit_count":               len(hits),
            "high_concern_hit_count":  len(high_concern_hits),
            "high_concern_markers_found": markers_found,
            "top_hits": [
                {
                    "subject":         h.subject_id,
                    "description":     h.subject_description[:80],
                    "pct_identity":    h.pct_identity,
                    "query_coverage":  h.query_coverage,
                    "evalue":          h.evalue,
                    "high_concern":    _marker_boost(h) > 1.0,
                }
                for h in sorted(hits, key=lambda h: h.evalue)[:5]
            ],
        },
    )
