"""Signal 5 — Prophage context (PHASTER API with local BLAST fallback)."""

import logging
import tempfile
import time
from pathlib import Path
from typing import Optional

import requests

from src.aggregator import SIGNAL_WEIGHTS
from src.blast import blast_score_from_hits, db_exists, run_blast, write_query_fasta
from src.models import HostProfile, PhasterRegion, PhasterResult, QuerySequence, SignalResult

logger = logging.getLogger(__name__)

WEIGHT = SIGNAL_WEIGHTS["prophage"]

PHASTER_URL = "https://phaster.ca/phaster_api"
_POLL_INTERVAL = 30         # seconds between status polls
_MAX_POLLS = 40             # 20 minutes maximum wall time
_REQUEST_TIMEOUT = 15       # per-request timeout in seconds

# Completeness → contribution to score (3 intact prophages → score 1.0)
_COMPLETENESS_WEIGHT = {"intact": 1.0, "questionable": 0.5, "incomplete": 0.25}


# ---------------------------------------------------------------------------
# PHASTER API helpers
# ---------------------------------------------------------------------------

def _submit(sequence: str, session: requests.Session) -> Optional[str]:
    """POST sequence to PHASTER.  Returns job_id or None."""
    try:
        resp = session.post(
            PHASTER_URL,
            data={"fasta_data": f">query\n{sequence}"},
            timeout=_REQUEST_TIMEOUT,
        )
        resp.raise_for_status()
        job_id = resp.json().get("job_id")
        if job_id:
            logger.info(f"PHASTER job submitted: {job_id}")
        return job_id
    except Exception as exc:
        logger.warning(f"PHASTER submission failed: {exc}")
        return None


def _parse_response(data: dict) -> PhasterResult:
    regions: list[PhasterRegion] = []
    for r in data.get("regions", []):
        try:
            regions.append(PhasterRegion(
                region_id=str(r.get("region", "")),
                completeness=r.get("completeness", "incomplete").lower(),
                start=int(r.get("region_start", 0)),
                end=int(r.get("region_end", 0)),
                gc_content=float(r.get("gc_percentage", 0.0)) / 100.0,
                num_cds=int(r.get("num_cds", 0)),
            ))
        except (ValueError, TypeError) as exc:
            logger.debug(f"Skipping malformed PHASTER region: {exc}")
    return PhasterResult(status="complete", regions=regions, raw_response=data)


def _poll(job_id: str, session: requests.Session) -> Optional[PhasterResult]:
    """Poll until the job is complete.  Returns PhasterResult or None."""
    for attempt in range(1, _MAX_POLLS + 1):
        time.sleep(_POLL_INTERVAL)
        try:
            resp = session.get(PHASTER_URL, params={"job_id": job_id},
                               timeout=_REQUEST_TIMEOUT)
            if resp.status_code == 503:
                logger.warning(f"PHASTER 503 on poll {attempt}; retrying...")
                continue
            resp.raise_for_status()
            data = resp.json()
            status = data.get("status", "").lower()
            logger.info(f"PHASTER poll {attempt}/{_MAX_POLLS}: {status}")
            if status in ("complete", "done"):
                return _parse_response(data)
            if "error" in status:
                logger.warning(f"PHASTER error: {data}")
                return None
        except Exception as exc:
            logger.warning(f"PHASTER poll {attempt} failed: {exc}")
    logger.warning("PHASTER timed out.")
    return None


def _score_phaster(result: PhasterResult) -> float:
    if not result.regions:
        return 0.0
    total = sum(_COMPLETENESS_WEIGHT.get(r.completeness, 0.1) for r in result.regions)
    return min(total / 3.0, 1.0)


# ---------------------------------------------------------------------------
# Signal runner
# ---------------------------------------------------------------------------

def run(
    query: QuerySequence,
    host: HostProfile,
    data_dir: Path,
    no_network: bool = False,
    threads: int = 4,
    blast_bin_dir: Optional[Path] = None,
    **kwargs,
) -> SignalResult:
    # --- Try PHASTER API ---
    phaster_result: Optional[PhasterResult] = None
    if not no_network:
        try:
            with requests.Session() as session:
                job_id = _submit(query.sequence, session)
                if job_id:
                    phaster_result = _poll(job_id, session)
        except Exception as exc:
            logger.warning(f"PHASTER API error: {exc}")

    if phaster_result is not None:
        return SignalResult(
            signal_name="prophage",
            score=_score_phaster(phaster_result),
            weight=WEIGHT,
            evidence={
                "source":       "phaster_api",
                "region_count": len(phaster_result.regions),
                "regions": [
                    {
                        "id":           r.region_id,
                        "completeness": r.completeness,
                        "start":        r.start,
                        "end":          r.end,
                        "gc_pct":       round(r.gc_content * 100, 1),
                        "num_cds":      r.num_cds,
                    }
                    for r in phaster_result.regions
                ],
            },
        )

    # --- Fallback: local BLAST ---
    db_path = data_dir / "phage" / "phage_genes"
    if db_exists(db_path, dbtype="nucl"):
        with tempfile.TemporaryDirectory() as tmp:
            query_fasta = write_query_fasta(query.sequence, query.identifier, Path(tmp))
            hits = run_blast(
                query_fasta, db_path, "blastn",
                evalue=1e-5, threads=threads,
                blast_bin_dir=blast_bin_dir,
                query_length=query.length,
            )
        return SignalResult(
            signal_name="prophage",
            score=blast_score_from_hits(hits),
            weight=WEIGHT,
            evidence={
                "source":    "local_blast",
                "hit_count": len(hits),
                "top_hits": [
                    {
                        "subject":        h.subject_id,
                        "description":    h.subject_description[:80],
                        "pct_identity":   h.pct_identity,
                        "query_coverage": h.query_coverage,
                        "evalue":         h.evalue,
                    }
                    for h in sorted(hits, key=lambda h: h.evalue)[:5]
                ],
            },
            warning="PHASTER API unavailable; using local phage gene BLAST DB.",
        )

    # --- Both unavailable ---
    return SignalResult(
        signal_name="prophage",
        score=None,
        weight=WEIGHT,
        evidence={},
        warning=(
            "PHASTER API unreachable and no local phage gene DB found. "
            "Run: python data/download_databases.py --phage  (or remove --no-network)"
        ),
        skipped=True,
    )
