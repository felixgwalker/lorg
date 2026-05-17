"""Signal 3 — Integron association (attC regex scan + INTEGRALL BLAST)."""

import logging
import re
import tempfile
from pathlib import Path
from typing import Optional

from src.aggregator import SIGNAL_WEIGHTS
from src.blast import blast_score_from_hits, db_exists, run_blast, write_query_fasta
from src.models import AttCSite, HostProfile, QuerySequence, SignalResult

logger = logging.getLogger(__name__)

WEIGHT = SIGNAL_WEIGHTS["integron"]

# attC bottom-strand consensus (simplified).
# Conserved left end:  GTT[Y][R][R]  (Y = C/T, R = A/G)
# Variable spacer:     50–200 bp
# Conserved right end: [R][R]AAC
# Real attC sites: 57–141 bp spacer.  This regex is intentionally permissive
# for screening; all hits are labelled "putative" in the report.
_ATTC = re.compile(r"GTT[CT][AG][AG].{50,200}[AG]{2}AAC", re.IGNORECASE)
_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _find_attc_sites(sequence: str) -> list[AttCSite]:
    sites: list[AttCSite] = []
    seq_len = len(sequence)
    rc = _reverse_complement(sequence)

    for strand, seq in (("+", sequence), ("-", rc)):
        for m in _ATTC.finditer(seq):
            matched = m.group()
            spacer_len = len(matched) - 12      # ~6 conserved bp each end
            if strand == "+":
                start, end = m.start(), m.end()
            else:
                start = seq_len - m.end()
                end = seq_len - m.start()
            sites.append(AttCSite(
                start=start,
                end=end,
                sequence=(matched[:20] + "...") if len(matched) > 20 else matched,
                spacer_length=max(spacer_len, 0),
                strand=strand,
            ))
    return sites


def run(
    query: QuerySequence,
    host: HostProfile,
    data_dir: Path,
    threads: int = 4,
    blast_bin_dir: Optional[Path] = None,
    **kwargs,
) -> SignalResult:
    # Phase 1: regex attC scan (always runs, no DB required)
    attc_sites = _find_attc_sites(query.sequence)
    attc_sub = min(len(attc_sites) / 3.0, 1.0)

    # Phase 2: INTEGRALL BLAST (runs only if DB present)
    db_path = data_dir / "integrall" / "integrall"
    blast_hits = []
    blast_available = db_exists(db_path, dbtype="nucl")

    if blast_available:
        with tempfile.TemporaryDirectory() as tmp:
            query_fasta = write_query_fasta(query.sequence, query.identifier, Path(tmp))
            blast_hits = run_blast(
                query_fasta, db_path, "blastn",
                evalue=1e-5, threads=threads,
                blast_bin_dir=blast_bin_dir,
                query_length=query.length,
            )
        blast_sub = blast_score_from_hits(blast_hits)
        score = 0.5 * attc_sub + 0.5 * blast_sub
    else:
        logger.info("INTEGRALL BLAST DB not found; integron signal uses attC regex only.")
        score = attc_sub

    warning = (
        "" if blast_available
        else "INTEGRALL BLAST DB absent — attC regex scan only. "
             "Putative hits may include false positives."
    )

    return SignalResult(
        signal_name="integron",
        score=score,
        weight=WEIGHT,
        evidence={
            "attc_sites_found": len(attc_sites),
            "attc_sites": [
                {
                    "start":          s.start,
                    "end":            s.end,
                    "strand":         s.strand,
                    "spacer_length":  s.spacer_length,
                    "sequence_prefix": s.sequence,
                }
                for s in attc_sites[:10]
            ],
            "attc_regex_only": not blast_available,
            "integrall_blast_hits": len(blast_hits),
            "top_blast_hits": [
                {
                    "subject":        h.subject_id,
                    "description":    h.subject_description[:80],
                    "pct_identity":   h.pct_identity,
                    "query_coverage": h.query_coverage,
                    "evalue":         h.evalue,
                }
                for h in sorted(blast_hits, key=lambda h: h.evalue)[:5]
            ],
        },
        warning=warning,
    )
