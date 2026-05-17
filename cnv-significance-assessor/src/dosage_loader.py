"""
Haploinsufficiency and triplosensitivity score loader.

Reads a CSV file and returns a mapping from gene name to available dosage
sensitivity scores.  Standard column names from gnomAD v4, ClinGen, and
the Collins et al. (2022) resource are recognised automatically via
case-insensitive alias matching.

Recognised score columns
------------------------
pLI      — probability of loss-of-function intolerance (gnomAD)
pHaplo   — probability of haploinsufficiency (Collins et al.)
pTriplo  — probability of triplosensitivity (Collins et al.)

Recognised gene name columns
-----------------------------
gene, gene_name, symbol, hgnc_symbol, genename, gene_id
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# ── Column alias tables (all lower-cased) ─────────────────────────────────
_GENE_ALIASES: set[str] = {
    "gene", "gene_name", "symbol", "hgnc_symbol", "genename", "gene_id", "name",
}
_PLI_ALIASES: set[str] = {
    "pli", "p_li", "p(li)", "gnomad_pli", "lof_pli",
}
_PHAPLO_ALIASES: set[str] = {
    "phaplo", "p_haplo", "haplo_score", "haploinsufficiency_score",
    "hi_score", "hi", "p_haploinsufficiency",
}
_PTRIPLO_ALIASES: set[str] = {
    "ptriplo", "p_triplo", "triplo_score", "triplosensitivity_score",
    "ts_score", "ts", "p_triplosensitivity",
}


def load_dosage_scores(path: Path) -> dict[str, dict[str, float]]:
    """
    Load dosage sensitivity scores from a CSV file.

    Args:
        path: Path to a CSV with a gene name column and any of pLI, pHaplo,
              pTriplo columns.  Column names are matched case-insensitively.

    Returns:
        Dict mapping upper-cased gene name → score dict.
        The score dict contains only the keys that were present and parseable
        in the input (keys: "pLI", "pHaplo", "pTriplo").

        Example::

            {
                "BRCA1": {"pHaplo": 0.97, "pTriplo": 0.12},
                "TP53":  {"pLI": 0.99, "pHaplo": 0.95},
            }

    Raises:
        FileNotFoundError: If the path does not exist.
        ValueError:        If no gene name column can be identified.
    """
    if not path.exists():
        raise FileNotFoundError(f"Dosage score file not found: {path}")

    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError(f"Dosage score file has no header row: {path}")

        # Map lower-cased header → original header string
        headers_lower: dict[str, str] = {
            h.lower().strip(): h for h in reader.fieldnames
        }

        gene_col    = _match_col(headers_lower, _GENE_ALIASES)
        pli_col     = _match_col(headers_lower, _PLI_ALIASES)
        phaplo_col  = _match_col(headers_lower, _PHAPLO_ALIASES)
        ptriplo_col = _match_col(headers_lower, _PTRIPLO_ALIASES)

        if gene_col is None:
            recognised = ", ".join(sorted(_GENE_ALIASES))
            raise ValueError(
                f"Cannot identify gene name column in {path.name}.  "
                f"Expected one of (case-insensitive): {recognised}."
            )

        score_cols: list[tuple[str, str | None]] = [
            ("pLI",    pli_col),
            ("pHaplo", phaplo_col),
            ("pTriplo", ptriplo_col),
        ]
        present_cols = [k for k, col in score_cols if col is not None]
        if not present_cols:
            logger.warning(
                "Dosage score file %s: none of pLI/pHaplo/pTriplo columns found. "
                "Dosage scores will be unavailable.",
                path.name,
            )

        result: dict[str, dict[str, float]] = {}
        n_rows = 0

        for row in reader:
            gene = row[gene_col].strip()
            if not gene:
                continue

            scores: dict[str, float] = {}
            for key, col in score_cols:
                if col is not None:
                    val = _parse_float(row.get(col, ""))
                    if val is not None:
                        scores[key] = val

            result[gene.upper()] = scores
            n_rows += 1

    logger.info(
        "Loaded dosage scores for %d genes from %s "
        "(columns present: %s).",
        n_rows, path.name, ", ".join(present_cols) if present_cols else "none",
    )
    return result


# ── Helpers ───────────────────────────────────────────────────────────────

def _match_col(headers_lower: dict[str, str], aliases: set[str]) -> str | None:
    """Return the original column name matching any alias, or None."""
    for alias in aliases:
        if alias in headers_lower:
            return headers_lower[alias]
    return None


def _parse_float(s: str) -> float | None:
    """Convert a string to float; return None for missing / non-numeric values."""
    s = s.strip()
    if not s or s.upper() in {".", "NA", "N/A", "NAN", "NULL", "NONE", ""}:
        return None
    try:
        return float(s)
    except ValueError:
        return None
