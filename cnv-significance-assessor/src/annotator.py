"""
CNV annotator — interval intersection and population frequency lookup.

For each CNV the annotator:
  1. Finds overlapping gene bodies and regulatory elements from GFF3.
  2. Retrieves the best-available dosage sensitivity score for each gene.
  3. Cross-references against the population CNV database (if provided).

Interval intersection is implemented with a chromsome-keyed index and a
sweep over sorted intervals; no external binary dependency is required.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

from src.cnv_parser import CNVRecord, parse_cnvs
from src.gff_parser import GeneRecord

logger = logging.getLogger(__name__)


@dataclass
class AnnotatedCNV:
    """A CNV enriched with gene overlap and population frequency data."""

    record: CNVRecord

    # Gene and regulatory overlap
    overlapping_genes: list[GeneRecord]       = field(default_factory=list)
    overlapping_regulatory: list[GeneRecord]  = field(default_factory=list)

    # Dosage sensitivity (max across overlapping genes)
    max_dosage_sensitivity: float = 0.0
    dosage_metric: str = "none"   # "pHaplo" | "pTriplo" | "pLI" | "none"

    # Population frequency
    pop_frequency: float | None = None    # None = not in database
    pop_match_count: int = 0              # number of overlapping population CNVs

    @property
    def n_genes(self) -> int:
        return len(self.overlapping_genes)

    @property
    def gene_names(self) -> list[str]:
        return [g.gene_name for g in self.overlapping_genes]


# ── Main annotation function ──────────────────────────────────────────────

def annotate_cnvs(
    cnvs: list[CNVRecord],
    genes: list[GeneRecord],
    dosage_scores: dict[str, dict[str, float]],
    pop_cnvs: list[CNVRecord] | None,
    overlap_fraction: float = 0.1,
) -> list[AnnotatedCNV]:
    """
    Annotate a list of CNVs with gene overlap and population frequency.

    Args:
        cnvs:             CNV records to annotate.
        genes:            Gene and regulatory records from GFF3.
        dosage_scores:    Gene-name → score dict from dosage_loader.
        pop_cnvs:         Population-frequency CNV database records, or None.
        overlap_fraction: Minimum fraction of the smaller interval that must
                          overlap for a gene to be considered "overlapping".

    Returns:
        List of AnnotatedCNV in the same order as *cnvs*.
    """
    gene_index       = _build_index(genes)
    pop_index        = _build_index(pop_cnvs) if pop_cnvs else {}

    results: list[AnnotatedCNV] = []
    for cnv in cnvs:
        ann = _annotate_one(cnv, gene_index, dosage_scores, pop_index, overlap_fraction)
        results.append(ann)

    n_with_genes = sum(1 for a in results if a.n_genes > 0)
    n_with_pop   = sum(1 for a in results if a.pop_frequency is not None)
    logger.info(
        "Annotation complete: %d/%d CNVs overlap ≥1 gene, "
        "%d/%d have population frequency data.",
        n_with_genes, len(results), n_with_pop, len(results),
    )
    return results


def load_population_cnvs(
    path: Path,
    min_size: int = 0,
) -> list[CNVRecord]:
    """
    Load population CNV database records (BED or VCF).

    Args:
        path:     Path to population CNV file.
        min_size: Minimum size filter (bp).

    Returns:
        Sorted list of CNVRecord.
    """
    logger.info("Loading population CNV database: %s", path)
    records = parse_cnvs(path, min_size=min_size)
    logger.info("  %d population CNV records loaded.", len(records))
    return records


# ── Internal helpers ──────────────────────────────────────────────────────

def _build_index(records: list) -> dict[str, list]:
    """
    Build a chromosome → sorted record list index.

    Works for any object with .chrom, .start, .end attributes.
    """
    index: dict[str, list] = {}
    for rec in records:
        index.setdefault(rec.chrom, []).append(rec)
    # Records within each chromosome are already sorted (both CNV and gene lists
    # come pre-sorted from their parsers), but we sort again to be defensive.
    for chrom in index:
        index[chrom].sort(key=lambda r: r.start)
    return index


def _annotate_one(
    cnv: CNVRecord,
    gene_index: dict[str, list[GeneRecord]],
    dosage_scores: dict[str, dict[str, float]],
    pop_index: dict[str, list[CNVRecord]],
    overlap_fraction: float,
) -> AnnotatedCNV:
    """Build an AnnotatedCNV for a single CNV record."""
    ann = AnnotatedCNV(record=cnv)

    # ── Gene and regulatory overlap ──────────────────────────────────────
    candidates = gene_index.get(cnv.chrom, [])
    for gene in _overlapping(cnv, candidates, overlap_fraction):
        if gene.feature_type == "gene":
            ann.overlapping_genes.append(gene)
        else:
            ann.overlapping_regulatory.append(gene)

    # ── Dosage sensitivity ───────────────────────────────────────────────
    is_del = cnv.cnv_type in {"DEL", "LOSS"}
    is_dup = cnv.cnv_type in {"DUP", "GAIN"}

    best_score = 0.0
    best_metric = "none"

    for gene in ann.overlapping_genes:
        scores = dosage_scores.get(gene.gene_name.upper(), {})

        if is_del:
            # Prefer pHaplo, fall back to pLI
            score = scores.get("pHaplo") or scores.get("pLI") or 0.0
            metric = "pHaplo" if "pHaplo" in scores else ("pLI" if "pLI" in scores else "none")
        elif is_dup:
            score  = scores.get("pTriplo", 0.0)
            metric = "pTriplo" if "pTriplo" in scores else "none"
        else:
            # INV / CNV / UNKNOWN: use max of any available score
            score  = max(scores.values(), default=0.0)
            metric = max(scores, key=lambda k: scores[k], default="none") if scores else "none"

        if score > best_score:
            best_score  = score
            best_metric = metric

    ann.max_dosage_sensitivity = best_score
    ann.dosage_metric = best_metric

    # ── Population frequency ─────────────────────────────────────────────
    pop_candidates = pop_index.get(cnv.chrom, [])
    if pop_candidates:
        ann.pop_frequency, ann.pop_match_count = _population_frequency(
            cnv, pop_candidates
        )

    return ann


def _overlapping(
    cnv: CNVRecord,
    candidates: list[GeneRecord],
    overlap_fraction: float,
) -> list[GeneRecord]:
    """
    Return genes that overlap the CNV by at least *overlap_fraction* of the
    shorter of the two intervals (reciprocal fraction logic).

    Candidates are assumed sorted by start; we stop once candidate.start
    exceeds cnv.end.
    """
    result: list[GeneRecord] = []
    for gene in candidates:
        if gene.start >= cnv.end:
            break       # sorted — no further candidates can overlap
        if gene.end <= cnv.start:
            continue    # gene ends before CNV starts

        overlap_len = min(cnv.end, gene.end) - max(cnv.start, gene.start)
        if overlap_len <= 0:
            continue

        shorter = min(cnv.end - cnv.start, gene.end - gene.start)
        if shorter > 0 and (overlap_len / shorter) >= overlap_fraction:
            result.append(gene)

    return result


def _population_frequency(
    cnv: CNVRecord,
    pop_records: list[CNVRecord],
) -> tuple[float | None, int]:
    """
    Estimate population frequency for a CNV from a frequency database.

    A population record is counted as "matching" when:
      - Same chromosome (already guaranteed by the caller's index lookup)
      - Reciprocal overlap ≥ 50 % (standard DGV/gnomAD-SV overlap criterion)
      - Same broad CNV type (DEL↔DEL, DUP↔DUP; INV/CNV match any)

    The returned frequency is the maximum AF INFO value observed in matching
    records, or 1/total_pop_size if no AF field is available (treated as
    "observed at least once" ≈ rare).

    When a population record has a name containing 'AF=' (VCF-derived), that
    value is extracted.  Otherwise, we return None for the frequency (≥1
    matching record means the variant is observed in the population but the
    exact frequency is unknown).

    Returns:
        Tuple of (frequency_or_None, match_count).
        frequency_or_None is None when no matching records exist.
    """
    _RECIP_OVERLAP = 0.5
    _FALLBACK_AF   = 1e-4     # treat "observed once, AF unknown" as 0.01%

    same_type_groups: dict[str, set[str]] = {
        "DEL":  {"DEL", "LOSS"},
        "DUP":  {"DUP", "GAIN"},
        "INV":  {"INV"},
        "CNV":  {"DEL", "DUP", "INV", "CNV", "LOSS", "GAIN"},
        "UNKNOWN": {"DEL", "DUP", "INV", "CNV", "LOSS", "GAIN", "UNKNOWN"},
    }
    allowed_types = same_type_groups.get(cnv.cnv_type, set())

    match_count = 0
    max_freq: float | None = None

    for pop in pop_records:
        if pop.start >= cnv.end:
            break
        if pop.end <= cnv.start:
            continue

        if pop.cnv_type not in allowed_types:
            continue

        overlap_len = min(cnv.end, pop.end) - max(cnv.start, pop.start)
        if overlap_len <= 0:
            continue

        cnv_len = cnv.end - cnv.start
        pop_len = pop.end - pop.start
        if cnv_len <= 0 or pop_len <= 0:
            continue

        recip = overlap_len / max(cnv_len, pop_len)
        if recip < _RECIP_OVERLAP:
            continue

        # Matching record found
        match_count += 1

        # Try to extract an AF value from the population record's ID field
        # (a common trick when converting DGV/gnomAD to BED — AF is embedded)
        af = _extract_af(pop.cnv_id)
        if af is None:
            af = _FALLBACK_AF

        if max_freq is None or af > max_freq:
            max_freq = af

    return max_freq, match_count


def _extract_af(name: str) -> float | None:
    """
    Extract an allele frequency from a CNV record name/ID string.

    Handles patterns like:
      - 'AF=0.0123'
      - 'gnomAD_DUP_AF_0.05'
      - plain float strings ('0.0042')
    """
    import re

    # Explicit AF= tag
    m = re.search(r"AF[=_]([0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)", name, re.IGNORECASE)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            pass

    # Name is itself a plain float
    try:
        val = float(name)
        if 0.0 <= val <= 1.0:
            return val
    except ValueError:
        pass

    return None
