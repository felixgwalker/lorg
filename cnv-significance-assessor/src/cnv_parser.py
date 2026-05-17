"""
CNV input parser — BED and VCF formats.

BED format (4+ columns, tab or space separated):
    chrom  start  end  name  [score  strand  cnv_type  ...]
    Columns 1–3 are mandatory.  Column 4 is used as the CNV ID.
    Column 7, if present, is interpreted as the CNV type.

VCF format (≥ 4.1 with structural variants):
    SVTYPE in INFO distinguishes deletions (DEL), duplications (DUP),
    inversions (INV), and generic copy-number variants (CNV).
    END and SVLEN INFO tags are used to derive the end coordinate.
    Records without a resolvable extent are silently skipped.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)

# CNV types recognised by the classifier
SUPPORTED_CNV_TYPES: frozenset[str] = frozenset({
    "DEL", "DUP", "INV", "CNV", "GAIN", "LOSS", "COMPLEX", "UNKNOWN",
})

# Canonical normalisations applied before SUPPORTED_CNV_TYPES lookup
_TYPE_ALIASES: dict[str, str] = {
    "DELETION":     "DEL",
    "DUPLICATION":  "DUP",
    "INVERSION":    "INV",
    "GAIN":         "DUP",
    "LOSS":         "DEL",
}


@dataclass
class CNVRecord:
    """A single CNV call in 0-based half-open coordinates."""

    cnv_id: str
    chrom: str
    start: int          # 0-based inclusive
    end: int            # 0-based exclusive  [start, end)
    cnv_type: str       # DEL | DUP | INV | CNV | COMPLEX | UNKNOWN
    size: int = field(init=False)

    def __post_init__(self) -> None:
        self.size = self.end - self.start
        upper = self.cnv_type.upper()
        self.cnv_type = _TYPE_ALIASES.get(upper, upper)
        if self.cnv_type not in SUPPORTED_CNV_TYPES:
            self.cnv_type = "UNKNOWN"


def parse_cnvs(path: Path, min_size: int = 0) -> list[CNVRecord]:
    """
    Auto-detect BED or VCF format and return filtered CNV records.

    Detection order:
      1. Extension (.vcf / .vcf.gz) → VCF
      2. First non-comment line starts with '##fileformat=VCF' → VCF
      3. Otherwise → BED

    Args:
        path:     Path to input file.
        min_size: Skip CNVs strictly smaller than this (bp).

    Returns:
        List of CNVRecord sorted by (chrom, start).

    Raises:
        ValueError: If the file is completely empty or unreadable.
        FileNotFoundError: If the path does not exist.
    """
    if not path.exists():
        raise FileNotFoundError(f"CNV input file not found: {path}")

    suffix = path.suffix.lower()
    if suffix in {".vcf", ".gz"} or _peek_vcf(path):
        records = _parse_vcf(path)
    else:
        records = _parse_bed(path)

    n_raw = len(records)
    if min_size > 0:
        records = [r for r in records if r.size >= min_size]
        n_dropped = n_raw - len(records)
        if n_dropped:
            logger.info("  Size filter (<  %d bp): removed %d CNVs.", min_size, n_dropped)

    logger.info("  %d CNVs retained after size filter.", len(records))
    return sorted(records, key=lambda r: (r.chrom, r.start))


# ── Format detection ──────────────────────────────────────────────────────

def _peek_vcf(path: Path) -> bool:
    """Return True if the first non-empty line looks like a VCF header."""
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    return line.startswith("##fileformat=VCF") or line.startswith("##VCF")
    except OSError:
        pass
    return False


# ── BED parser ────────────────────────────────────────────────────────────

def _parse_bed(path: Path) -> list[CNVRecord]:
    """
    Parse CNV calls from a BED-like file.

    Column layout (1-indexed):
        1  chrom
        2  start  (0-based)
        3  end    (0-based, exclusive)
        4  name   → CNV ID                       (optional)
        5  score                                  (ignored)
        6  strand                                 (ignored)
        7  cnv_type (DEL/DUP/INV/…)              (optional)

    Lines beginning with '#', 'track', or 'browser' are skipped.
    """
    records: list[CNVRecord] = []
    n_bad = 0

    with open(path, encoding="utf-8", errors="replace") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue

            cols = line.split("\t")
            if len(cols) < 3:
                cols = line.split()
            if len(cols) < 3:
                logger.debug("BED line %d: fewer than 3 fields, skipping.", lineno)
                n_bad += 1
                continue

            chrom = cols[0]
            try:
                start = int(cols[1])
                end   = int(cols[2])
            except ValueError:
                logger.debug("BED line %d: non-integer start/end, skipping.", lineno)
                n_bad += 1
                continue

            if end <= start:
                logger.debug("BED line %d: end ≤ start, skipping.", lineno)
                n_bad += 1
                continue

            cnv_id   = cols[3].strip() if len(cols) > 3 else f"cnv_{lineno}"
            cnv_type = "UNKNOWN"

            if len(cols) > 6:
                cnv_type = cols[6].strip()
            if cnv_type == "UNKNOWN":
                cnv_type = _type_from_name(cnv_id)

            records.append(CNVRecord(
                cnv_id=cnv_id,
                chrom=chrom,
                start=start,
                end=end,
                cnv_type=cnv_type,
            ))

    if n_bad:
        logger.warning("BED: skipped %d malformed lines.", n_bad)
    logger.info("Parsed %d CNVs from BED: %s", len(records), path.name)
    return records


# ── VCF parser ────────────────────────────────────────────────────────────

def _parse_vcf(path: Path) -> list[CNVRecord]:
    """
    Parse structural variant records from a VCF file.

    Records are included when SVTYPE ∈ SUPPORTED_CNV_TYPES.  The end
    coordinate is derived from the END INFO tag; if absent, from
    abs(SVLEN) + POS.  Records with no resolvable extent are skipped.

    VCF POS is 1-based; converted to 0-based on output.
    """
    records: list[CNVRecord] = []
    n_skipped = 0

    with open(path, encoding="utf-8", errors="replace") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) < 8:
                logger.debug("VCF line %d: fewer than 8 fields, skipping.", lineno)
                n_skipped += 1
                continue

            chrom  = cols[0]
            try:
                pos = int(cols[1])      # 1-based
            except ValueError:
                n_skipped += 1
                continue

            vcf_id = cols[2] if cols[2] != "." else f"sv_{lineno}"
            info   = _parse_info(cols[7])

            svtype = _TYPE_ALIASES.get(info.get("SVTYPE", "").upper(),
                                       info.get("SVTYPE", "UNKNOWN").upper())
            if svtype not in SUPPORTED_CNV_TYPES:
                n_skipped += 1
                continue

            end = _resolve_end(info, pos)
            if end is None:
                logger.debug(
                    "VCF line %d (%s): cannot resolve END; skipping.", lineno, vcf_id
                )
                n_skipped += 1
                continue

            start = pos - 1     # convert 1-based → 0-based

            records.append(CNVRecord(
                cnv_id=vcf_id,
                chrom=chrom,
                start=start,
                end=end,
                cnv_type=svtype,
            ))

    if n_skipped:
        logger.info("VCF: skipped %d records (missing SVTYPE/END or unsupported type).", n_skipped)
    logger.info("Parsed %d CNVs from VCF: %s", len(records), path.name)
    return records


def _resolve_end(info: dict[str, str], pos: int) -> int | None:
    """Derive end coordinate from INFO END or SVLEN + POS."""
    if "END" in info:
        try:
            return int(info["END"])
        except ValueError:
            pass
    if "SVLEN" in info:
        try:
            return pos + abs(int(info["SVLEN"]))
        except ValueError:
            pass
    return None


def _parse_info(info_str: str) -> dict[str, str]:
    """Parse VCF INFO field → {key: value}; flag-only keys get value 'True'."""
    result: dict[str, str] = {}
    if info_str in {".", ""}:
        return result
    for token in info_str.split(";"):
        if "=" in token:
            k, _, v = token.partition("=")
            result[k] = v
        else:
            result[token] = "True"
    return result


# ── Helpers ───────────────────────────────────────────────────────────────

def _type_from_name(name: str) -> str:
    """Infer CNV type from a free-text name string."""
    upper = name.upper()
    for keyword, ctype in [
        ("DELETION", "DEL"), ("DEL", "DEL"),
        ("DUPLICATION", "DUP"), ("DUP", "DUP"), ("GAIN", "DUP"),
        ("LOSS", "DEL"),
        ("INVERSION", "INV"), ("INV", "INV"),
        ("CNV", "CNV"),
    ]:
        if keyword in upper:
            return ctype
    return "UNKNOWN"
