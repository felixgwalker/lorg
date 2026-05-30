"""
CNV input parser — BED and VCF formats.

BED format (4+ columns, tab or comma separated):
    chrom  start  end  name  [score  strand  cnv_type  ...]
    Columns 1–3 are mandatory.  Column 4 is used as the CNV ID.
    Column 7, if present, is interpreted as the CNV type.

VCF format (≥ 4.1 with structural variants):
    SVTYPE in INFO distinguishes deletions (DEL), duplications (DUP),
    inversions (INV), and generic copy-number variants (CNV).
    END and SVLEN INFO tags are used to derive the end coordinate.
    Records without a resolvable extent are silently skipped.

Public helpers
--------------
parse_cnv_bed(path)   — parse a BED-style file (tab or comma) → list[dict]
make_demo_cnvs()      — return a list of synthetic CNV dicts for demo mode
parse_cnvs(path)      — auto-detect BED/VCF and return list[CNVRecord]
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


def parse_cnv_bed(path: Path) -> list[dict]:
    """
    Parse a BED-style CNV file and return a list of plain dicts.

    Accepts tab-separated or comma-separated files.  Header lines beginning
    with '#', 'track', or 'browser' are skipped automatically.

    Expected columns (0-indexed):
        0  chrom
        1  start   (0-based integer)
        2  end     (0-based integer, exclusive)
        3  name    → cnv_id                      (optional, default "cnv_<n>")
        4  score                                  (ignored)
        5  strand                                 (ignored)
        6  cnv_type (DEL/DUP/…)                  (optional)

    The returned dicts have keys:
        chrom, start, end, cnv_id, cnv_type, size_bp

    Args:
        path: Path to BED or CSV CNV file.

    Returns:
        List of CNV dicts, one per valid line.

    Raises:
        FileNotFoundError: If *path* does not exist.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"CNV BED file not found: {path}")

    records: list[dict] = []
    n_bad = 0

    with open(path, encoding="utf-8", errors="replace") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue

            # Detect delimiter: prefer tab, fall back to comma
            if "\t" in line:
                cols = line.split("\t")
            else:
                cols = line.split(",")

            if len(cols) < 3:
                logger.debug("BED line %d: fewer than 3 fields — skipping.", lineno)
                n_bad += 1
                continue

            chrom = cols[0].strip()
            try:
                start = int(cols[1].strip())
                end   = int(cols[2].strip())
            except ValueError:
                logger.debug("BED line %d: non-integer start/end — skipping.", lineno)
                n_bad += 1
                continue

            if end <= start:
                logger.debug("BED line %d: end ≤ start — skipping.", lineno)
                n_bad += 1
                continue

            cnv_id = cols[3].strip() if len(cols) > 3 and cols[3].strip() else f"cnv_{lineno}"
            cnv_type = "UNKNOWN"

            if len(cols) > 6 and cols[6].strip():
                cnv_type = cols[6].strip()
            if cnv_type == "UNKNOWN":
                cnv_type = _type_from_name(cnv_id)

            # Normalise aliases
            upper = cnv_type.upper()
            cnv_type = _TYPE_ALIASES.get(upper, upper)
            if cnv_type not in SUPPORTED_CNV_TYPES:
                cnv_type = "UNKNOWN"

            records.append({
                "chrom":    chrom,
                "start":    start,
                "end":      end,
                "cnv_id":   cnv_id,
                "cnv_type": cnv_type,
                "size_bp":  end - start,
            })

    if n_bad:
        logger.warning("parse_cnv_bed: skipped %d malformed lines in %s.", n_bad, path.name)
    logger.info("parse_cnv_bed: parsed %d CNVs from %s.", len(records), path.name)
    return records


def make_demo_cnvs() -> list[dict]:
    """
    Return a list of synthetic CNV dicts for demo / smoke-test mode.

    The demo set covers a range of sizes, types, and chromosomes so that
    every significance tier (LIKELY_BENIGN, VUS, LIKELY_PATHOGENIC) is
    represented after scoring.

    Returns:
        List of dicts with keys: chrom, start, end, cnv_id, cnv_type, size_bp.
    """
    raw = [
        # (chrom, start,        end,          cnv_id,             cnv_type)
        ("chr1",  1_000_000,    1_050_000,    "demo_DEL_small",   "DEL"),   # 50 kb   — benign size
        ("chr1",  5_000_000,    5_800_000,    "demo_DEL_BRCA2",   "DEL"),   # 800 kb  — overlaps BRCA2
        ("chr2",  25_000_000,   30_000_000,   "demo_DUP_large",   "DUP"),   # 5 Mb    — large DUP
        ("chr3",  120_000_000,  121_000_000,  "demo_DEL_TP53",    "DEL"),   # 1 Mb    — overlaps TP53
        ("chr4",  10_000_000,   10_020_000,   "demo_DUP_tiny",    "DUP"),   # 20 kb   — benign size
        ("chr5",  50_000_000,   50_600_000,   "demo_DEL_PTEN",    "DEL"),   # 600 kb  — overlaps PTEN
        ("chr6",  30_000_000,   30_500_000,   "demo_DUP_mid",     "DUP"),   # 500 kb  — mid-size DUP
        ("chr7",  117_000_000,  117_800_000,  "demo_DEL_CFTR",    "DEL"),   # 800 kb  — overlaps CFTR
        ("chr9",  21_000_000,   21_200_000,   "demo_DEL_CDKN2A",  "DEL"),   # 200 kb  — overlaps CDKN2A
        ("chr10", 89_000_000,   89_400_000,   "demo_DEL_PTEN_2",  "DEL"),   # 400 kb  — overlaps PTEN
        ("chr13", 32_000_000,   33_000_000,   "demo_DEL_BRCA2_2", "DEL"),   # 1 Mb    — overlaps BRCA2
        ("chr13", 48_000_000,   49_500_000,   "demo_DEL_RB1",     "DEL"),   # 1.5 Mb  — overlaps RB1
        ("chr17", 7_000_000,    7_700_000,    "demo_DEL_TP53_2",  "DEL"),   # 700 kb  — overlaps TP53
        ("chr17", 43_000_000,   44_200_000,   "demo_DEL_BRCA1",   "DEL"),   # 1.2 Mb  — overlaps BRCA1
        ("chr17", 29_000_000,   32_000_000,   "demo_DUP_NF1",     "DUP"),   # 3 Mb    — overlaps NF1
        ("chr22", 22_000_000,   22_500_000,   "demo_DEL_NF2",     "DEL"),   # 500 kb  — overlaps NF2
        ("chrX",  30_000_000,   30_300_000,   "demo_DUP_chrX",    "DUP"),   # 300 kb
        ("chr1",  200_000_000,  206_000_000,  "demo_DEL_6Mb",     "DEL"),   # 6 Mb    — very large
        ("chr8",  128_000_000,  128_030_000,  "demo_DEL_MYC",     "DEL"),   # 30 kb   — small near MYC
        ("chr11", 5_000_000,    5_200_000,    "demo_DEL_KCNQ1",   "DEL"),   # 200 kb
    ]
    records = []
    for chrom, start, end, cnv_id, cnv_type in raw:
        records.append({
            "chrom":    chrom,
            "start":    start,
            "end":      end,
            "cnv_id":   cnv_id,
            "cnv_type": cnv_type,
            "size_bp":  end - start,
        })
    logger.info("make_demo_cnvs: generated %d synthetic CNV records.", len(records))
    return records


def cnv_dicts_to_records(cnv_dicts: list[dict]) -> list[CNVRecord]:
    """Convert a list of CNV dicts (from parse_cnv_bed / make_demo_cnvs) to CNVRecord objects."""
    records = []
    for d in cnv_dicts:
        records.append(CNVRecord(
            cnv_id=d["cnv_id"],
            chrom=d["chrom"],
            start=d["start"],
            end=d["end"],
            cnv_type=d["cnv_type"],
        ))
    return sorted(records, key=lambda r: (r.chrom, r.start))


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
