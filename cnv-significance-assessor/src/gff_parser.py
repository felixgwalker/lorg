"""
GFF3 annotation parser.

Extracts gene bodies and regulatory elements for CNV interval intersection.

GFF3 column layout (1-based, closed coordinates):
    seqname  source  feature  start  end  score  strand  frame  attributes

All coordinates are converted to 0-based half-open [start, end) on output
to match the BED/VCF convention used throughout the pipeline.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class GeneRecord:
    """A single gene body or regulatory element from a GFF3 file."""

    chrom: str
    start: int          # 0-based inclusive
    end: int            # 0-based exclusive  [start, end)
    gene_id: str
    gene_name: str
    feature_type: str   # "gene" | "regulatory"
    strand: str         # "+" | "-" | "."


def parse_gff3(
    path: Path,
    gene_features: set[str] | None = None,
    regulatory_features: set[str] | None = None,
) -> list[GeneRecord]:
    """
    Parse GFF3 and return gene and regulatory element records.

    Only features whose type appears in *gene_features* or
    *regulatory_features* are retained.  All others are silently skipped.

    Args:
        path:                Path to GFF3 file.
        gene_features:       GFF3 feature types treated as gene bodies.
                             Defaults to config.GFF3_GENE_FEATURES.
        regulatory_features: GFF3 feature types treated as regulatory elements.
                             Defaults to config.GFF3_REGULATORY_FEATURES.

    Returns:
        List of GeneRecord sorted by (chrom, start).

    Raises:
        FileNotFoundError: If the path does not exist.
        ValueError:        If no records matching the wanted features are found.
    """
    from src.config import GFF3_GENE_FEATURES, GFF3_REGULATORY_FEATURES

    if gene_features is None:
        gene_features = GFF3_GENE_FEATURES
    if regulatory_features is None:
        regulatory_features = GFF3_REGULATORY_FEATURES

    wanted = gene_features | regulatory_features

    if not path.exists():
        raise FileNotFoundError(f"GFF3 annotation file not found: {path}")

    records: list[GeneRecord] = []
    n_skipped = 0
    _warned_attrs: set[str] = set()

    with open(path, encoding="utf-8", errors="replace") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) < 9:
                logger.debug("GFF3 line %d: fewer than 9 columns, skipping.", lineno)
                n_skipped += 1
                continue

            feature = cols[2]
            if feature not in wanted:
                continue

            chrom = cols[0]
            try:
                gff_start = int(cols[3])    # 1-based inclusive
                gff_end   = int(cols[4])    # 1-based inclusive
            except ValueError:
                logger.debug("GFF3 line %d: non-integer coordinates, skipping.", lineno)
                n_skipped += 1
                continue

            if gff_end < gff_start:
                logger.debug("GFF3 line %d: end < start, skipping.", lineno)
                n_skipped += 1
                continue

            # GFF3 [1-based, closed] → [0-based, half-open)
            start = gff_start - 1
            end   = gff_end

            strand = cols[6] if cols[6] in {"+", "-"} else "."
            attrs  = _parse_attributes(cols[8])

            gene_id   = (attrs.get("ID")
                         or attrs.get("gene_id")
                         or f"feature_{lineno}")
            gene_name = (attrs.get("Name")
                         or attrs.get("gene_name")
                         or attrs.get("gene")
                         or gene_id)

            # Warn once if neither ID nor Name is found (likely unusual GFF3 dialect)
            if "ID" not in attrs and "Name" not in attrs and chrom not in _warned_attrs:
                logger.warning(
                    "GFF3: no 'ID' or 'Name' attribute on line %d (feature '%s'). "
                    "Using positional fallback IDs.", lineno, feature,
                )
                _warned_attrs.add(chrom)

            ftype = "gene" if feature in gene_features else "regulatory"

            records.append(GeneRecord(
                chrom=chrom,
                start=start,
                end=end,
                gene_id=gene_id,
                gene_name=gene_name,
                feature_type=ftype,
                strand=strand,
            ))

    if n_skipped:
        logger.debug("GFF3: skipped %d lines (malformed or wrong column count).", n_skipped)

    records.sort(key=lambda r: (r.chrom, r.start))

    n_gene = sum(1 for r in records if r.feature_type == "gene")
    n_reg  = sum(1 for r in records if r.feature_type == "regulatory")
    logger.info(
        "Parsed %d annotation records from GFF3 (%d gene, %d regulatory): %s",
        len(records), n_gene, n_reg, path.name,
    )
    return records


def _parse_attributes(attr_str: str) -> dict[str, str]:
    """
    Parse a GFF3 attribute string into a key → value dict.

    Handles:
      - Standard 'key=value' pairs separated by ';'
      - Bare flags (no '=') are stored with empty-string value
      - Percent-decoded values (basic URL-decoding for %2C, %3B, %3D)
    """
    result: dict[str, str] = {}
    for token in attr_str.split(";"):
        token = token.strip()
        if not token:
            continue
        if "=" in token:
            k, _, v = token.partition("=")
            result[k.strip()] = _url_decode(v.strip())
        else:
            result[token] = ""
    return result


def _url_decode(s: str) -> str:
    """Minimal percent-decoding for the characters most common in GFF3 attributes."""
    return (s
            .replace("%2C", ",")
            .replace("%3B", ";")
            .replace("%3D", "=")
            .replace("%25", "%"))
