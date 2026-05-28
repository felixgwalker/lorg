import re
from typing import List, Dict, Any, Optional


# Feature type categories we care about
_GENE_TYPES = {"gene", "ncRNA_gene", "pseudogene"}
_EXON_TYPES = {"exon", "CDS"}
_REGULATORY_TYPES = {"promoter", "enhancer", "CTCF_binding_site", "insulator",
                     "transcription_factor_binding_site", "regulatory_region",
                     "DNase_I_hypersensitive_site", "open_chromatin_region"}


def annotate_locus(chrom: str, start: int, end: int,
                   gff_path: Optional[str] = None) -> List[Dict[str, Any]]:
    """Find GFF3 features overlapping the given locus interval.

    Looks for genes, exons/CDS, and regulatory elements.  Returns a list of
    GenomicFeature dicts, each with keys:
        feature_type, chrom, start, end, strand, name, attributes

    If gff_path is None or the file has no overlapping features, returns an
    empty list.
    """
    if not gff_path:
        return []

    features: List[Dict[str, Any]] = []
    try:
        with open(gff_path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 9:
                    continue
                feat_chrom = parts[0]
                if feat_chrom != chrom:
                    continue
                try:
                    feat_start = int(parts[3])
                    feat_end = int(parts[4])
                except ValueError:
                    continue
                # Overlap check (GFF3 coords are 1-based inclusive)
                if feat_end < start + 1 or feat_start > end:
                    continue
                feat_type = parts[2]
                strand = parts[6]
                raw_attrs = parts[8]

                # Determine broad category
                if feat_type in _GENE_TYPES:
                    category = "gene"
                elif feat_type in _EXON_TYPES:
                    category = "exon"
                elif feat_type in _REGULATORY_TYPES:
                    category = "regulatory"
                else:
                    category = feat_type

                # Parse Name/ID from attributes
                name = feat_type
                for token in raw_attrs.split(";"):
                    token = token.strip()
                    if token.startswith("Name=") or token.startswith("gene_name="):
                        name = token.split("=", 1)[1]
                        break
                    if token.startswith("ID=") and name == feat_type:
                        name = token.split("=", 1)[1]

                features.append({
                    "feature_type": category,
                    "chrom": feat_chrom,
                    "start": feat_start,
                    "end": feat_end,
                    "strand": strand,
                    "name": name,
                    "attributes": raw_attrs,
                })
    except OSError:
        return []

    return features


def parse_gff3(gff_path):
    genes = []
    with open(gff_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            if parts[2] == "gene":
                try:
                    genes.append({
                        "chrom": parts[0],
                        "start": int(parts[3]),
                        "end": int(parts[4]),
                        "strand": parts[6],
                    })
                except ValueError:
                    continue
    return genes


def gene_density_score(genes, chrom, window_start, window_end, flank=500000):
    region_start = window_start - flank
    region_end = window_end + flank
    count = 0
    for g in genes:
        if g["chrom"] == chrom and g["start"] < region_end and g["end"] > region_start:
            count += 1
    if count == 0:
        return 30
    elif count <= 2:
        return 20
    elif count <= 5:
        return 10
    else:
        return 0


def _scan_regulatory_motifs(seq):
    if not seq:
        return []
    positions = []
    for m in re.finditer(r"ATG", seq.upper()):
        positions.append(m.start())
    for m in re.finditer(r"TATAAA", seq.upper()):
        positions.append(m.start())
    return positions


def regulatory_proximity_score(sequences, chrom, window_start, window_end, flank=50000):
    seq = sequences.get(chrom, "")
    if not seq:
        return 15

    left_start = max(0, window_start - flank)
    left_seq = seq[left_start:window_start]
    right_end = min(len(seq), window_end + flank)
    right_seq = seq[window_end:right_end]

    left_hits = _scan_regulatory_motifs(left_seq)
    right_hits = _scan_regulatory_motifs(right_seq)

    min_dist = float("inf")
    for pos in left_hits:
        dist = (window_start - left_start) - pos
        if 0 < dist < min_dist:
            min_dist = dist
    for pos in right_hits:
        dist = pos
        if 0 < dist < min_dist:
            min_dist = dist

    if min_dist == float("inf") or min_dist > 50000:
        return 30
    elif min_dist > 10000:
        return 10
    else:
        return 0


def demo_genes():
    return [
        {"chrom": "chr1", "start": 10000, "end": 12000, "strand": "+"},
        {"chrom": "chr1", "start": 80000, "end": 85000, "strand": "-"},
        {"chrom": "chr1", "start": 90000, "end": 95000, "strand": "+"},
        {"chrom": "chr1", "start": 96000, "end": 98000, "strand": "+"},
        {"chrom": "chr1", "start": 99000, "end": 103000, "strand": "-"},
        {"chrom": "chr1", "start": 105000, "end": 110000, "strand": "+"},
        {"chrom": "chr2", "start": 20000, "end": 25000, "strand": "+"},
    ]
