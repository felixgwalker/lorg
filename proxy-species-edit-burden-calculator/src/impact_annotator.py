"""Impact annotation for proxy-species edit-burden variants.

Classifies each variant as one of:
    coding_synonymous        — SNP in CDS with same amino-acid (heuristic)
    coding_nonsynonymous     — SNP in CDS with amino-acid change
    coding_nonsense          — SNP in CDS creating a stop codon (heuristic)
    splice_site              — within 2 bp of an exon boundary
    regulatory               — within 1 kb upstream of a TSS (gene start)
    intronic                 — inside a gene but not in CDS/exon
    intergenic               — outside all annotated genes

If no GFF3 is supplied, all variants are classified as intergenic.
"""

from __future__ import annotations

IMPACT_SCORES = {
    "coding_nonsense": 12,
    "coding_nonsynonymous": 10,
    "splice_site": 8,
    "coding_synonymous": 3,
    "regulatory": 4,
    "intronic": 2,
    "intergenic": 1,
}

# Legacy names kept for backward compatibility with existing pipeline code.
IMPACT_SCORES["exonic_nonsynonymous"] = 10
IMPACT_SCORES["exonic_synonymous"] = 3


# ---------------------------------------------------------------------------
# Genetic-code helpers for synonymous / nonsense heuristics
# ---------------------------------------------------------------------------

_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def _snp_coding_impact(ref_base: str, alt_base: str, codon_pos: int,
                       context: str) -> str:
    """Return 'coding_synonymous', 'coding_nonsynonymous', or 'coding_nonsense'.

    Uses a simple heuristic: extract a codon from *context* (±10 bp window),
    replace the appropriate base, and compare amino acids.  Falls back to
    nonsynonymous if the context is too short to form a complete codon.
    """
    # codon_pos: 0, 1, or 2 — position within the codon.
    # Attempt to read the codon from context centred on the variant.
    ctx = context.upper()
    mid = len(ctx) // 2  # approximate position of the SNP
    # Snap to codon frame start.
    codon_start = mid - codon_pos
    if codon_start < 0 or codon_start + 3 > len(ctx):
        return "coding_nonsynonymous"
    ref_codon = ctx[codon_start: codon_start + 3]
    alt_codon = list(ref_codon)
    alt_codon[codon_pos] = alt_base.upper()
    alt_codon_str = "".join(alt_codon)
    ref_aa = _CODON_TABLE.get(ref_codon, "?")
    alt_aa = _CODON_TABLE.get(alt_codon_str, "?")
    if alt_aa == "*":
        return "coding_nonsense"
    if ref_aa == alt_aa:
        return "coding_synonymous"
    return "coding_nonsynonymous"


# ---------------------------------------------------------------------------
# GFF3 parser
# ---------------------------------------------------------------------------

def _parse_gff3(gff_path: str) -> dict:
    """Return {chrom: {"genes": [...], "exons": [...], "cds": [...]}}."""
    by_chrom: dict = {}

    def _ensure(chrom: str) -> dict:
        if chrom not in by_chrom:
            by_chrom[chrom] = {"genes": [], "exons": [], "cds": []}
        return by_chrom[chrom]

    with open(gff_path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom, _, ftype, start_s, end_s, _, strand, _, attrs = parts[:9]
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                continue

            rec = {"chrom": chrom, "start": start, "end": end,
                   "strand": strand, "attrs": attrs}

            bucket = _ensure(chrom)
            ftype_lower = ftype.lower()
            if ftype_lower in ("gene", "mrna", "transcript"):
                bucket["genes"].append(rec)
            elif ftype_lower == "exon":
                bucket["exons"].append(rec)
            elif ftype_lower == "cds":
                bucket["cds"].append(rec)

    return by_chrom


# ---------------------------------------------------------------------------
# Overlap helpers
# ---------------------------------------------------------------------------

def _overlaps(pos: int, feat: dict, margin: int = 0) -> bool:
    return (feat["start"] - margin) <= pos <= (feat["end"] + margin)


def _find_overlapping(pos: int, features: list, margin: int = 0):
    """Return first feature overlapping *pos* ± *margin*, or None."""
    for feat in features:
        if _overlaps(pos, feat, margin):
            return feat
    return None


def _is_splice_site(pos: int, exons: list, margin: int = 2) -> bool:
    """True if *pos* is within *margin* bp of any exon boundary but not inside."""
    for exon in exons:
        near_start = abs(pos - exon["start"]) <= margin
        near_end = abs(pos - exon["end"]) <= margin
        inside = exon["start"] <= pos <= exon["end"]
        if (near_start or near_end) and not inside:
            return True
    return False


def _is_regulatory(pos: int, genes: list, upstream_bp: int = 1000) -> bool:
    """True if *pos* is within *upstream_bp* bp upstream of a gene TSS."""
    for gene in genes:
        if gene["strand"] == "-":
            tss = gene["end"]
            if tss < pos <= tss + upstream_bp:
                return True
        else:
            tss = gene["start"]
            if tss - upstream_bp <= pos < tss:
                return True
    return False


# ---------------------------------------------------------------------------
# Public API: annotate_impact()
# ---------------------------------------------------------------------------

_SPLICE_MARGIN = 2       # bp from exon boundary → splice_site
_UPSTREAM_BP = 1_000     # bp upstream of TSS → regulatory


def annotate_impact(variants, gff_path=None):
    """Annotate each variant dict (or Variant dataclass) with impact category.

    Parameters
    ----------
    variants : list
        List of variant dicts (from align_genomes) or Variant dataclass instances
        from classify_variants().
    gff_path : str or None
        Path to a GFF3 annotation file.  If None or the file cannot be parsed,
        all variants are classified as intergenic.

    Returns
    -------
    list
        The same list, with 'impact_category' and 'impact_score' set on each
        element (dict key or dataclass attribute).
    """
    by_chrom: dict = {}
    if gff_path:
        try:
            by_chrom = _parse_gff3(gff_path)
        except Exception:
            by_chrom = {}

    annotated = []
    for v in variants:
        # Support both dict-style (legacy pipeline) and Variant dataclass.
        _is_dc = hasattr(v, "__dataclass_fields__")

        if _is_dc:
            pos = v.position
            chrom = v.chrom
            vtype = v.type          # "SNP", "SMALL_INS", "SMALL_DEL", etc.
            ref_b = v.ref_allele
            alt_b = v.alt_allele
            ctx = v.context_seq
        else:
            pos = v.get("pos", 0)
            chrom = v.get("chrom", "")
            vtype = v.get("type", "SNV")
            ref_b = v.get("ref", "N") or "N"
            alt_b = v.get("alt", "N") or "N"
            ctx = ""

        chrom_data = by_chrom.get(chrom, {})
        genes = chrom_data.get("genes", [])
        exons = chrom_data.get("exons", [])
        cds_feats = chrom_data.get("cds", [])

        if not by_chrom:
            # No annotation available — use lightweight heuristics.
            vc = (v.variant_class if _is_dc else v.get("variant_class", ""))
            if vc.startswith("SV"):
                category = "coding_nonsynonymous"
            elif vc.startswith("LARGE"):
                category = "intronic"
            else:
                category = "intergenic"
        elif _is_splice_site(pos, exons, _SPLICE_MARGIN):
            category = "splice_site"
        elif _find_overlapping(pos, cds_feats):
            # Variant is inside a CDS.
            if vtype in ("SNP", "SNV") and len(ref_b) == 1 and len(alt_b) == 1:
                # Determine codon position from distance to CDS start.
                cds_feat = _find_overlapping(pos, cds_feats)
                codon_pos = (pos - cds_feat["start"]) % 3
                category = _snp_coding_impact(ref_b, alt_b, codon_pos, ctx)
            else:
                # Indel in CDS — frameshift if not divisible by 3.
                ref_l = len(ref_b) if ref_b else 0
                alt_l = len(alt_b) if alt_b else 0
                delta = abs(ref_l - alt_l)
                if delta % 3 != 0:
                    category = "coding_nonsynonymous"  # frameshift
                else:
                    category = "coding_nonsynonymous"  # in-frame indel
        elif _find_overlapping(pos, exons):
            # Exon but not annotated as CDS (e.g. UTR).
            category = "coding_nonsynonymous"
        elif _find_overlapping(pos, genes):
            # Inside a gene but not in an exon or CDS.
            category = "intronic"
        elif _is_regulatory(pos, genes, _UPSTREAM_BP):
            category = "regulatory"
        else:
            category = "intergenic"

        score = IMPACT_SCORES.get(category, 1)

        if _is_dc:
            v.impact_category = category
            v.impact_score = score
        else:
            v["impact_category"] = category
            v["impact_score"] = score

        annotated.append(v)

    return annotated
