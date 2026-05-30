"""Parse VCF for synonymous variants, or generate synthetic demo variants."""

from __future__ import annotations

import random
from dataclasses import dataclass, field

from src.codon_tables import HUMAN_CODON_FREQ


@dataclass
class SynonymousVariant:
    """A single synonymous variant with codon context."""
    variant_id: str
    chrom: str
    pos: int
    ref_allele: str
    alt_allele: str
    ref_codon: str
    alt_codon: str
    exon_pos: int
    exon_length: int
    distance_to_donor: int
    distance_to_acceptor: int
    gene: str = ""
    transcript: str = ""
    extra: dict = field(default_factory=dict)


def parse_vcf(vcf_path: str) -> list[SynonymousVariant]:
    """Parse a VCF file and extract synonymous variant records."""
    variants: list[SynonymousVariant] = []
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos_str, vid, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
            info = fields[7]
            info_dict = _parse_info(info)
            if "CSQ" not in info_dict and "ANN" not in info_dict:
                if info_dict.get("SVTYPE"):
                    continue
                ref_codon, alt_codon = _infer_codons(ref, alt, info_dict)
                if ref_codon is None:
                    continue
            else:
                ref_codon, alt_codon = _extract_csq_codons(info_dict)
                if ref_codon is None:
                    continue

            pos = int(pos_str)
            exon_length = int(info_dict.get("EXON_LEN", 300))
            exon_pos = int(info_dict.get("EXON_POS", pos % exon_length))
            v = SynonymousVariant(
                variant_id=vid if vid != "." else f"{chrom}:{pos}:{ref}>{alt}",
                chrom=chrom,
                pos=pos,
                ref_allele=ref,
                alt_allele=alt,
                ref_codon=ref_codon.upper(),
                alt_codon=alt_codon.upper(),
                exon_pos=exon_pos,
                exon_length=exon_length,
                distance_to_donor=exon_length - exon_pos,
                distance_to_acceptor=exon_pos,
                gene=info_dict.get("GENE", ""),
                transcript=info_dict.get("TRANSCRIPT", ""),
            )
            variants.append(v)
    return variants


def generate_demo_variants(n: int = 10, seed: int = 42) -> list[SynonymousVariant]:
    """Generate synthetic synonymous variants for demo mode."""
    rng = random.Random(seed)
    syn_pairs = _build_synonymous_pairs()
    chroms = [f"chr{i}" for i in range(1, 6)]
    variants = []
    genes = ["BRCA1", "TP53", "CFTR", "LDLR", "APOE", "MYH7", "PKD1", "NF1", "RYR1", "TSC2"]
    for i in range(n):
        ref_codon, alt_codon = rng.choice(syn_pairs)
        chrom = rng.choice(chroms)
        pos = rng.randint(1_000_000, 50_000_000)
        exon_length = rng.randint(120, 600)
        exon_pos = rng.randint(3, exon_length - 3)
        v = SynonymousVariant(
            variant_id=f"DEMO_{i+1:03d}",
            chrom=chrom,
            pos=pos,
            ref_allele=ref_codon[2],
            alt_allele=alt_codon[2],
            ref_codon=ref_codon,
            alt_codon=alt_codon,
            exon_pos=exon_pos,
            exon_length=exon_length,
            distance_to_donor=exon_length - exon_pos,
            distance_to_acceptor=exon_pos,
            gene=genes[i % len(genes)],
            transcript=f"NM_{100000 + i:06d}.1",
        )
        variants.append(v)
    return variants


def _build_synonymous_pairs() -> list[tuple[str, str]]:
    """Build all synonymous codon pairs from the standard genetic code.

    Uses the locally-defined STANDARD_GENETIC_CODE to avoid any external
    dependency (e.g. BioPython).  Stop codons are excluded.
    """
    from src.codon_tables import STANDARD_GENETIC_CODE
    aa_to_codons: dict[str, list[str]] = {}
    for codon, aa in STANDARD_GENETIC_CODE.items():
        if aa == "*":
            continue  # skip stop codons
        aa_to_codons.setdefault(aa, []).append(codon)
    pairs: list[tuple[str, str]] = []
    for codons in aa_to_codons.values():
        if len(codons) < 2:
            continue
        for i, c1 in enumerate(codons):
            for c2 in codons[i + 1:]:
                if c1 != c2:
                    pairs.append((c1, c2))
    return pairs


def _parse_info(info: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for token in info.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            result[k] = v
        else:
            result[token] = "1"
    return result


def _infer_codons(ref: str, alt: str, info: dict) -> tuple[str | None, str | None]:
    rc = info.get("REF_CODON")
    ac = info.get("ALT_CODON")
    if rc and ac and len(rc) == 3 and len(ac) == 3:
        return rc, ac
    return None, None


def _extract_csq_codons(info: dict) -> tuple[str | None, str | None]:
    csq = info.get("CSQ", info.get("ANN", ""))
    for field in csq.split(","):
        parts = field.split("|")
        for j, p in enumerate(parts):
            if len(p) == 7 and "/" in p:
                halves = p.split("/")
                if len(halves) == 2 and len(halves[0]) == 3 and len(halves[1]) == 3:
                    return halves[0], halves[1]
    return None, None
