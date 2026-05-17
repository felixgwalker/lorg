"""Parse FASTA / VCF input and resolve host organism GC content."""

import logging
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

from src.models import InputFormat, QuerySequence, HostProfile

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Sequence input
# ---------------------------------------------------------------------------

def parse_fasta(path: Path) -> QuerySequence:
    """Parse the first record from a FASTA file."""
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError(f"No sequences found in {path}")
    if len(records) > 1:
        logger.warning(
            f"{path.name} contains {len(records)} records; using only the first "
            f"({records[0].id}). Run the pipeline once per sequence for multi-FASTA files."
        )
    rec = records[0]
    seq_str = str(rec.seq).upper()
    return QuerySequence(
        sequence=seq_str,
        identifier=rec.id,
        length=len(seq_str),
        gc_content=gc_fraction(rec.seq),
        source_format=InputFormat.FASTA,
    )


def parse_vcf(vcf_path: Path, ref_path: Path) -> QuerySequence:
    """
    Apply SNPs and small indels from a VCF to a reference FASTA, then return
    the mutated sequence as a QuerySequence.

    Limitations (v1): SNPs and small indels only; structural variants are not
    supported.
    """
    try:
        import vcf as pyvcf
    except ImportError:
        raise ImportError("VCF input requires PyVCF3: pip install pyvcf3")

    ref_records = {r.id: r for r in SeqIO.parse(str(ref_path), "fasta")}

    variants: dict[str, list] = {}
    with open(vcf_path) as fh:
        reader = pyvcf.Reader(fh)
        for rec in reader:
            variants.setdefault(rec.CHROM, []).append(rec)

    if not variants:
        raise ValueError("No variants found in VCF file.")

    chrom = next(iter(variants))
    if chrom not in ref_records:
        raise ValueError(
            f"Chromosome '{chrom}' in VCF not found in reference FASTA. "
            "Ensure sequence identifiers match."
        )

    ref_seq = list(str(ref_records[chrom].seq).upper())
    offset = 0

    for variant in sorted(variants[chrom], key=lambda v: v.POS):
        pos = variant.POS - 1 + offset          # VCF is 1-based
        ref_allele = str(variant.REF)
        alt_allele = str(variant.ALT[0]) if variant.ALT else ref_allele

        actual_ref = "".join(ref_seq[pos: pos + len(ref_allele)])
        if actual_ref != ref_allele:
            logger.warning(
                f"REF mismatch at {chrom}:{variant.POS} — "
                f"expected {ref_allele}, found {actual_ref}. Skipping."
            )
            continue

        ref_seq[pos: pos + len(ref_allele)] = list(alt_allele)
        offset += len(alt_allele) - len(ref_allele)

    from Bio.Seq import Seq
    mutated = "".join(ref_seq)
    return QuerySequence(
        sequence=mutated,
        identifier=f"{chrom}_mutated",
        length=len(mutated),
        gc_content=gc_fraction(Seq(mutated)),
        source_format=InputFormat.VCF,
    )


def parse_input(
    input_path: Path,
    input_format: InputFormat,
    ref_path: Optional[Path] = None,
) -> QuerySequence:
    if input_format == InputFormat.FASTA:
        return parse_fasta(input_path)
    if input_format == InputFormat.VCF:
        if ref_path is None:
            raise ValueError("--ref is required when --input-format is vcf")
        return parse_vcf(input_path, ref_path)
    raise ValueError(f"Unsupported input format: {input_format}")


# ---------------------------------------------------------------------------
# Host resolution
# ---------------------------------------------------------------------------

def resolve_host(
    identifier: str,
    user_gc: Optional[float] = None,
    entrez_email: str = "",
    entrez_api_key: str = "",
) -> HostProfile:
    """
    Resolve host GC content.

    Priority:
      1. User-supplied --host-gc (no network needed)
      2. NCBI Entrez lookup via assembly docsum
      3. NCBI Entrez FASTA fetch (slow fallback)
    """
    if user_gc is not None:
        return HostProfile(identifier=identifier, gc_content=user_gc, source="user_supplied")

    if not entrez_email:
        raise ValueError(
            "Either --host-gc or --entrez-email must be provided. "
            "Use --host-gc 0.52 to supply GC content directly."
        )

    from Bio import Entrez
    Entrez.email = entrez_email
    if entrez_api_key:
        Entrez.api_key = entrez_api_key

    logger.info(f"Fetching host GC for {identifier!r} from NCBI Entrez...")

    # Try lightweight docsum first
    gc = _gc_from_docsum(identifier)
    if gc is not None:
        return HostProfile(identifier=identifier, gc_content=gc, source="ncbi_entrez_docsum")

    # Fall back to FASTA fetch and compute
    gc = _gc_from_fasta_fetch(identifier)
    if gc is not None:
        return HostProfile(identifier=identifier, gc_content=gc, source="ncbi_entrez_fasta")

    raise RuntimeError(
        f"Could not retrieve GC content for {identifier!r} from NCBI. "
        "Use --host-gc to supply it directly."
    )


def _gc_from_docsum(identifier: str) -> Optional[float]:
    import re
    from Bio import Entrez
    try:
        handle = Entrez.esummary(db="assembly", term=identifier)
        summary = Entrez.read(handle)
        handle.close()
        meta = str(summary.get("DocumentSummarySet", {})
                   .get("DocumentSummary", [{}])[0]
                   .get("Meta", ""))
        match = re.search(r"GC[^0-9]*([0-9]+(?:\.[0-9]+)?)", meta)
        if match:
            return float(match.group(1)) / 100.0
    except Exception as exc:
        logger.debug(f"Docsum lookup failed: {exc}")
    return None


def _gc_from_fasta_fetch(identifier: str) -> Optional[float]:
    from Bio import Entrez
    from Bio.Seq import Seq
    try:
        handle = Entrez.efetch(
            db="nucleotide", id=identifier, rettype="fasta", retmode="text"
        )
        seq_str = "".join(str(r.seq) for r in SeqIO.parse(handle, "fasta"))
        handle.close()
        if seq_str:
            return gc_fraction(Seq(seq_str))
    except Exception as exc:
        logger.debug(f"FASTA fetch failed: {exc}")
    return None
