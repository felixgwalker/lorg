"""VCF writer for proxy-species edit-burden variants.

Exports variants as VCF 4.2.  Supports both:
  - write_vcf(variants, output_path, sample_name)   ← spec API
  - write_vcf(variants, output_dir, ...)             ← legacy pipeline API
"""

from __future__ import annotations

import os
from datetime import date


def write_vcf(
    variants,
    output_path: str,
    sample_name: str = "SAMPLE",
    *,
    # Legacy keyword args kept for backward compatibility with pipeline.py.
    proxy_name: str | None = None,
    target_name: str | None = None,
) -> str:
    """Write a VCF 4.2 file for *variants*.

    Parameters
    ----------
    variants : list
        List of variant dicts (legacy pipeline) or Variant dataclass instances.
    output_path : str
        Either a full file path (ending in .vcf) **or** a directory.
        When a directory is given, the file is written as
        ``<output_path>/variants.vcf`` for backward compatibility.
    sample_name : str
        Sample column name used in the VCF header.

    Returns
    -------
    str
        Absolute path to the written VCF file.
    """
    # Resolve final file path.
    if os.path.isdir(output_path) or (
        not output_path.endswith(".vcf") and not os.path.splitext(output_path)[1]
    ):
        path = os.path.join(output_path, "variants.vcf")
    else:
        path = output_path

    today = date.today().strftime("%Y%m%d")
    ref_label = proxy_name or sample_name

    # Collect contig names from variants.
    chroms_seen: list[str] = []
    seen_set: set[str] = set()
    for v in variants:
        chrom = (v.chrom if hasattr(v, "chrom") else v.get("chrom", "chrUn"))
        if chrom not in seen_set:
            seen_set.add(chrom)
            chroms_seen.append(chrom)

    with open(path, "w", newline="\n") as fh:
        # --- Meta-information lines ---
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##fileDate={today}\n")
        fh.write("##source=proxy-species-edit-burden-calculator\n")
        fh.write(f"##reference={ref_label}\n")

        for chrom in chroms_seen:
            fh.write(f"##contig=<ID={chrom}>\n")

        fh.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant class">\n')
        fh.write('##INFO=<ID=IMPACT,Number=1,Type=String,Description="Impact category">\n')
        fh.write('##INFO=<ID=WEIGHT,Number=1,Type=Integer,Description="Edit burden weight">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

        # --- Header line ---
        fh.write(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
        )

        # --- Variant records ---
        for v in variants:
            _is_dc = hasattr(v, "__dataclass_fields__")

            if _is_dc:
                chrom = v.chrom or "chrUn"
                pos = v.position + 1          # convert to 1-based
                ref = v.ref_allele or "N"
                alt = v.alt_allele or "."
                vc = v.variant_class or "SNV"
                impact = v.impact_category or "intergenic"
                weight = v.weight
            else:
                chrom = v.get("chrom", "chrUn")
                pos = v.get("pos", 0) + 1     # convert to 1-based
                ref = v.get("ref", "N") or "N"
                alt = v.get("alt", ".") or "."
                vc = v.get("variant_class", "SNV")
                impact = v.get("impact_category", "intergenic")
                weight = v.get("weight", 1)

            # Ensure REF and ALT are valid VCF tokens (not bare dots unless
            # they are the missing-value placeholder in the ALT column only).
            if ref == ".":
                ref = "N"

            info = f"SVTYPE={vc};IMPACT={impact};WEIGHT={weight}"
            gt = "0/1"

            fh.write(
                f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\tGT\t{gt}\n"
            )

    return path
