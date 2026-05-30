"""Parse VCF GT fields to per-individual genotype arrays, or generate demo data."""

from __future__ import annotations

import random
from dataclasses import dataclass


@dataclass
class GenotypeData:
    """Per-individual genotype data across all SNPs."""
    individual_id: str
    chrom: str
    positions: list[int]
    is_homozygous: list[bool]


def parse_vcf(
    vcf_path: str,
    target_chrom: str | None = None,
    min_gq: int = 20,
    min_dp: int = 5,
) -> list[GenotypeData]:
    """Parse a VCF file with GT fields into per-individual GenotypeData.

    Parameters
    ----------
    vcf_path:
        Path to the VCF file (plain text or gzip-compressed).
    target_chrom:
        If provided, only parse records on this chromosome.
    min_gq:
        Minimum genotype quality (GQ field). Sites with GQ < min_gq for a
        given sample are treated as missing (non-homozygous).
    min_dp:
        Minimum read depth (DP field). Sites with DP < min_dp are treated as
        missing.
    """
    import gzip

    samples: list[str] = []
    # records: list of (chrom, pos, list-of-raw-sample-columns)
    records: list[tuple[str, int, list[str], list[str]]] = []

    opener = gzip.open if vcf_path.endswith(".gz") else open

    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                samples = cols[9:]
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 + len(samples):
                continue
            chrom = cols[0]
            if target_chrom is not None and chrom != target_chrom:
                continue
            pos = int(cols[1])
            fmt_fields = cols[8].split(":")
            sample_cols = cols[9:9 + len(samples)]
            records.append((chrom, pos, fmt_fields, sample_cols))

    if not samples or not records:
        return []

    chroms_seen: list[str] = []
    seen_set: set[str] = set()
    for r in records:
        if r[0] not in seen_set:
            chroms_seen.append(r[0])
            seen_set.add(r[0])

    result: list[GenotypeData] = []
    for chrom in chroms_seen:
        chrom_records = [(r[1], r[2], r[3]) for r in records if r[0] == chrom]
        for si, sample_id in enumerate(samples):
            positions: list[int] = []
            homo_flags: list[bool] = []
            for pos, fmt_fields, sample_cols in chrom_records:
                s_col = sample_cols[si]
                parts = s_col.split(":")

                # Resolve FORMAT field indices
                gt_idx = fmt_fields.index("GT") if "GT" in fmt_fields else 0
                gq_idx = fmt_fields.index("GQ") if "GQ" in fmt_fields else None
                dp_idx = fmt_fields.index("DP") if "DP" in fmt_fields else None

                # Filter by GQ
                if gq_idx is not None and gq_idx < len(parts):
                    try:
                        if int(parts[gq_idx]) < min_gq:
                            continue
                    except ValueError:
                        continue  # missing GQ value

                # Filter by DP
                if dp_idx is not None and dp_idx < len(parts):
                    try:
                        if int(parts[dp_idx]) < min_dp:
                            continue
                    except ValueError:
                        continue  # missing DP value

                gt = parts[gt_idx] if gt_idx < len(parts) else "."
                positions.append(pos)
                homo_flags.append(_is_homozygous(gt))

            result.append(GenotypeData(
                individual_id=sample_id,
                chrom=chrom,
                positions=positions,
                is_homozygous=homo_flags,
            ))
    return result


def make_demo_genotype_data(
    seed: int = 42,
) -> list[GenotypeData]:
    """Generate synthetic GenotypeData with realistic ROH patterns.

    Produces three individuals:
      - IND_01: highly inbred (large long ROH, multiple medium ROH per chrom)
      - IND_02: moderately inbred (several medium ROH)
      - IND_03: outbred baseline (only short ROH / background)

    Returns a flat list of GenotypeData, one entry per (individual, chrom).
    """
    return generate_demo_genotypes(n_individuals=3, n_snps_per_chrom=300,
                                   n_chroms=5, seed=seed)


def generate_demo_genotypes(
    n_individuals: int = 3,
    n_snps_per_chrom: int = 300,
    n_chroms: int = 5,
    seed: int = 42,
) -> list[GenotypeData]:
    """Generate synthetic genotype data with embedded ROH regions.

    Individual 0 is highly inbred: large long ROH on every chromosome.
    Individual 1 is moderately inbred: medium ROH on most chromosomes.
    Individual 2 is outbred: only low-level background homozygosity.
    """
    rng = random.Random(seed)
    result: list[GenotypeData] = []
    # Use human-like chromosome lengths (autosomes 1–5 scaled down)
    chrom_lengths = [248_956_422, 242_193_529, 198_295_559, 190_214_555, 181_538_259]
    chrom_lengths = chrom_lengths[:n_chroms]

    for ind_idx in range(n_individuals):
        ind_id = f"IND_{ind_idx + 1:02d}"
        for chrom_idx in range(n_chroms):
            chrom = f"chr{chrom_idx + 1}"
            chrom_len = chrom_lengths[chrom_idx]
            positions = sorted(rng.sample(range(500_000, chrom_len - 500_000), n_snps_per_chrom))
            # Base background heterozygosity: ~15 % homozygous by chance
            is_homo = [rng.random() < 0.15 for _ in positions]

            if ind_idx == 0:
                # Highly inbred individual: inject 2–3 long ROH (>1 Mb) per chrom
                n_roh = rng.randint(2, 3)
                for _ in range(n_roh):
                    frac_start = rng.uniform(0.05, 0.70)
                    frac_len = rng.uniform(0.06, 0.20)   # 6–20 % of chrom → multi-Mb
                    roh_start = int(chrom_len * frac_start)
                    roh_end = int(chrom_len * (frac_start + frac_len))
                    for i, p in enumerate(positions):
                        if roh_start <= p <= roh_end:
                            is_homo[i] = rng.random() < 0.98

            elif ind_idx == 1:
                # Moderately inbred: 1–2 medium ROH (100 kb–1 Mb) per chrom
                n_roh = rng.randint(1, 2)
                for _ in range(n_roh):
                    frac_start = rng.uniform(0.05, 0.80)
                    # 0.5–3 % of chrom ≈ a few hundred kb to ~1 Mb
                    frac_len = rng.uniform(0.005, 0.030)
                    roh_start = int(chrom_len * frac_start)
                    roh_end = int(chrom_len * (frac_start + frac_len))
                    for i, p in enumerate(positions):
                        if roh_start <= p <= roh_end:
                            is_homo[i] = rng.random() < 0.97

            # ind_idx == 2: outbred — leave at background level

            result.append(GenotypeData(
                individual_id=ind_id,
                chrom=chrom,
                positions=positions,
                is_homozygous=is_homo,
            ))
    return result


def _is_homozygous(gt: str) -> bool:
    for sep in ("/", "|"):
        if sep in gt:
            alleles = gt.split(sep)
            non_missing = [a for a in alleles if a != "."]
            if len(non_missing) >= 2:
                return len(set(non_missing)) == 1
            return False
    return False
