"""Parse VCF GT fields to per-individual genotype arrays, or generate demo data."""

from __future__ import annotations

import random
from dataclasses import dataclass, field


@dataclass
class GenotypeData:
    """Per-individual genotype data across all SNPs."""
    individual_id: str
    chrom: str
    positions: list[int]
    is_homozygous: list[bool]


def parse_vcf(vcf_path: str) -> list[GenotypeData]:
    """Parse a VCF file with GT fields into per-individual GenotypeData."""
    samples: list[str] = []
    records: list[tuple[str, int, list[str]]] = []

    with open(vcf_path) as fh:
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
            pos = int(cols[1])
            fmt = cols[8].split(":")
            gt_idx = fmt.index("GT") if "GT" in fmt else 0
            gts = []
            for s_col in cols[9:9 + len(samples)]:
                gt = s_col.split(":")[gt_idx]
                gts.append(gt)
            records.append((chrom, pos, gts))

    if not samples or not records:
        return []

    chroms_seen = sorted({r[0] for r in records})
    result: list[GenotypeData] = []
    for chrom in chroms_seen:
        chrom_records = [(r[1], r[2]) for r in records if r[0] == chrom]
        for si, sample_id in enumerate(samples):
            positions = [p for p, _ in chrom_records]
            homo = [_is_homozygous(gts[si]) for _, gts in chrom_records]
            result.append(GenotypeData(
                individual_id=sample_id,
                chrom=chrom,
                positions=positions,
                is_homozygous=homo,
            ))
    return result


def generate_demo_genotypes(
    n_individuals: int = 3,
    n_snps_per_chrom: int = 200,
    n_chroms: int = 5,
    seed: int = 42,
) -> list[GenotypeData]:
    """Generate synthetic genotype data with embedded ROH regions."""
    rng = random.Random(seed)
    result: list[GenotypeData] = []
    chrom_lengths = [50_000_000, 45_000_000, 40_000_000, 35_000_000, 30_000_000]

    for ind_idx in range(n_individuals):
        ind_id = f"IND_{ind_idx + 1:02d}"
        for chrom_idx in range(n_chroms):
            chrom = f"chr{chrom_idx + 1}"
            chrom_len = chrom_lengths[chrom_idx]
            positions = sorted(rng.sample(range(100_000, chrom_len), n_snps_per_chrom))
            is_homo = [rng.random() < 0.15 for _ in positions]
            roh_start_frac = rng.uniform(0.1, 0.6)
            roh_len_frac = rng.uniform(0.05, 0.20)
            roh_start_pos = int(chrom_len * roh_start_frac)
            roh_end_pos = int(chrom_len * (roh_start_frac + roh_len_frac))
            for i, p in enumerate(positions):
                if roh_start_pos <= p <= roh_end_pos:
                    is_homo[i] = rng.random() < 0.97
            if ind_idx == 0 and chrom_idx == 0:
                long_start = int(chrom_len * 0.05)
                long_end = int(chrom_len * 0.30)
                for i, p in enumerate(positions):
                    if long_start <= p <= long_end:
                        is_homo[i] = rng.random() < 0.98
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
