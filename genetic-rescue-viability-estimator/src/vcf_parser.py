from dataclasses import dataclass, field
from typing import List
import numpy as np


@dataclass
class GenotypeData:
    individual_ids: List[str]
    positions: List[int]
    is_homozygous: np.ndarray   # bool array shape (n_individuals, n_sites)
    allele_freqs: List[float]
    chrom: str


def parse_vcf(vcf_path: str) -> GenotypeData:
    """Read a VCF file and return a GenotypeData dataclass.

    Extracts per-individual diploid genotypes for all bi-allelic SNP sites that
    have a GT field.  Missing genotypes ('./.') are treated as heterozygous so
    they do not artificially inflate ROH lengths.
    """
    sample_names: List[str] = []
    positions: List[int] = []
    chrom = ""
    raw_gts: List[List[int]] = []   # list of per-site vectors (0=het, 1=hom)
    alt_counts: List[int] = []      # alt allele counts for AF
    header_found = False

    with open(vcf_path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.split("\t")
                sample_names = fields[9:]
                header_found = True
                continue
            if not header_found:
                continue
            fields = line.split("\t")
            n_samples = len(sample_names)
            if len(fields) < 9 + n_samples:
                continue

            chrom_field = fields[0]
            pos = int(fields[1])
            fmt = fields[8].split(":")
            gt_idx = fmt.index("GT") if "GT" in fmt else 0

            if chrom == "":
                chrom = chrom_field

            site_hom = []
            alt_count = 0
            total_alleles = 0
            for i, _name in enumerate(sample_names):
                sample_field = fields[9 + i].split(":")
                gt_raw = sample_field[gt_idx] if gt_idx < len(sample_field) else "./."
                sep = "|" if "|" in gt_raw else "/"
                alleles = gt_raw.split(sep)
                if len(alleles) == 2 and alleles[0] != "." and alleles[1] != ".":
                    a0, a1 = alleles[0], alleles[1]
                    is_hom = int(a0 == a1)
                    site_hom.append(is_hom)
                    alt_count += int(a0 != "0") + int(a1 != "0")
                    total_alleles += 2
                else:
                    site_hom.append(0)  # treat missing as het

            positions.append(pos)
            raw_gts.append(site_hom)
            af = alt_count / total_alleles if total_alleles > 0 else 0.0
            alt_counts.append(af)

    n_sites = len(positions)
    n_indiv = len(sample_names)

    if n_sites == 0 or n_indiv == 0:
        return GenotypeData(
            individual_ids=sample_names,
            positions=positions,
            is_homozygous=np.zeros((n_indiv, n_sites), dtype=bool),
            allele_freqs=alt_counts,
            chrom=chrom,
        )

    # raw_gts is (n_sites, n_individuals) → transpose to (n_individuals, n_sites)
    is_hom_array = np.array(raw_gts, dtype=bool).T

    return GenotypeData(
        individual_ids=sample_names,
        positions=positions,
        is_homozygous=is_hom_array,
        allele_freqs=alt_counts,
        chrom=chrom,
    )


def make_demo_vcf_data(
    n_individuals: int = 20,
    n_snps: int = 1000,
    seed: int = 42,
) -> GenotypeData:
    """Generate synthetic GenotypeData with realistic ROH patterns.

    Five 'inbred' individuals have elevated runs of homozygosity (FROH ~0.35),
    the remaining fifteen have lower background inbreeding (FROH ~0.10).
    """
    rng = np.random.default_rng(seed)

    individual_ids = [f"IND_{i+1:03d}" for i in range(n_individuals)]
    positions = sorted(rng.integers(1, 2_700_000, size=n_snps).tolist())
    chrom = "chr1"

    is_homozygous = np.zeros((n_individuals, n_snps), dtype=bool)

    for i in range(n_individuals):
        # Baseline per-site homozygosity probability
        if i < 5:
            # High-inbreeding individuals: intersperse long ROH segments
            base_p_hom = 0.35
            n_roh_segments = rng.integers(5, 12)
        else:
            base_p_hom = 0.10
            n_roh_segments = rng.integers(0, 4)

        # Random background homozygosity
        hom_row = rng.random(n_snps) < base_p_hom

        # Overlay ROH segments (contiguous runs of homozygosity)
        for _ in range(n_roh_segments):
            seg_start = rng.integers(0, n_snps - 50)
            seg_len = rng.integers(30, min(150, n_snps - seg_start))
            hom_row[seg_start:seg_start + seg_len] = True

        is_homozygous[i] = hom_row

    # Allele frequencies: random MAF in [0.05, 0.5]
    allele_freqs = rng.uniform(0.05, 0.5, size=n_snps).tolist()

    return GenotypeData(
        individual_ids=individual_ids,
        positions=positions,
        is_homozygous=is_homozygous,
        allele_freqs=allele_freqs,
        chrom=chrom,
    )


def parse_vcf_legacy(vcf_path):
    """Legacy list-of-dicts interface used by the pipeline."""
    individuals = []
    genotypes_by_indiv = {}
    header_found = False
    sample_names = []

    with open(vcf_path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.split("\t")
                sample_names = fields[9:]
                for name in sample_names:
                    genotypes_by_indiv[name] = []
                header_found = True
                continue
            if not header_found:
                continue
            fields = line.split("\t")
            if len(fields) < 9 + len(sample_names):
                continue
            fmt = fields[8].split(":")
            gt_idx = fmt.index("GT") if "GT" in fmt else 0
            for i, name in enumerate(sample_names):
                sample_field = fields[9 + i].split(":")
                gt_raw = sample_field[gt_idx] if gt_idx < len(sample_field) else "./."
                sep = "|" if "|" in gt_raw else "/"
                alleles = gt_raw.split(sep)
                genotypes_by_indiv[name].append(alleles)

    results = []
    for name in sample_names:
        gts = genotypes_by_indiv[name]
        n_sites = len(gts)
        if n_sites == 0:
            results.append({"individual": name, "F_initial": 0.5, "H_initial": 0.5})
            continue
        het_count = sum(
            1 for a in gts
            if len(a) == 2 and a[0] != "." and a[1] != "." and a[0] != a[1]
        )
        H = het_count / n_sites
        F = 1 - H
        results.append({"individual": name, "F_initial": F, "H_initial": H})

    return results


def generate_synthetic_population(n_individuals=20, n_snps=500, rng=None):
    if rng is None:
        rng = np.random.default_rng(42)

    base_F = 0.35
    results = []
    for i in range(n_individuals):
        indiv_F = np.clip(base_F + rng.normal(0, 0.05), 0.05, 0.95)
        H = 1 - indiv_F
        results.append({
            "individual": f"IND_{i+1:03d}",
            "F_initial": round(float(indiv_F), 4),
            "H_initial": round(float(H), 4),
        })
    return results
