import numpy as np


def parse_vcf(vcf_path):
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
