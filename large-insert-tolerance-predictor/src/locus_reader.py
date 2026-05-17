import os


def parse_bed(bed_path):
    loci = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                parts = line.split()
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
            loci.append({"chrom": chrom, "start": start, "end": end, "name": name})
    return loci


def demo_loci():
    return [
        {"chrom": "chr1", "start": 50000, "end": 51000, "name": "locus_A"},
        {"chrom": "chr1", "start": 100000, "end": 101000, "name": "locus_B"},
        {"chrom": "chr2", "start": 75000, "end": 76000, "name": "locus_C"},
    ]
