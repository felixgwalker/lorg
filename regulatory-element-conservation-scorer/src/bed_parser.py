def parse_bed(bed_path):
    elements = []
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
            elements.append({"id": name, "chrom": chrom, "start": start, "end": end})
    return elements


def demo_elements():
    elements = []
    specs = [
        ("chr1", 1000, 1200, "enhancer_01"),
        ("chr1", 5000, 5150, "promoter_01"),
        ("chr1", 12000, 12300, "enhancer_02"),
        ("chr1", 25000, 25100, "silencer_01"),
        ("chr2", 3000, 3250, "enhancer_03"),
        ("chr2", 8000, 8200, "promoter_02"),
        ("chr2", 15000, 15180, "insulator_01"),
        ("chr2", 30000, 30120, "enhancer_04"),
        ("chr3", 2000, 2300, "promoter_03"),
        ("chr3", 9000, 9200, "enhancer_05"),
        ("chr3", 18000, 18150, "silencer_02"),
        ("chr3", 27000, 27250, "enhancer_06"),
    ]
    for chrom, start, end, name in specs:
        elements.append({"id": name, "chrom": chrom, "start": start, "end": end})
    return elements
