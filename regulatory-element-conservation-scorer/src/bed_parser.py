def parse_bed(bed_path):
    """Read BED3/BED6 format.

    Returns a list of RegionRecord dicts with keys:
        chrom, start, end, name, score, strand, id
    BED3 files will have name=chrom:start-end, score=0, strand='.'.
    """
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
            score = float(parts[4]) if len(parts) > 4 else 0.0
            strand = parts[5] if len(parts) > 5 else "."
            elements.append({
                "id": name,
                "chrom": chrom,
                "start": start,
                "end": end,
                "name": name,
                "score": score,
                "strand": strand,
            })
    return elements


def make_demo_bed():
    """Return 20 synthetic regulatory elements across 3 chromosomes.

    Each record is a RegionRecord dict with keys:
        chrom, start, end, name, score, strand, id
    """
    specs = [
        # chr1 — 7 elements
        ("chr1", 1000,  1200,  "enhancer_01",   850, "+"),
        ("chr1", 5000,  5150,  "promoter_01",   920, "+"),
        ("chr1", 12000, 12300, "enhancer_02",   780, "-"),
        ("chr1", 25000, 25100, "silencer_01",   710, "-"),
        ("chr1", 38000, 38250, "enhancer_07",   670, "+"),
        ("chr1", 51000, 51180, "promoter_04",   900, "+"),
        ("chr1", 64000, 64300, "insulator_02",  760, "."),
        # chr2 — 7 elements
        ("chr2", 3000,  3250,  "enhancer_03",   820, "+"),
        ("chr2", 8000,  8200,  "promoter_02",   940, "+"),
        ("chr2", 15000, 15180, "insulator_01",  800, "."),
        ("chr2", 30000, 30120, "enhancer_04",   620, "-"),
        ("chr2", 44000, 44220, "enhancer_08",   690, "+"),
        ("chr2", 58000, 58300, "silencer_03",   730, "-"),
        ("chr2", 72000, 72150, "promoter_05",   910, "+"),
        # chr3 — 6 elements
        ("chr3", 2000,  2300,  "promoter_03",   960, "+"),
        ("chr3", 9000,  9200,  "enhancer_05",   840, "+"),
        ("chr3", 18000, 18150, "silencer_02",   700, "-"),
        ("chr3", 27000, 27250, "enhancer_06",   790, "+"),
        ("chr3", 36000, 36200, "insulator_03",  750, "."),
        ("chr3", 45000, 45180, "enhancer_09",   660, "-"),
    ]
    elements = []
    for chrom, start, end, name, score, strand in specs:
        elements.append({
            "id": name,
            "chrom": chrom,
            "start": start,
            "end": end,
            "name": name,
            "score": float(score),
            "strand": strand,
        })
    return elements


# Keep legacy alias so pipeline.py continues to work unchanged.
def demo_elements():
    """Legacy alias for make_demo_bed() — returns 20 synthetic elements."""
    return make_demo_bed()
