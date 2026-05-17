WEIGHTS = {
    "SNV": 1,
    "SMALL_INS": 3,
    "SMALL_DEL": 3,
    "LARGE_INS": 10,
    "LARGE_DEL": 10,
    "SV_INS": 50,
    "SV_DEL": 50,
}


def classify_variant(variant):
    vtype = variant["type"]
    length = variant.get("length", 1)

    if vtype == "SNV":
        return "SNV"
    elif vtype == "INS":
        if length <= 50:
            return "SMALL_INS"
        elif length <= 500:
            return "LARGE_INS"
        else:
            return "SV_INS"
    elif vtype == "DEL":
        if length <= 50:
            return "SMALL_DEL"
        elif length <= 500:
            return "LARGE_DEL"
        else:
            return "SV_DEL"
    return "SNV"


def classify_all_variants(variants):
    for v in variants:
        v["variant_class"] = classify_variant(v)
        v["weight"] = WEIGHTS.get(v["variant_class"], 1)
    return variants


def compute_burden(variants, total_genome_bp):
    class_counts = {}
    for v in variants:
        vc = v.get("variant_class", "SNV")
        class_counts[vc] = class_counts.get(vc, 0) + 1

    total_edits = sum(class_counts.values())
    weighted_burden = sum(
        count * WEIGHTS.get(vc, 1) for vc, count in class_counts.items()
    )

    total_mb = max(total_genome_bp / 1_000_000, 1e-6)
    normalized_burden = weighted_burden / total_mb

    return {
        "class_counts": class_counts,
        "total_edits": total_edits,
        "weighted_burden": weighted_burden,
        "normalized_burden_per_mb": round(normalized_burden, 2),
        "total_genome_bp": total_genome_bp,
    }
