import json
import os
import pandas as pd


def write_burden_summary(burden, output_dir):
    path = os.path.join(output_dir, "edit_burden_summary.csv")
    rows = []
    weights = {"SNV": 1, "SMALL_INS": 3, "SMALL_DEL": 3,
                "LARGE_INS": 10, "LARGE_DEL": 10, "SV_INS": 50, "SV_DEL": 50}

    for vc, count in burden["class_counts"].items():
        w = weights.get(vc, 1)
        rows.append({
            "variant_class": vc,
            "count": count,
            "weight": w,
            "weighted_contribution": count * w,
        })

    rows.append({
        "variant_class": "TOTAL",
        "count": burden["total_edits"],
        "weight": "",
        "weighted_contribution": burden["weighted_burden"],
    })

    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)
    return path


def write_prioritized_edits(variants, output_dir, top_n=None):
    path = os.path.join(output_dir, "prioritized_edits.csv")
    sorted_vars = sorted(variants, key=lambda v: -v.get("impact_score", 0))
    if top_n:
        sorted_vars = sorted_vars[:top_n]

    rows = []
    for v in sorted_vars:
        rows.append({
            "chrom": v.get("chrom", ""),
            "pos": v.get("pos", 0),
            "ref": v.get("ref", ""),
            "alt": v.get("alt", ""),
            "type": v.get("type", ""),
            "variant_class": v.get("variant_class", ""),
            "impact_category": v.get("impact_category", ""),
            "impact_score": v.get("impact_score", 0),
            "weight": v.get("weight", 1),
        })

    df = pd.DataFrame(rows) if rows else pd.DataFrame(columns=[
        "chrom", "pos", "ref", "alt", "type", "variant_class",
        "impact_category", "impact_score", "weight"
    ])
    df.to_csv(path, index=False)
    return path


def write_burden_json(burden, output_dir):
    path = os.path.join(output_dir, "edit_burden.json")
    with open(path, "w") as fh:
        json.dump(burden, fh, indent=2)
    return path
