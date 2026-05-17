import os
import csv


def write_tolerance_scores(all_results, output_dir):
    path = os.path.join(output_dir, "tolerance_scores.csv")
    fieldnames = [
        "locus_name", "chrom", "start", "end",
        "gene_density_score", "regulatory_score",
        "repeat_score", "complexity_score", "total_score",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in all_results:
            writer.writerow({k: row[k] for k in fieldnames})
    return path


def write_ranked_sites(all_results, output_dir, top_n=10):
    sorted_results = sorted(all_results, key=lambda r: r["total_score"], reverse=True)
    top = sorted_results[:top_n]
    path = os.path.join(output_dir, "ranked_sites.csv")
    fieldnames = [
        "rank", "locus_name", "chrom", "start", "end",
        "gene_density_score", "regulatory_score",
        "repeat_score", "complexity_score", "total_score", "rationale",
    ]

    def make_rationale(row):
        parts = []
        if row["gene_density_score"] >= 20:
            parts.append("low gene density")
        if row["regulatory_score"] >= 20:
            parts.append("far from regulatory elements")
        if row["repeat_score"] >= 15:
            parts.append("low repeat content")
        if row["complexity_score"] >= 15:
            parts.append("high sequence complexity")
        return "; ".join(parts) if parts else "moderate tolerance"

    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for rank, row in enumerate(top, 1):
            out = {k: row[k] for k in fieldnames if k in row}
            out["rank"] = rank
            out["rationale"] = make_rationale(row)
            writer.writerow(out)
    return path
