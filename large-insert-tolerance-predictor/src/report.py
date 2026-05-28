import os
import csv


# ---------------------------------------------------------------------------
# Composite scoring report (new score_tolerance() model)
# ---------------------------------------------------------------------------

_COMPOSITE_FIELDS = [
    "locus_name", "chrom", "start", "end", "insert_size_bp",
    "chromatin_score", "sequence_complexity_score", "gene_density_score",
    "size_penalty", "composite_tolerance", "tolerance_tier",
    "n_genes", "n_exons", "n_regulatory",
]


def _composite_rationale(row):
    parts = []
    if float(row.get("chromatin_score", 0)) >= 0.7:
        parts.append("favourable chromatin context")
    if float(row.get("sequence_complexity_score", 0)) >= 0.7:
        parts.append("low structural motif burden")
    if float(row.get("gene_density_score", 0)) >= 0.7:
        parts.append("no overlapping genes")
    if float(row.get("size_penalty", 0)) >= 0.7:
        parts.append("insert within safe size range")
    tier = row.get("tolerance_tier", "")
    if not parts:
        parts.append(f"{tier} tolerance (no strongly favourable features)")
    return "; ".join(parts)


def write_composite_scores(locus_results, output_dir):
    """Write per-locus composite tolerance scores to CSV.

    locus_results is a list of dicts as built by run_pipeline (composite path).
    """
    path = os.path.join(output_dir, "composite_tolerance_scores.csv")
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_COMPOSITE_FIELDS,
                                extrasaction="ignore")
        writer.writeheader()
        for row in locus_results:
            writer.writerow({k: row.get(k, "") for k in _COMPOSITE_FIELDS})
    return path


def write_composite_ranked(locus_results, output_dir, top_n=10):
    """Write ranked composite sites with rationale."""
    sorted_results = sorted(
        locus_results,
        key=lambda r: float(r.get("composite_tolerance", 0)),
        reverse=True,
    )
    top = sorted_results[:top_n]
    path = os.path.join(output_dir, "composite_ranked_sites.csv")
    ranked_fields = ["rank", "rationale"] + _COMPOSITE_FIELDS
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=ranked_fields, extrasaction="ignore")
        writer.writeheader()
        for rank, row in enumerate(top, 1):
            out = {k: row.get(k, "") for k in _COMPOSITE_FIELDS}
            out["rank"] = rank
            out["rationale"] = _composite_rationale(row)
            writer.writerow(out)
    return path


# ---------------------------------------------------------------------------
# Legacy window-scan report (original additive scoring model)
# ---------------------------------------------------------------------------

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
