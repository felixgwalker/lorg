import json
import os
import pandas as pd


def write_dnds_table(tested_genes, output_dir):
    path = os.path.join(output_dir, "dnds_table.csv")
    rows = []
    for g in tested_genes:
        rows.append({
            "gene": g["gene"],
            "mean_dN": round(g["mean_dN"], 6),
            "mean_dS": round(g["mean_dS"], 6),
            "omega": round(g["omega"], 4),
            "lrt_pval": round(g["lrt_pval"], 6),
            "adjusted_pval": round(g.get("adjusted_pval", g["lrt_pval"]), 6),
            "selection_signal": g["selection_signal"],
            "significant": g.get("significant", False),
        })
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)
    return path


def write_selected_genes(tested_genes, output_dir):
    path = os.path.join(output_dir, "selected_genes.csv")
    rows = [
        {"gene": g["gene"], "omega": round(g["omega"], 4),
         "adjusted_pval": round(g.get("adjusted_pval", 1.0), 6)}
        for g in tested_genes
        if g.get("significant", False)
    ]
    df = pd.DataFrame(rows) if rows else pd.DataFrame(columns=["gene", "omega", "adjusted_pval"])
    df.to_csv(path, index=False)
    return path


def write_summary_json(tested_genes, output_dir):
    path = os.path.join(output_dir, "summary.json")
    n_genes = len(tested_genes)
    n_significant = sum(1 for g in tested_genes if g.get("significant", False))
    omegas = [g["omega"] for g in tested_genes]
    summary = {
        "n_genes_analyzed": n_genes,
        "n_significant_positive_selection": n_significant,
        "mean_omega": round(float(sum(omegas) / max(n_genes, 1)), 4),
        "max_omega": round(max(omegas) if omegas else 0.0, 4),
        "significant_genes": [g["gene"] for g in tested_genes if g.get("significant", False)],
    }
    with open(path, "w") as fh:
        json.dump(summary, fh, indent=2)
    return path
