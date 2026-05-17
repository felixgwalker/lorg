import os
import numpy as np

from .alignment_reader import load_gene_alignments, generate_synthetic_alignments
from .dnds_calculator import compute_gene_dnds
from .lrt_tester import test_all_genes
from .fdr_corrector import apply_fdr
from .report import write_dnds_table, write_selected_genes, write_summary_json
from .plot import plot_dnds


def run_pipeline(args):
    os.makedirs(args.output_dir, exist_ok=True)
    rng = np.random.default_rng(getattr(args, "seed", 42))

    if args.demo:
        gene_alignments = generate_synthetic_alignments(
            n_genes=20, n_species=5, seq_len=300, n_positive=3, rng=rng
        )
    else:
        if not args.input:
            raise ValueError("Provide --input or use --demo")
        gene_alignments = load_gene_alignments(args.input)

    dnds_results = {}
    for gene_name, sequences in gene_alignments.items():
        dnds_results[gene_name] = compute_gene_dnds(sequences)

    tested_genes = test_all_genes(dnds_results)
    tested_genes = apply_fdr(tested_genes)

    dnds_path = write_dnds_table(tested_genes, args.output_dir)
    selected_path = write_selected_genes(tested_genes, args.output_dir)
    summary_path = write_summary_json(tested_genes, args.output_dir)

    plot_paths = []
    if not getattr(args, "no_plot", False):
        fmt = getattr(args, "plot_format", "png")
        p = plot_dnds(tested_genes, args.output_dir, fmt=fmt)
        plot_paths.append(p)

    return {
        "gene_alignments": gene_alignments,
        "dnds_results": dnds_results,
        "tested_genes": tested_genes,
        "outputs": {
            "dnds_table": dnds_path,
            "selected_genes": selected_path,
            "summary_json": summary_path,
            "plots": plot_paths,
        },
    }


def main():
    import sys
    import argparse

    parser = argparse.ArgumentParser(
        description="Positive Selection Signal Detector"
    )
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--input", type=str, default=None)
    parser.add_argument("--output-dir", dest="output_dir", type=str, default="output_detector")
    parser.add_argument("--no-plot", dest="no_plot", action="store_true")
    parser.add_argument("--plot-format", dest="plot_format", choices=["png", "svg"], default="png")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    if not args.demo and args.input is None:
        parser.error("Provide --input or use --demo.")

    result = run_pipeline(args)
    n_sig = sum(1 for g in result["tested_genes"] if g.get("significant", False))
    print(f"Analyzed {len(result['tested_genes'])} genes")
    print(f"Significant positive selection signals: {n_sig}")
    print(f"Outputs written to: {args.output_dir}")
    return 0
