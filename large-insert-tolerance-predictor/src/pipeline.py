import os
import sys

from .locus_reader import parse_bed, demo_loci
from .sequence_analyzer import make_demo_sequences
from .feature_annotator import parse_gff3, demo_genes
from .tolerance_scorer import scan_locus
from .report import write_tolerance_scores, write_ranked_sites
from .plot import plot_locus_context


def _load_fasta(fasta_path):
    sequences = {}
    current_chrom = None
    parts = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if current_chrom is not None:
                    sequences[current_chrom] = "".join(parts)
                current_chrom = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())
    if current_chrom is not None:
        sequences[current_chrom] = "".join(parts)
    return sequences


def run_pipeline(args):
    os.makedirs(args.output_dir, exist_ok=True)

    if args.demo:
        loci = demo_loci()
        sequences = make_demo_sequences()
        genes = demo_genes()
    else:
        if not args.bed or not args.fasta:
            print("ERROR: --bed and --fasta are required unless --demo is used.", file=sys.stderr)
            return {"status": "error", "message": "missing inputs"}
        loci = parse_bed(args.bed)
        sequences = _load_fasta(args.fasta)
        genes = []
        if args.gff3:
            genes = parse_gff3(args.gff3)

    all_results = []
    for locus in loci:
        results = scan_locus(sequences, genes, locus)
        all_results.extend(results)

    scores_path = write_tolerance_scores(all_results, args.output_dir)
    ranked_path = write_ranked_sites(all_results, args.output_dir)

    plot_paths = []
    if not getattr(args, "no_plot", False):
        plot_paths = plot_locus_context(all_results, loci, args.output_dir)

    print(f"Scored {len(all_results)} windows across {len(loci)} loci.")
    print(f"Outputs written to: {args.output_dir}")

    return {
        "status": "ok",
        "n_loci": len(loci),
        "n_windows": len(all_results),
        "outputs": [scores_path, ranked_path] + plot_paths,
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Large Insert Tolerance Predictor")
    parser.add_argument("--bed", help="BED file of target loci")
    parser.add_argument("--fasta", help="Reference FASTA file")
    parser.add_argument("--gff3", help="Optional GFF3 annotation file")
    parser.add_argument("--output-dir", default="output", dest="output_dir")
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--no-plot", action="store_true", dest="no_plot")
    args = parser.parse_args()
    result = run_pipeline(args)
    return 0 if result.get("status") == "ok" else 1
