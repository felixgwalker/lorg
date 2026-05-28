import os
import sys

from .locus_reader import parse_bed, demo_loci, make_demo_locus
from .sequence_analyzer import make_demo_sequences, analyze_sequence
from .feature_annotator import parse_gff3, demo_genes, annotate_locus
from .tolerance_scorer import scan_locus, score_tolerance
from .report import (
    write_tolerance_scores, write_ranked_sites,
    write_composite_scores, write_composite_ranked,
)
from .plot import plot_locus_context, plot_composite_scores, plot_score_components


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


def _extract_seq(sequences, chrom, start, end, flank_bp=5000):
    """Return the sequence for a locus plus flanking context."""
    seq = sequences.get(chrom, "")
    chrom_len = len(seq)
    ctx_start = max(0, start - flank_bp)
    ctx_end = min(chrom_len, end + flank_bp)
    return seq[ctx_start:ctx_end]


def run_pipeline(args):
    os.makedirs(args.output_dir, exist_ok=True)

    insert_size_bp = getattr(args, "insert_size", 10000)

    if args.demo:
        loci = demo_loci()
        sequences = make_demo_sequences()
        genes = demo_genes()
        gff_path = None
    else:
        if not args.bed or not args.fasta:
            print("ERROR: --bed and --fasta are required unless --demo is used.",
                  file=sys.stderr)
            return {"status": "error", "message": "missing inputs"}
        loci = parse_bed(args.bed)
        sequences = _load_fasta(args.fasta)
        genes = []
        gff_path = getattr(args, "gff3", None)
        if gff_path:
            genes = parse_gff3(gff_path)

    # -----------------------------------------------------------------------
    # Composite scoring path (new model)
    # -----------------------------------------------------------------------
    composite_results = []
    for locus in loci:
        chrom = locus["chrom"]
        start = locus["start"]
        end = locus["end"]

        context_seq = _extract_seq(sequences, chrom, start, end)
        seq_features = analyze_sequence(context_seq)

        genomic_features = annotate_locus(chrom, start, end, gff_path)

        scores = score_tolerance(seq_features, genomic_features, insert_size_bp)

        row = {
            "locus_name": locus["name"],
            "chrom": chrom,
            "start": start,
            "end": end,
            "insert_size_bp": insert_size_bp,
            "n_genes": sum(1 for f in genomic_features if f["feature_type"] == "gene"),
            "n_exons": sum(1 for f in genomic_features if f["feature_type"] == "exon"),
            "n_regulatory": sum(1 for f in genomic_features
                                if f["feature_type"] == "regulatory"),
        }
        row.update(scores)
        composite_results.append(row)

    composite_scores_path = write_composite_scores(composite_results, args.output_dir)
    composite_ranked_path = write_composite_ranked(composite_results, args.output_dir)

    # -----------------------------------------------------------------------
    # Legacy window-scan path (original additive model)
    # -----------------------------------------------------------------------
    all_results = []
    for locus in loci:
        results = scan_locus(sequences, genes, locus)
        all_results.extend(results)

    scores_path = write_tolerance_scores(all_results, args.output_dir)
    ranked_path = write_ranked_sites(all_results, args.output_dir)

    # -----------------------------------------------------------------------
    # Plots
    # -----------------------------------------------------------------------
    plot_paths = []
    if not getattr(args, "no_plot", False):
        plot_paths += plot_composite_scores(composite_results, args.output_dir)
        plot_paths += plot_score_components(composite_results, args.output_dir)
        plot_paths += plot_locus_context(all_results, loci, args.output_dir)

    n_loci = len(loci)
    n_windows = len(all_results)
    print(f"Scored {n_loci} loci ({n_windows} scan windows).")
    print(f"Outputs written to: {args.output_dir}")

    for row in composite_results:
        tier = row["tolerance_tier"].upper()
        score = row["composite_tolerance"]
        print(f"  {row['locus_name']}: composite={score:.3f}  tier={tier}")

    all_outputs = (
        [composite_scores_path, composite_ranked_path, scores_path, ranked_path]
        + plot_paths
    )

    return {
        "status": "ok",
        "n_loci": n_loci,
        "n_windows": n_windows,
        "composite_results": composite_results,
        "outputs": all_outputs,
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Large Insert Tolerance Predictor")
    parser.add_argument("--bed", help="BED file of target loci")
    parser.add_argument("--fasta", help="Reference FASTA file")
    parser.add_argument("--gff3", help="Optional GFF3 annotation file")
    parser.add_argument("--output-dir", default="output", dest="output_dir")
    parser.add_argument("--insert-size", type=int, default=10000, dest="insert_size",
                        help="Insert size in bp (default: 10000)")
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--no-plot", action="store_true", dest="no_plot")
    args = parser.parse_args()
    result = run_pipeline(args)
    return 0 if result.get("status") == "ok" else 1
