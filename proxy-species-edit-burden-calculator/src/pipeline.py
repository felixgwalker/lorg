import os
import numpy as np

from .fasta_parser import parse_fasta, generate_demo_genomes
from .sequence_aligner import align_genomes
from .variant_classifier import classify_all_variants, compute_burden
from .impact_annotator import annotate_impact
from .vcf_writer import write_vcf
from .report import write_burden_summary, write_prioritized_edits, write_burden_json
from .plot import plot_edit_burden


def run_pipeline(args):
    os.makedirs(args.output_dir, exist_ok=True)
    rng = np.random.default_rng(getattr(args, "seed", 42))

    if args.demo:
        proxy_seqs, target_seqs = generate_demo_genomes(
            n_chroms=3, chrom_len=10000, n_snvs=50, n_indels=10, rng=rng
        )
        n_samples = 200
    else:
        if not args.proxy or not args.target:
            raise ValueError("Provide --proxy and --target FASTA files, or use --demo")
        proxy_seqs = parse_fasta(args.proxy)
        target_seqs = parse_fasta(args.target)
        n_samples = getattr(args, "n_samples", 1000)

    variants, chrom_lengths = align_genomes(
        proxy_seqs, target_seqs, window_size=500, n_samples=n_samples, rng=rng
    )

    variants = classify_all_variants(variants)

    gff_path = getattr(args, "gff", None)
    variants = annotate_impact(variants, gff_path=gff_path)

    total_genome_bp = sum(chrom_lengths.values())
    burden = compute_burden(variants, total_genome_bp)

    vcf_path = write_vcf(variants, args.output_dir)
    burden_csv = write_burden_summary(burden, args.output_dir)
    prioritized_path = write_prioritized_edits(variants, args.output_dir)
    burden_json = write_burden_json(burden, args.output_dir)

    plot_paths = []
    if not getattr(args, "no_plot", False):
        fmt = getattr(args, "plot_format", "png")
        p = plot_edit_burden(variants, burden, chrom_lengths, args.output_dir, fmt=fmt)
        plot_paths.append(p)

    return {
        "proxy_seqs": proxy_seqs,
        "target_seqs": target_seqs,
        "variants": variants,
        "chrom_lengths": chrom_lengths,
        "burden": burden,
        "outputs": {
            "vcf": vcf_path,
            "burden_csv": burden_csv,
            "prioritized_csv": prioritized_path,
            "burden_json": burden_json,
            "plots": plot_paths,
        },
    }


def main():
    import sys
    from run_calculator import build_parser, validate_args
    parser = build_parser()
    args = parser.parse_args()
    validate_args(args, parser)
    result = run_pipeline(args)
    b = result["burden"]
    print(f"Total variants called: {b['total_edits']}")
    print(f"Weighted edit burden: {b['weighted_burden']}")
    print(f"Normalized burden: {b['normalized_burden_per_mb']:.2f} per Mb")
    print(f"Outputs written to: {args.output_dir}")
    return 0
