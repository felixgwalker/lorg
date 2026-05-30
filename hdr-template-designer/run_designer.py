"""CLI entry point for hdr-template-designer."""

import argparse
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run_designer",
        description="Design HDR donor templates for CRISPR edits.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--reference", type=Path, metavar="FASTA",
                   help="Reference sequence FASTA file")
    p.add_argument("--chrom", metavar="CHROM",
                   help="Chromosome/contig name to use from reference")
    p.add_argument("--cut-pos", type=int, metavar="INT",
                   help="Cut site position (0-based, in reference coordinates)")
    p.add_argument("--edit-type", choices=["snv", "ins", "del"], default="snv",
                   help="Type of edit to introduce")
    p.add_argument("--ref-allele", metavar="STR",
                   help="Reference allele at cut position")
    p.add_argument("--alt-allele", metavar="STR",
                   help="Alternate (desired) allele")
    p.add_argument("--arm-len", type=int, default=100,
                   help="Homology arm length in bp (each side)")
    p.add_argument("--silent-pam", action="store_true",
                   help="Introduce synonymous PAM-disrupting mutation in right arm")
    p.add_argument("--output-dir", type=Path, default=Path("results"),
                   metavar="PATH", help="Output directory")
    p.add_argument("--no-plot", action="store_true", help="Skip template diagram")
    p.add_argument("--plot-format", choices=["png", "svg"], default="png",
                   help="Output plot format")
    p.add_argument("--demo", action="store_true",
                   help="Run demo: 500bp synthetic reference, C->T SNV at pos 250, 100bp arms")
    return p


def validate_args(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    if args.demo:
        return
    if args.reference is None:
        parser.error("--reference is required unless --demo is set")
    if not args.reference.exists():
        parser.error(f"Reference file not found: {args.reference}")
    if args.chrom is None:
        parser.error("--chrom is required unless --demo is set")
    if args.cut_pos is None:
        parser.error("--cut-pos is required unless --demo is set")
    if args.ref_allele is None:
        parser.error("--ref-allele is required unless --demo is set")
    if args.alt_allele is None:
        parser.error("--alt-allele is required unless --demo is set")
    if args.arm_len < 10:
        parser.error("--arm-len must be at least 10 bp")


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    validate_args(args, parser)

    from src.reference_reader import load_fasta, make_demo_locus
    from src.pipeline import run_pipeline

    if args.demo:
        ref_seq, cut_pos = make_demo_locus()
        edit_type, ref_allele, alt_allele, silent_pam = "snv", "C", "T", True
        print(
            f"Demo mode: 1000 bp synthetic locus | C->T SNP at cut pos {cut_pos} "
            f"| PAM disruption ON"
        )
    else:
        genome = load_fasta(args.reference)
        if args.chrom not in genome:
            parser.error(f"Chromosome '{args.chrom}' not found; available: {list(genome.keys())}")
        ref_seq = genome[args.chrom]
        cut_pos = args.cut_pos
        edit_type = args.edit_type
        ref_allele = args.ref_allele
        alt_allele = args.alt_allele
        silent_pam = args.silent_pam
        print(f"Reference: {args.chrom} ({len(ref_seq):,} bp) | cut={cut_pos} | {ref_allele}->{alt_allele}")

    result = run_pipeline(
        ref_seq=ref_seq,
        cut_pos=cut_pos,
        edit_type=edit_type,
        ref_allele=ref_allele,
        alt_allele=alt_allele,
        arm_len=args.arm_len,
        silent_pam=silent_pam,
        output_dir=args.output_dir,
        no_plot=args.no_plot,
        plot_fmt=args.plot_format,
    )

    print(f"Template: {result['template_len']} bp  "
          f"(left={result['left_arm_len']} | edit={result['edit_seq']} | right={result['right_arm_len']})")
    print(f"PAM:      {result['pam_note']}")
    print(f"QC flags: {result['qc_flags']}  |  QC passed: {result['qc_passed']}")
    if result.get("qc_warnings"):
        for w in result["qc_warnings"]:
            print(f"  [WARN] {w}")
    print(f"PAM mutation sites found: {result['pam_mutations_found']}")
    print(f"  FASTA:              {result['fasta']}")
    print(f"  Annotated sequence: {result['annotated_sequence']}")
    print(f"  QC report:          {result['qc_report']}")
    print(f"  Arm variants:       {result['arm_variants']}")
    if result["plot"]:
        print(f"  Template diagram:   {result['plot']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
