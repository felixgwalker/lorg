"""Pipeline orchestrator for hdr-template-designer."""

import sys
from pathlib import Path


def run_pipeline(
    ref_seq: str,
    cut_pos: int,
    edit_type: str,
    ref_allele: str,
    alt_allele: str,
    arm_len: int = 100,
    silent_pam: bool = False,
    output_dir: Path = Path("results"),
    no_plot: bool = False,
    plot_fmt: str = "png",
) -> dict:
    from .edit_applicator import build_template
    from .pam_disruptor import disrupt_pam
    from .qc_checker import run_all_checks, arm_length_variants
    from .report import write_template_fasta, write_qc_report, write_arm_variants
    from .plot import plot_template_diagram

    output_dir.mkdir(parents=True, exist_ok=True)

    template = build_template(ref_seq, cut_pos, arm_len, arm_len, edit_type, ref_allele, alt_allele)

    pam_note = "PAM disruption not requested"
    if silent_pam:
        new_right, disrupted, pam_note = disrupt_pam(template["right_arm"])
        if disrupted:
            template["right_arm"] = new_right
            template["full_template"] = template["left_arm"] + template["edit_seq"] + new_right

    qc_checks = run_all_checks(template["left_arm"], template["right_arm"])
    variants = arm_length_variants(ref_seq, cut_pos)

    fasta_path = write_template_fasta(template, output_dir)
    qc_path = write_qc_report(qc_checks, output_dir)
    variants_path = write_arm_variants(variants, output_dir)

    plot_path = None
    if not no_plot:
        plot_path = plot_template_diagram(template, qc_checks, output_dir, fmt=plot_fmt)

    n_flags = sum(1 for c in qc_checks if c["flag"])
    return {
        "template_len": len(template["full_template"]),
        "left_arm_len": template["left_arm_len"],
        "right_arm_len": template["right_arm_len"],
        "edit_seq": template["edit_seq"],
        "pam_note": pam_note,
        "qc_flags": n_flags,
        "fasta": str(fasta_path),
        "qc_report": str(qc_path),
        "arm_variants": str(variants_path),
        "plot": str(plot_path) if plot_path else None,
    }


def main() -> int:
    import argparse
    from .reference_reader import load_fasta, make_demo_reference

    parser = argparse.ArgumentParser(
        prog="hdr-template-designer",
        description="Design HDR donor templates for CRISPR edits.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--reference", type=Path, metavar="FASTA",
                        help="Reference sequence FASTA")
    parser.add_argument("--chrom", metavar="CHROM",
                        help="Chromosome/sequence name in reference")
    parser.add_argument("--cut-pos", type=int, metavar="INT",
                        help="Cut site position (0-based)")
    parser.add_argument("--edit-type", choices=["snv", "ins", "del"], default="snv",
                        help="Type of edit")
    parser.add_argument("--ref-allele", metavar="STR",
                        help="Reference allele (e.g. C for SNV)")
    parser.add_argument("--alt-allele", metavar="STR",
                        help="Alternate allele (e.g. T for C->T SNV)")
    parser.add_argument("--arm-len", type=int, default=100,
                        help="Homology arm length (bp each side)")
    parser.add_argument("--silent-pam", action="store_true",
                        help="Introduce silent PAM-disrupting mutation in right arm")
    parser.add_argument("--output-dir", type=Path, default=Path("results"),
                        metavar="PATH", help="Output directory")
    parser.add_argument("--no-plot", action="store_true", help="Skip diagram plot")
    parser.add_argument("--plot-format", choices=["png", "svg"], default="png")
    parser.add_argument("--demo", action="store_true",
                        help="Run demo: synthetic 500bp reference, C->T SNV at pos 250")
    args = parser.parse_args()

    if not args.demo:
        required = ["reference", "chrom", "cut_pos", "ref_allele", "alt_allele"]
        missing = [r for r in required if getattr(args, r) is None]
        if missing:
            parser.error(f"Required args missing: {', '.join('--' + r.replace('_', '-') for r in missing)}")

    if args.demo:
        chrom_name, ref_seq = make_demo_reference()
        cut_pos, edit_type, ref_allele, alt_allele = 250, "snv", "C", "T"
        silent_pam = True
        print(f"Demo mode: 500bp synthetic reference, C->T SNV at pos {cut_pos}, PAM disruption ON")
    else:
        genome = load_fasta(args.reference)
        if args.chrom not in genome:
            parser.error(f"Chromosome '{args.chrom}' not found in reference")
        ref_seq = genome[args.chrom]
        cut_pos = args.cut_pos
        edit_type = args.edit_type
        ref_allele = args.ref_allele
        alt_allele = args.alt_allele
        silent_pam = args.silent_pam

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
    print(f"PAM: {result['pam_note']}")
    print(f"QC flags: {result['qc_flags']}")
    print(f"  FASTA:          {result['fasta']}")
    print(f"  QC report:      {result['qc_report']}")
    print(f"  Arm variants:   {result['arm_variants']}")
    if result["plot"]:
        print(f"  Template diagram: {result['plot']}")
    return 0
