#!/usr/bin/env python3
"""CLI for the Codon Optimisation Engine.

Usage:
    python run_engine.py proteins.fa --host "Homo sapiens" --output-dir results/
    python run_engine.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_engine.py",
        description=(
            "Codon Optimisation Engine — optimises protein-coding sequences for "
            "expression in a target host by replacing codons using the host codon "
            "usage table, maximising CAI and avoiding restriction sites."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("proteins", metavar="FASTA", nargs="?", default=None,
                        help="Protein FASTA of sequences to optimise. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--host", default="Homo sapiens", metavar="ORGANISM",
                        help="Target host organism for codon usage.")
    parser.add_argument("--strategy",
                        choices=["most_frequent", "harmonised", "random_weighted", "CAI_maximised"],
                        default="CAI_maximised", help="Optimisation strategy.")
    parser.add_argument("--avoid-sites", nargs="*", metavar="ENZYME",
                        default=None, help="Restriction enzyme sites to avoid (e.g. EcoRI BamHI).")
    parser.add_argument("--gc-min", type=float, default=0.40, metavar="FLOAT",
                        help="Minimum target GC content.")
    parser.add_argument("--gc-max", type=float, default=0.65, metavar="FLOAT",
                        help="Maximum target GC content.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.proteins is None:
        _fatal("Provide a protein FASTA or use --demo.")
    if args.proteins and not Path(args.proteins).exists():
        _fatal(f"Protein FASTA not found: {args.proteins}")
    if args.gc_min >= args.gc_max:
        _fatal("--gc-min must be less than --gc-max.")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            protein_fasta=Path(args.proteins) if args.proteins else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            host_organism=args.host,
            strategy=args.strategy,
            avoid_restriction_sites=args.avoid_sites,
            target_gc_min=args.gc_min,
            target_gc_max=args.gc_max,
            no_plot=args.no_plot,
        )
    except NotImplementedError:
        print("ERROR: This tool is not yet implemented.", file=sys.stderr)
        return 1
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        logging.exception("Unhandled exception")
        return 1
    _print_summary(result)
    return 0


def _print_summary(result: dict) -> None:
    seqs = result.get("optimised_sequences", [])
    print(f"\nCodon Optimisation Engine v{result['pipeline_version']}")
    print(f"Sequences optimised: {len(seqs)}")
    out = result.get("output_files", {})
    print("Output files:")
    for k, v in out.items():
        print(f"  {k:20s}: {v}")
    print()


def _fatal(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
