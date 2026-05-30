#!/usr/bin/env python3
"""CLI for the Ancestral State Reconstructor.

Usage:
    python run_reconstructor.py alignment.fa --phylogeny tree.nwk --output-dir results/
    python run_reconstructor.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_reconstructor.py",
        description=(
            "Ancestral State Reconstructor — reconstructs ancestral nucleotide or amino "
            "acid sequences at internal phylogenetic nodes using marginal ML, parsimony, "
            "or Bayesian inference."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("alignment", metavar="FASTA", nargs="?", default=None,
                        help="FASTA multiple sequence alignment. Omit with --demo.")
    parser.add_argument("--phylogeny", metavar="NWK", default=None,
                        help="Newick phylogenetic tree.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--method", choices=["parsimony", "maximum_likelihood", "bayesian"],
                        default="maximum_likelihood", help="Reconstruction method.")
    parser.add_argument("--substitution-model", default="GTR+G", metavar="MODEL",
                        help="Substitution model for ML reconstruction (e.g. GTR+G, JC).")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.alignment is None:
        _fatal("Provide an alignment FASTA or use --demo.")
    if args.alignment and not Path(args.alignment).exists():
        _fatal(f"Alignment not found: {args.alignment}")
    if args.phylogeny and not Path(args.phylogeny).exists():
        _fatal(f"Phylogeny not found: {args.phylogeny}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            alignment=Path(args.alignment) if args.alignment else None,
            phylogeny=Path(args.phylogeny) if args.phylogeny else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            method=args.method,
            substitution_model=args.substitution_model,
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
    seqs = result.get("ancestral_sequences", [])
    print(f"\nAncestral State Reconstructor v{result['pipeline_version']}")
    print(f"Ancestral sequences reconstructed: {len(seqs)}")
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
