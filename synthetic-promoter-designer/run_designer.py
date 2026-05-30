#!/usr/bin/env python3
"""CLI for the Synthetic Promoter Designer.

Usage:
    python run_designer.py --tf-list tfs.txt --pwm-db jaspar.meme --output-dir results/
    python run_designer.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_designer.py",
        description=(
            "Synthetic Promoter Designer — assembles synthetic promoter sequences "
            "from core elements and specified TFBS arrangements, scoring predicted "
            "strength and generating multiple design variants."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--tf-list", metavar="FILE", default=None,
                        help="Text file of TF names to include as binding sites.")
    parser.add_argument("--pwm-db", metavar="MEME/JASPAR", default=None,
                        help="JASPAR or MEME-format PWM database.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--promoter-type",
                        choices=["constitutive", "inducible", "tissue_specific",
                                 "cell_cycle_regulated"],
                        default="constitutive", help="Promoter class.")
    parser.add_argument("--n-designs", type=int, default=10, metavar="INT",
                        help="Number of distinct promoter designs to generate.")
    parser.add_argument("--promoter-length", type=int, default=200, metavar="INT",
                        help="Total promoter length in base pairs.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.tf_list is None:
        _fatal("Provide --tf-list or use --demo.")
    for p in [args.tf_list, args.pwm_db]:
        if p and not Path(p).exists():
            _fatal(f"File not found: {p}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            tf_list=Path(args.tf_list) if args.tf_list else None,
            pwm_db=Path(args.pwm_db) if args.pwm_db else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            promoter_type=args.promoter_type,
            n_designs=args.n_designs,
            promoter_length=args.promoter_length,
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
    designs = result.get("designs", [])
    print(f"\nSynthetic Promoter Designer v{result['pipeline_version']}")
    print(f"Promoter designs generated: {len(designs)}")
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
