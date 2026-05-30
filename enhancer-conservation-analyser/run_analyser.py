#!/usr/bin/env python3
"""CLI for the Enhancer Conservation Analyser.

Usage:
    python run_analyser.py --enhancers enhancers.bed --conservation phastcons.bw --output-dir results/
    python run_analyser.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_analyser.py",
        description=(
            "Enhancer Conservation Analyser — retrieves phastCons/PhyloP scores and "
            "multiple alignment coverage for enhancer elements, classifying each as "
            "highly conserved, moderately conserved, or lineage-specific."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--enhancers", metavar="BED", default=None,
                        help="BED of enhancer elements.")
    parser.add_argument("--conservation", metavar="BIGWIG", default=None,
                        help="PhastCons or PhyloP bigWig file.")
    parser.add_argument("--maf", metavar="MAF", default=None,
                        help="Multiple alignment MAF file.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--min-phastcons", type=float, default=0.4, metavar="FLOAT",
                        help="Minimum mean phastCons to call moderate conservation.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.enhancers is None:
        _fatal("Provide --enhancers or use --demo.")
    for p in [args.enhancers, args.conservation, args.maf]:
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
            enhancers=Path(args.enhancers) if args.enhancers else None,
            conservation_bigwig=Path(args.conservation) if args.conservation else None,
            alignment_maf=Path(args.maf) if args.maf else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            min_phastcons=args.min_phastcons,
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
    conserved = result.get("conserved_enhancers", [])
    print(f"\nEnhancer Conservation Analyser v{result['pipeline_version']}")
    print(f"Elements assessed: {result.get('n_elements_assessed', 'N/A')}")
    print(f"Elements classified: {len(conserved)}")
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
