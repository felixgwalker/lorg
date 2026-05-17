"""Command-line interface for hgt-risk-assessor."""

import argparse
import logging
import sys
from pathlib import Path

from src.models import InputFormat


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="hgt-risk-assessor",
        description=(
            "Sequence-level horizontal gene transfer risk scoring "
            "for engineered organisms."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument("--input", required=True, type=Path, metavar="FILE",
                   help="Input sequence file (FASTA or VCF)")
    p.add_argument("--host", required=True, metavar="HOST_ID",
                   help="Host organism (NCBI accession or common identifier)")
    p.add_argument("--input-format", choices=["fasta", "vcf"], default="fasta",
                   metavar="{fasta,vcf}",
                   help="Input file format")
    p.add_argument("--ref", type=Path, metavar="FASTA",
                   help="Reference FASTA (required when --input-format vcf)")
    p.add_argument("--host-gc", type=float, metavar="FLOAT",
                   help="Host GC content as a fraction 0.0–1.0 (skips Entrez lookup)")

    p.add_argument("--data-dir", type=Path, default=Path("data"), metavar="PATH",
                   help="Directory containing BLAST databases")
    p.add_argument("--output-dir", type=Path, default=Path("results"), metavar="PATH",
                   help="Directory for output reports")
    p.add_argument("--output-prefix", metavar="STR",
                   help="Prefix for output filenames (default: {query_id}_{timestamp})")

    p.add_argument("--threads", type=int, default=4,
                   help="BLAST thread count")
    p.add_argument("--no-network", action="store_true",
                   help="Disable PHASTER API and Entrez lookups (offline mode)")
    p.add_argument("--entrez-email", metavar="EMAIL",
                   help="NCBI Entrez email (required unless --host-gc or --no-network)")
    p.add_argument("--entrez-api-key", metavar="KEY",
                   help="NCBI Entrez API key (raises rate limit to 10 req/s)")
    p.add_argument("--blast-bin-dir", type=Path, metavar="PATH",
                   help="Directory containing BLAST+ binaries if not on PATH")
    p.add_argument("--skip-signal", nargs="*", default=[],
                   choices=["gc_content", "is_proximity", "integron",
                             "conjugative", "prophage"],
                   metavar="SIGNAL",
                   help="Signal(s) to skip")

    # Three-layer model options
    p.add_argument("--weight-profile", default="default",
                   choices=["default", "environmental", "clinical_amr"],
                   help="Weight profile for the three-layer HGT Risk Index")
    p.add_argument("--donor-taxon", metavar="TAXON",
                   help="Donor organism name for taxonomic distance feature (optional)")
    p.add_argument("--recipient-taxon", metavar="TAXON",
                   help="Recipient organism name (defaults to --host if omitted)")

    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    p.add_argument("--version", action="version", version="%(prog)s 0.2.0")

    return p


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
        datefmt="%H:%M:%S",
    )

    # Input validation
    if args.input_format == "vcf" and args.ref is None:
        parser.error("--ref is required when --input-format is vcf")

    if (not args.no_network
            and args.host_gc is None
            and not args.entrez_email):
        parser.error(
            "Provide one of: --host-gc FLOAT, --entrez-email EMAIL, or --no-network. "
            "Example: --host-gc 0.52"
        )

    if args.host_gc is not None and not (0.0 <= args.host_gc <= 1.0):
        parser.error("--host-gc must be in the range 0.0–1.0")

    from src.pipeline import run
    try:
        result = run(
            input_path=args.input,
            host_id=args.host,
            input_format=InputFormat(args.input_format),
            ref_path=args.ref,
            host_gc=args.host_gc,
            data_dir=args.data_dir,
            output_dir=args.output_dir,
            output_prefix=args.output_prefix,
            threads=args.threads,
            no_network=args.no_network,
            entrez_email=args.entrez_email or "",
            entrez_api_key=args.entrez_api_key or "",
            blast_bin_dir=args.blast_bin_dir,
            skip_signals=args.skip_signal or [],
            weight_profile=args.weight_profile,
            donor_taxon=args.donor_taxon,
            recipient_taxon=args.recipient_taxon,
        )
        agg = result.aggregation
        tl  = result.three_layer

        print(f"\n  ── Flat signal model ──────────────────────")
        print(f"  Risk Index : {agg.risk_index:.3f}  ({agg.risk_level.value})")
        if agg.skipped_signals:
            print(f"  Skipped    : {', '.join(agg.skipped_signals)}")

        if tl:
            print(f"\n  ── Three-layer model ({tl.weight_profile_name}) ─────")
            print(f"  HGT Risk Index : {tl.hgt_risk_index:.3f}  ({tl.score_band.value})")
            print(f"  Transfer       : {tl.transfer_layer.layer_score:.3f}")
            print(f"  Establishment  : {tl.establishment_layer.layer_score:.3f}")
            print(f"  Consequence    : {tl.consequence_layer.layer_score:.3f}")
            print(f"  Completeness   : {tl.overall_completeness:.0%}")
            print(f"\n  {tl.explanation}")
        print()
    except Exception as exc:
        logging.error(str(exc))
        sys.exit(1)


if __name__ == "__main__":
    main()
