#!/usr/bin/env python3
"""
CLI orchestrator for the CRISPR Base Editor Window Visualiser.

Run from the project root (crispr-base-editor-window-visualiser/):

    python run_visualiser.py GUIDE TARGET EDITOR [options]

Basic usage:
    python run_visualiser.py \\
        GCACTGACCTGAGTTCAGTG \\
        GCACTGACCTGAGTTCAGTGNGG \\
        ABE8e

Full usage example:
    python run_visualiser.py \\
        GCACTGACCTGAGTTCAGTG \\
        GCACTGACCTGAGTTCAGTGNGG \\
        BE4max \\
        --output results/thylacine_guide1 \\
        --target-position 5 \\
        --bystander-threshold 0.05

Custom window:
    python run_visualiser.py \\
        ATCGATCGATCGATCGATCG \\
        ATCGATCGATCGATCGATCGNGG \\
        custom \\
        --editor-class ABE \\
        --window-start 3 \\
        --window-end 8

List all built-in editors:
    python run_visualiser.py --list-editors
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_visualiser.py",
        description=(
            "CRISPR Base Editor Window Visualiser — maps a base editor's "
            "activity window onto a guide RNA–target duplex and produces a "
            "colour-coded diagram, editability table, and bystander warnings."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Built-in editors: ABE7.10, ABE8e, BE3, BE4max, evoAPOBEC, AncBE4max\n"
            "Use 'custom' as EDITOR together with --editor-class, --window-start, "
            "--window-end for a non-standard window.\n\n"
            "Position numbering: 1 = PAM-distal (5' end of protospacer), "
            "20 = PAM-proximal."
        ),
    )

    # ── Positional arguments ──────────────────────────────────────────────
    parser.add_argument(
        "guide_rna",
        metavar="GUIDE",
        nargs="?",
        help=(
            "20-nt guide RNA protospacer sequence, 5'→3' "
            "(DNA or RNA alphabet; U and T are treated as equivalent).  "
            "Omit when using --list-editors."
        ),
    )
    parser.add_argument(
        "target_dna",
        metavar="TARGET",
        nargs="?",
        help=(
            "Target DNA non-template strand, 5'→3', including PAM "
            "(e.g. 20-nt protospacer + NGG).  Minimum length: 20 nt."
        ),
    )
    parser.add_argument(
        "editor",
        metavar="EDITOR",
        nargs="?",
        help=(
            "Base editor identifier.  Supported: ABE7.10, ABE8e, BE3, BE4max, "
            "evoAPOBEC, AncBE4max, or 'custom'."
        ),
    )

    # ── Output ───────────────────────────────────────────────────────────
    parser.add_argument(
        "-o", "--output",
        metavar="PREFIX",
        default="be_output",
        dest="output_prefix",
        help=(
            "Output file prefix, including directory (default: be_output).  "
            "Four files are written: "
            "<prefix>_duplex.png, <prefix>_duplex.svg, "
            "<prefix>_editability.csv, <prefix>_outcomes.csv, "
            "<prefix>_bystander_warnings.txt."
        ),
    )

    # ── Editor customisation ─────────────────────────────────────────────
    parser.add_argument(
        "--editor-class",
        choices=["ABE", "CBE"],
        metavar="CLASS",
        help="Editor class for 'custom' editor: ABE (A→G) or CBE (C→T).",
    )
    parser.add_argument(
        "--window-start",
        type=int,
        metavar="INT",
        help=(
            "Override or define window start position (1-indexed, PAM-distal).  "
            "Required for 'custom'; optional override for built-in editors."
        ),
    )
    parser.add_argument(
        "--window-end",
        type=int,
        metavar="INT",
        help=(
            "Override or define window end position (1-indexed, PAM-distal).  "
            "Required for 'custom'; optional override for built-in editors."
        ),
    )

    # ── Analysis parameters ──────────────────────────────────────────────
    parser.add_argument(
        "--target-position",
        type=int,
        metavar="INT",
        help=(
            "1-indexed protospacer position of the intended primary edit "
            "(PAM-distal = 1).  When omitted the highest-efficiency editable "
            "base in the window is selected automatically."
        ),
    )
    parser.add_argument(
        "--bystander-threshold",
        type=float,
        metavar="FLOAT",
        default=0.10,
        help=(
            "Absolute editing frequency (0–1) above which a bystander position "
            "is classified HIGH risk (default: 0.10 = 10%%)."
        ),
    )

    # ── Misc ──────────────────────────────────────────────────────────────
    parser.add_argument(
        "--no-diagram",
        action="store_true",
        help="Skip PNG/SVG diagram rendering (useful in headless environments).",
    )
    parser.add_argument(
        "--list-editors",
        action="store_true",
        help="Print all built-in base editor profiles and exit.",
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress console summary.  Output files are still written.",
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging verbosity (default: INFO).",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 1.0.0",
    )

    return parser


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_args(args: argparse.Namespace) -> None:
    if args.list_editors:
        return

    if not args.guide_rna:
        _fatal("GUIDE is required.  Pass --list-editors to see available editors.")
    if not args.target_dna:
        _fatal("TARGET is required.")
    if not args.editor:
        _fatal("EDITOR is required.")

    if len(args.guide_rna.replace("U", "T")) != 20:
        _fatal(
            f"GUIDE must be exactly 20 nt (protospacer only); "
            f"got {len(args.guide_rna)} nt."
        )
    if len(args.target_dna) < 20:
        _fatal(
            f"TARGET must be at least 20 nt; got {len(args.target_dna)} nt."
        )

    if args.editor.lower() == "custom":
        if not args.editor_class:
            _fatal("--editor-class (ABE or CBE) is required for 'custom' editor.")
        if args.window_start is None or args.window_end is None:
            _fatal("--window-start and --window-end are required for 'custom' editor.")

    if args.window_start is not None or args.window_end is not None:
        ws = args.window_start or 1
        we = args.window_end or 20
        if not (1 <= ws <= we <= 20):
            _fatal(
                f"Window positions must satisfy 1 ≤ window-start ≤ window-end ≤ 20; "
                f"got {ws}–{we}."
            )

    if args.target_position is not None and not (1 <= args.target_position <= 20):
        _fatal(f"--target-position must be 1–20; got {args.target_position}.")

    if not (0.0 < args.bystander_threshold < 1.0):
        _fatal(
            f"--bystander-threshold must be in (0.0, 1.0); "
            f"got {args.bystander_threshold}."
        )


# ---------------------------------------------------------------------------
# Editor resolution
# ---------------------------------------------------------------------------

def resolve_editor(args: argparse.Namespace):
    """Return a BaseEditorProfile from the parsed arguments."""
    from src.config import EDITOR_PROFILES, build_custom_profile, BaseEditorProfile
    import dataclasses

    name = args.editor.lower()

    if name == "custom":
        return build_custom_profile(
            editor_class=args.editor_class,
            window_start=args.window_start,
            window_end=args.window_end,
        )

    # Case-insensitive lookup
    key = next((k for k in EDITOR_PROFILES if k.lower() == name), None)
    if key is None:
        _fatal(
            f"Unknown editor '{args.editor}'.  "
            f"Supported: {', '.join(EDITOR_PROFILES)}, custom."
        )

    editor = EDITOR_PROFILES[key]

    # Optional window override for built-in editors
    if args.window_start is not None or args.window_end is not None:
        ws = args.window_start if args.window_start is not None else editor.window_start
        we = args.window_end   if args.window_end   is not None else editor.window_end
        from src.config import _gaussian_profile
        editor = dataclasses.replace(
            editor,
            window_start=ws,
            window_end=we,
            efficiency_profile=_gaussian_profile(
                (ws + we) / 2,
                max((we - ws) / 2.5, 0.5),
                ws, we,
            ),
            description=editor.description + f"  [window overridden: {ws}–{we}]",
        )

    return editor


# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------

def print_summary(result: dict, quiet: bool) -> None:
    if quiet:
        return

    positions   = result["positions"]
    editor      = result["editor"]
    pam_seq     = result["pam_seq"]
    out_files   = result["output_files"]
    by_lines    = result["bystander_warnings"]

    protospacer = "".join(p.protospacer_base for p in positions)
    guide_rna   = "".join(p.guide_base for p in positions)

    print()
    print("=" * 62)
    print(f"  CRISPR Base Editor Window Visualiser  —  {editor.name}")
    print("=" * 62)
    print(f"  Guide RNA  : 5'-{guide_rna}-3'")
    print(f"  Target DNA : 5'-{protospacer}[{pam_seq}]-3'")
    print(f"  Editor     : {editor.name}")
    print(f"  Window     : positions {editor.window_start}–{editor.window_end}")
    print(f"  Edit type  : {editor.target_base}->{editor.product_base}")
    print()

    primary = [p for p in positions if p.is_primary_target]
    bystand = [p for p in positions if p.is_bystander]

    if primary:
        p = primary[0]
        print(
            f"  Primary target  : position {p.position} ({p.protospacer_base})  "
            f"— predicted {p.absolute_efficiency * 100:.1f}% editing"
        )
    else:
        print(f"  Primary target  : none found in window")

    if bystand:
        print(f"  Bystander risks : {len(bystand)} position(s)")
        for p in bystand:
            pct = p.absolute_efficiency * 100
            print(f"      pos {p.position:2d}  {p.protospacer_base}  {pct:.1f}%")
    else:
        print("  Bystander risks : none")

    print()
    print("  Output files:")
    labels = {
        "duplex_png":      "  Duplex PNG         ->",
        "duplex_svg":      "  Duplex SVG         ->",
        "editability_csv": "  Editability CSV    ->",
        "outcomes_csv":    "  Outcomes CSV       ->",
        "bystander_txt":   "  Bystander TXT      ->",
    }
    for key, label in labels.items():
        if key in out_files:
            print(f"  {label} {out_files[key]}")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args   = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s: %(message)s",
        stream=sys.stderr,
    )

    if args.list_editors:
        from src.config import EDITOR_PROFILES
        print("\nBuilt-in base editor profiles:\n")
        header = f"  {'Name':<12}  {'Class':<5}  {'Window':<8}  Description"
        print(header)
        print("  " + "-" * (len(header) - 2))
        for name, profile in EDITOR_PROFILES.items():
            print(
                f"  {name:<12}  {profile.editor_class:<5}  "
                f"{profile.window_start}–{profile.window_end:<6}  "
                f"{profile.description}"
            )
        print(
            "\n  custom        ABE/CBE  user-defined  "
            "requires --editor-class --window-start --window-end\n"
        )
        return 0

    validate_args(args)
    editor = resolve_editor(args)

    from src.pipeline import run_pipeline

    try:
        result = run_pipeline(
            guide_rna=args.guide_rna,
            target_dna=args.target_dna,
            editor=editor,
            output_prefix=args.output_prefix,
            target_position=args.target_position,
            bystander_threshold=args.bystander_threshold,
            no_diagram=args.no_diagram,
        )
    except ValueError as exc:
        _fatal(str(exc))
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: Unexpected failure — {exc}", file=sys.stderr)
        logging.exception("Unhandled exception in pipeline")
        return 1

    print_summary(result, args.quiet)
    return 0


def _fatal(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
