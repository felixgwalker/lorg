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

Named-argument form:
    python run_visualiser.py \\
        --guide GCACTGACCTGAGTTCAGTG \\
        --sequence GCACTGACCTGAGTTCAGTGNGG \\
        --editor-type ABE8e \\
        --output-dir results/

Custom window:
    python run_visualiser.py \\
        ATCGATCGATCGATCGATCG \\
        ATCGATCGATCGATCGATCGNGG \\
        custom \\
        --editor-class ABE \\
        --window-start 3 \\
        --window-end 8

Demo mode (hardcoded example, no sequences required):
    python run_visualiser.py --demo
    python run_visualiser.py --demo --output-dir results/demo/

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
            "Built-in editors: ABE7.10, ABE8e, BE3, BE4max, evoAPOBEC, AncBE4max, CBE-NG, ABE-NG\n"
            "  Standard editors (ABE/CBE): require NGG PAM; window positions 4-8.\n"
            "  NG-variant editors (CBE-NG/ABE-NG): accept relaxed NG PAM; window 4-8.\n"
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
            "Omit when using --list-editors or --demo.  "
            "Alternatively use --guide."
        ),
    )
    parser.add_argument(
        "target_dna",
        metavar="TARGET",
        nargs="?",
        help=(
            "Target DNA non-template strand, 5'→3', including PAM "
            "(e.g. 20-nt protospacer + NGG).  Minimum length: 20 nt.  "
            "Alternatively use --sequence."
        ),
    )
    parser.add_argument(
        "editor",
        metavar="EDITOR",
        nargs="?",
        help=(
            "Base editor identifier.  Supported: ABE7.10, ABE8e, BE3, BE4max, "
            "evoAPOBEC, AncBE4max, CBE-NG, ABE-NG, or 'custom'.  "
            "Alternatively use --editor-type."
        ),
    )

    # ── Named sequence/editor arguments (alternative to positional) ───────
    parser.add_argument(
        "--guide",
        metavar="SEQ",
        dest="guide_named",
        default=None,
        help=(
            "20-nt guide RNA protospacer, 5'→3'.  "
            "Alternative to the GUIDE positional argument."
        ),
    )
    parser.add_argument(
        "--sequence",
        metavar="SEQ",
        dest="sequence_named",
        default=None,
        help=(
            "Target DNA sequence including PAM, 5'→3'.  "
            "Alternative to the TARGET positional argument."
        ),
    )
    parser.add_argument(
        "--editor-type",
        metavar="EDITOR",
        dest="editor_type_named",
        default=None,
        help=(
            "Base editor identifier.  "
            "Alternative to the EDITOR positional argument.  "
            "Supported: ABE7.10, ABE8e, BE3, BE4max, "
            "evoAPOBEC, AncBE4max, CBE-NG, ABE-NG, custom."
        ),
    )

    # ── Output ───────────────────────────────────────────────────────────
    parser.add_argument(
        "-o", "--output",
        metavar="PREFIX",
        default=None,
        dest="output_prefix",
        help=(
            "Output file prefix, including directory (default: be_output).  "
            "Files written: "
            "<prefix>_duplex.png, <prefix>_duplex.svg, <prefix>_duplex.html, "
            "<prefix>_editability.csv, <prefix>_outcomes.csv, "
            "<prefix>_bystander_warnings.txt, <prefix>_target_summary.tsv, "
            "<prefix>_window_coords.tsv, <prefix>_bystander_predictions.tsv."
        ),
    )
    parser.add_argument(
        "--output-dir",
        metavar="DIR",
        dest="output_dir",
        default=None,
        help=(
            "Output directory.  When provided, the output prefix is set to "
            "<DIR>/be_output (or <DIR>/demo_output in demo mode).  "
            "Takes precedence over -o/--output when both are given."
        ),
    )

    # ── Demo mode ─────────────────────────────────────────────────────────
    parser.add_argument(
        "--demo",
        action="store_true",
        help=(
            "Run in demo mode using a hardcoded example guide and target "
            "sequence (CBE BE4max, guide with a C at position 5).  "
            "No GUIDE/TARGET/EDITOR arguments are required."
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

def merge_named_args(args: argparse.Namespace) -> None:
    """
    Merge --guide/--sequence/--editor-type named args into the positional slots
    if the positional slots are empty.  Named args take precedence when both
    positional and named forms are provided.

    Also resolves --output-dir into output_prefix.
    """
    if args.guide_named:
        args.guide_rna = args.guide_named
    if args.sequence_named:
        args.target_dna = args.sequence_named
    if args.editor_type_named:
        args.editor = args.editor_type_named

    # Resolve output_prefix from --output-dir / -o / default
    if args.output_dir:
        import os
        stem = "demo_output" if getattr(args, "demo", False) else "be_output"
        args.output_prefix = os.path.join(args.output_dir, stem)
    elif args.output_prefix is None:
        args.output_prefix = "be_output"


def validate_args(args: argparse.Namespace) -> None:
    if args.list_editors or getattr(args, "demo", False):
        return

    if not args.guide_rna:
        _fatal("GUIDE is required.  Pass --list-editors to see available editors, or --demo to run a demo.")
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
        "duplex_png":               "  Duplex PNG              ->",
        "duplex_svg":               "  Duplex SVG              ->",
        "duplex_html":              "  Duplex HTML             ->",
        "editability_csv":          "  Editability CSV         ->",
        "outcomes_csv":             "  Outcomes CSV            ->",
        "bystander_txt":            "  Bystander TXT           ->",
        "target_summary_tsv":       "  Target summary TSV      ->",
        "window_coords_tsv":        "  Window coords TSV       ->",
        "bystander_predictions_tsv":"  Bystander pred. TSV     ->",
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

    # Merge --guide/--sequence/--editor-type into positional slots and
    # resolve --output-dir into output_prefix.
    merge_named_args(args)

    if args.list_editors:
        from src.config import EDITOR_PROFILES
        print("\nBuilt-in base editor profiles:\n")
        header = f"  {'Name':<12}  {'Class':<5}  {'Window':<8}  {'PAM':<5}  Description"
        print(header)
        print("  " + "-" * (len(header) - 2))
        from src.config import get_pam_requirement
        for name, profile in EDITOR_PROFILES.items():
            pam_req = get_pam_requirement(name)
            print(
                f"  {name:<12}  {profile.editor_class:<5}  "
                f"{profile.window_start}–{profile.window_end:<6}  "
                f"{pam_req:<5}  {profile.description}"
            )
        print(
            "\n  custom        ABE/CBE  user-defined  NGG    "
            "requires --editor-class --window-start --window-end\n"
        )
        return 0

    if getattr(args, "demo", False):
        return _run_demo(args)

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


# ---------------------------------------------------------------------------
# Demo mode
# ---------------------------------------------------------------------------

# Hardcoded demo sequences.
#
# CBE demo (BE4max, C→T):
#   Guide: AAAACGCAAGTGCGTTCAGT
#   Positions (1-indexed, PAM-distal):
#     pos 5 = C  (primary target, highest efficiency in window 4-8)
#     pos 7 = C  (bystander)
#   Window 4–8 contains: A(4), C(5), G(6), C(7), A(8)
#
# ABE demo (ABE8e, A→G):
#   Guide: TTTTAGTAACGCTTACAGTG
#   Positions:
#     pos 5 = A  (primary target)
#     pos 7 = A  (bystander)
#   Window 4–8 contains: A(4), A(5), G(6), A(7), C(8)

_DEMO_GUIDE  = "AAAACGCAAGTGCGTTCAGT"   # C at pos 5 (primary) and pos 7 (bystander)
_DEMO_TARGET = "AAAACGCAAGTGCGTTCAGTNGG"
_DEMO_EDITOR = "BE4max"

_DEMO_GUIDE_ABE  = "TTTAAGTAACGCTTACAGTG"  # A at pos 4,5,8 — primary=5, bystanders=4,8
_DEMO_TARGET_ABE = "TTTAAGTAACGCTTACAGTGNGG"
_DEMO_EDITOR_ABE = "ABE8e"


def _run_demo(args: argparse.Namespace) -> int:
    """Run the full pipeline on a hardcoded demo example and print results."""
    from src.config import EDITOR_PROFILES
    from src.pipeline import run_pipeline

    print()
    print("=" * 64)
    print("  CRISPR Base Editor Window Visualiser — DEMO MODE")
    print("=" * 64)
    print()
    print("  Demo 1: CBE (BE4max)  -- C->T editing")
    print(f"    Guide  : 5'-{_DEMO_GUIDE}-3'")
    print(f"    Target : 5'-{_DEMO_TARGET}-3'")
    print(f"    Editor : {_DEMO_EDITOR}")
    print()

    output_prefix = args.output_prefix if args.output_prefix else "demo_output"
    output_prefix_cbe = output_prefix + "_cbe"
    output_prefix_abe = output_prefix + "_abe"

    editor_cbe = EDITOR_PROFILES[_DEMO_EDITOR]
    editor_abe = EDITOR_PROFILES[_DEMO_EDITOR_ABE]

    try:
        result_cbe = run_pipeline(
            guide_rna=_DEMO_GUIDE,
            target_dna=_DEMO_TARGET,
            editor=editor_cbe,
            output_prefix=output_prefix_cbe,
            no_diagram=args.no_diagram,
        )
        print_summary(result_cbe, quiet=False)
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR in CBE demo: {exc}", file=sys.stderr)
        return 1

    print()
    print("  Demo 2: ABE (ABE8e)  -- A->G editing")
    print(f"    Guide  : 5'-{_DEMO_GUIDE_ABE}-3'")
    print(f"    Target : 5'-{_DEMO_TARGET_ABE}-3'")
    print(f"    Editor : {_DEMO_EDITOR_ABE}")
    print(f"    (A at positions 4, 5, 8 in window 4-8; pos 5 = primary)")
    print()

    try:
        result_abe = run_pipeline(
            guide_rna=_DEMO_GUIDE_ABE,
            target_dna=_DEMO_TARGET_ABE,
            editor=editor_abe,
            output_prefix=output_prefix_abe,
            no_diagram=args.no_diagram,
        )
        print_summary(result_abe, quiet=False)
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR in ABE demo: {exc}", file=sys.stderr)
        return 1

    return 0


def _fatal(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
