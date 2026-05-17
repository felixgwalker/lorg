"""
Pipeline orchestrator for the CRISPR Base Editor Window Visualiser.

Calls analysis, diagram, and writer modules in order and returns a structured
result dict.  Also exposes main() so that main.py can delegate to it.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional

from src import __version__
from src.analyser import analyse_sequence
from src.config import (
    DEFAULT_BYSTANDER_THRESHOLD,
    DEFAULT_OUTPUT_PREFIX,
    BaseEditorProfile,
    EDITOR_PROFILES,
    build_custom_profile,
)
from src.diagram import generate_duplex_diagram
from src.writer import (
    write_bystander_warnings,
    write_editability_csv,
    write_outcomes_csv,
)

logger = logging.getLogger(__name__)


def run_pipeline(
    guide_rna: str,
    target_dna: str,
    editor: BaseEditorProfile,
    output_prefix: str = DEFAULT_OUTPUT_PREFIX,
    target_position: Optional[int] = None,
    bystander_threshold: float = DEFAULT_BYSTANDER_THRESHOLD,
    no_diagram: bool = False,
) -> dict:
    """
    Execute the complete Base Editor Window Visualiser pipeline.

    Steps:
        1. Analyse    — map activity window onto duplex, classify positions.
        2. Diagram    — render colour-coded PNG + SVG duplex diagram.
        3. Editability CSV — per-position annotation table.
        4. Outcomes CSV    — predicted edit products and frequencies.
        5. Bystander TXT   — warning report for off-target edits.

    Args:
        guide_rna:            20-nt protospacer (5'→3', DNA or RNA alphabet).
        target_dna:           Non-template strand incl. PAM (5'→3', ≥ 20 nt).
        editor:               BaseEditorProfile describing window and efficiency.
        output_prefix:        Path prefix for all output files.
        target_position:      1-indexed position of the intended primary edit.
        bystander_threshold:  Absolute editing frequency flagging bystanders HIGH.
        no_diagram:           If True skip matplotlib rendering.

    Returns:
        Structured dict with keys:
            "positions"           → list[PositionInfo]
            "editor"              → BaseEditorProfile
            "pam_seq"             → str
            "output_files"        → dict[str, str]
            "bystander_warnings"  → list[str]
            "pipeline_version"    → str
    """
    Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)

    guide  = guide_rna.upper().replace("U", "T").strip()
    target = target_dna.upper().replace("U", "T").strip()
    pam_seq = target[20:23] if len(target) >= 23 else "NGG"

    logger.info("=== CRISPR Base Editor Window Visualiser v%s ===", __version__)
    logger.info("Editor  : %s  (%s)", editor.name, editor.description)
    logger.info("Window  : positions %d–%d (PAM-distal indexing)", editor.window_start, editor.window_end)
    logger.info("Guide   : 5'-%s-3'", guide)
    logger.info("Target  : 5'-%s[%s]-3'", target[:20], pam_seq)

    # Step 1 — Analyse
    logger.info("Step 1/4 — Analysing sequence...")
    positions = analyse_sequence(guide, target, editor, target_position=target_position)

    n_primary   = sum(1 for p in positions if p.is_primary_target)
    n_bystander = sum(1 for p in positions if p.is_bystander)
    logger.info(
        "  Primary targets: %d  |  Bystander risks: %d",
        n_primary, n_bystander,
    )

    output_files: dict[str, str] = {}

    # Step 2 — Diagram
    if no_diagram:
        logger.info("Step 2/4 — Diagram skipped (--no-diagram).")
    else:
        logger.info("Step 2/4 — Rendering duplex diagram...")
        try:
            png_path, svg_path = generate_duplex_diagram(
                positions, editor, pam_seq, output_prefix
            )
            output_files["duplex_png"] = str(png_path)
            output_files["duplex_svg"] = str(svg_path)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Diagram rendering failed: %s.  Continuing.", exc)

    # Step 3 — Editability CSV
    logger.info("Step 3/4 — Writing editability table...")
    edit_path = write_editability_csv(positions, output_prefix)
    output_files["editability_csv"] = str(edit_path)

    # Step 4 — Outcomes CSV + bystander report (combined step for user clarity)
    logger.info("Step 4/4 — Writing outcomes and bystander report...")
    out_path = write_outcomes_csv(positions, editor, output_prefix)
    output_files["outcomes_csv"] = str(out_path)

    by_path, by_lines = write_bystander_warnings(
        positions, editor, bystander_threshold, output_prefix
    )
    output_files["bystander_txt"] = str(by_path)

    logger.info("Done.  Output prefix: %s", output_prefix)

    return {
        "positions":          positions,
        "editor":             editor,
        "pam_seq":            pam_seq,
        "output_files":       output_files,
        "bystander_warnings": by_lines,
        "pipeline_version":   __version__,
    }


def main(argv: list[str] | None = None) -> int:
    """
    Entry point delegated to from main.py — imports and invokes run_visualiser.main().
    """
    import importlib.util

    rv_path = Path(__file__).parent.parent / "run_visualiser.py"
    spec = importlib.util.spec_from_file_location("run_visualiser", rv_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_visualiser.py at {rv_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
