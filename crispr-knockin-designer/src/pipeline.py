"""CRISPR Knock-in Designer — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    locus_fasta: Path | None,
    insert_fasta: Path | None,
    output_dir: Path,
    demo: bool = False,
    pam: str = "NGG",
    arm_length: int = 100,
    max_cut_distance: int = 30,
    no_plot: bool = False,
) -> dict:
    """Design guide RNA and HDR donor template for CRISPR knock-in.

    Finds the PAM site(s) within `max_cut_distance` bp of the desired insertion
    point, selects the highest-scoring guide, and constructs an HDR donor
    with `arm_length` bp homology arms flanking the cut site.  Silent PAM
    mutations are introduced in the donor to prevent re-cutting after successful
    integration.  Outputs a guide TSV and donor FASTA suitable for synthesis.

    Args:
        locus_fasta: Target locus FASTA (≥500 bp). None if demo.
        insert_fasta: Insert sequence FASTA. None if demo.
        output_dir: Directory for guide TSV, donor FASTA, and optional schematic.
        demo: Use synthetic sequences.
        pam: PAM motif.
        arm_length: Homology arm length in bp.
        max_cut_distance: Max bp from insertion point to cut site.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'guide' (dict), 'donor_sequence' (str),
        'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== CRISPR Knock-in Designer v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_designer.py"
    spec = importlib.util.spec_from_file_location("run_designer", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_designer.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
