"""Contig Scaffolding Helper — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    contigs_fasta: Path | None,
    links_file: Path | None,
    output_dir: Path,
    demo: bool = False,
    evidence_type: str = "paired_reads",
    min_link_support: int = 3,
    no_plot: bool = False,
) -> dict:
    """Order and orient contigs into scaffolds using linking evidence.

    Uses paired-end read links, Hi-C contacts, or reference-guided scaffolding
    to determine contig order and orientation, joining them with estimated gap
    sizes.

    Args:
        contigs_fasta: Assembly contigs FASTA. None if demo.
        links_file: BAM or Hi-C contact file providing linking evidence. None if demo.
        output_dir: Directory for output FASTA and optional plot.
        demo: Run on synthetic data without real inputs.
        evidence_type: Scaffolding evidence (paired_reads, hi_c, reference_guided).
        min_link_support: Minimum links to join two contigs.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'scaffolds', 'unplaced_contigs', 'n_contigs_placed', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Contig Scaffolding Helper v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_helper.py"
    spec = importlib.util.spec_from_file_location("run_helper", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_helper.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
