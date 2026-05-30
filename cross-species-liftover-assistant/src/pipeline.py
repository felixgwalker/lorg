"""Cross-Species Liftover Assistant — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    input_bed: Path | None,
    chain_file: Path | None,
    source_fasta: Path | None,
    target_fasta: Path | None,
    output_dir: Path,
    demo: bool = False,
    min_identity: float = 70.0,
    synteny_fallback: bool = True,
    no_plot: bool = False,
) -> dict:
    """Lift genomic coordinates from one species to another using a chain file or BLAST alignment.

    Maps BED intervals from a source genome to a target genome using pairwise
    alignment chains (liftOver-compatible), with a synteny-based fallback for
    coordinates in unconventional regions.

    Args:
        input_bed: BED of intervals to lift over. None if demo.
        chain_file: UCSC-format chain file for the source→target pair. None if demo.
        source_fasta: Source species FASTA (for BLAST fallback). Optional.
        target_fasta: Target species FASTA (for BLAST fallback). Optional.
        output_dir: Directory for output BED and optional plot.
        demo: Run on synthetic data without real inputs.
        min_identity: Minimum alignment identity to accept a liftover.
        synteny_fallback: Use synteny-based mapping if chain file fails.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'records', 'n_input', 'n_success', 'n_failed', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Cross-Species Liftover Assistant v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_assistant.py"
    spec = importlib.util.spec_from_file_location("run_assistant", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_assistant.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
