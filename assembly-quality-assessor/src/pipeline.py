"""Assembly Quality Assessor — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    assembly_fasta: Path | None,
    output_dir: Path,
    demo: bool = False,
    no_plot: bool = False,
) -> dict:
    """Assess genome assembly quality from a FASTA file.

    Computes N50, N90, L50, L90, total length, GC content, gap content,
    and ambiguous base count, classifying the assembly as reference-quality,
    chromosome-level, scaffold-level, or contig-level.

    Args:
        assembly_fasta: Assembly FASTA file. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'stats', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Assembly Quality Assessor v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_assessor.py"
    spec = importlib.util.spec_from_file_location("run_assessor", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_assessor.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
