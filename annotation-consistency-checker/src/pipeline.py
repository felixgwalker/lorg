"""Annotation Consistency Checker — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    annotation: Path | None,
    fasta: Path | None,
    output_dir: Path,
    demo: bool = False,
    strict: bool = False,
    no_plot: bool = False,
) -> dict:
    """Check a gene annotation GTF for internal consistency and structural errors.

    Validates coordinate hierarchies, detects duplicate IDs, checks for
    overlapping features on the same strand, and verifies that chromosome names
    match the reference FASTA.

    Args:
        annotation: GTF annotation to check. None if demo.
        fasta: Reference genome FASTA for chromosome name validation. None optional.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        strict: Fail on any issue (rather than warning).
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'issues', 'n_genes_checked', 'n_transcripts_checked', 'passes_qc',
        'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Annotation Consistency Checker v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_checker.py"
    spec = importlib.util.spec_from_file_location("run_checker", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_checker.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
