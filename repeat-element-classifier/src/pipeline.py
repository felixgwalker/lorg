"""Repeat Element Classifier — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    assembly_fasta: Path | None,
    repeat_masker_out: Path | None,
    output_dir: Path,
    demo: bool = False,
    min_length: int = 50,
    no_plot: bool = False,
) -> dict:
    """Classify and summarise repeat elements in a genome assembly.

    Parses RepeatMasker output (or runs de novo repeat identification) to
    classify repeats by class and family, compute the repeat landscape by
    divergence, and summarise the repeat content of the assembly.

    Args:
        assembly_fasta: Assembly FASTA to annotate. None if demo.
        repeat_masker_out: Pre-existing RepeatMasker .out file. None if not available.
        output_dir: Directory for output BED/TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        min_length: Minimum element length to include.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'elements', 'landscape', 'total_repeat_fraction', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Repeat Element Classifier v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_classifier.py"
    spec = importlib.util.spec_from_file_location("run_classifier", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_classifier.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
