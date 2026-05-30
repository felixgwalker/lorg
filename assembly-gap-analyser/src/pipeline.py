"""Assembly Gap Analyser — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    assembly_fasta: Path | None,
    gene_annotation: Path | None,
    reference_annotation: Path | None,
    output_dir: Path,
    demo: bool = False,
    min_gap_length: int = 10,
    no_plot: bool = False,
) -> dict:
    """Analyse gaps in a genome assembly, classifying their context and likely impact.

    Identifies all N-runs in the assembly, classifies gap type (contig/scaffold/
    centromere/telomere), and determines whether gaps interrupt gene models or
    fall between synteny-supported gene pairs.

    Args:
        assembly_fasta: Assembly FASTA. None if demo.
        gene_annotation: Gene annotation GTF for gap-gene overlap analysis. None optional.
        reference_annotation: Reference species annotation for synteny-gap analysis. None optional.
        output_dir: Directory for output BED/TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        min_gap_length: Minimum N run length to report.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'gaps', 'summary', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Assembly Gap Analyser v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_analyser.py"
    spec = importlib.util.spec_from_file_location("run_analyser", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_analyser.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
