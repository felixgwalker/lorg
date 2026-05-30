"""Genome Completeness Estimator — pipeline orchestrator."""

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
    lineage: str = "vertebrata_odb10",
    mode: str = "genome",
    no_plot: bool = False,
) -> dict:
    """Estimate genome assembly completeness using BUSCO conserved gene benchmarks.

    Searches for a set of BUSCO conserved orthologs from a specified lineage
    database in the assembly, classifying each as complete (single/duplicated),
    fragmented, or missing.

    Args:
        assembly_fasta: Assembly FASTA (or protein FASTA for --mode proteins). None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        lineage: BUSCO lineage dataset (e.g. vertebrata_odb10, metazoa_odb10).
        mode: Analysis mode (genome, proteins, transcriptome).
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'stats', 'busco_results', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Genome Completeness Estimator v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_estimator.py"
    spec = importlib.util.spec_from_file_location("run_estimator", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_estimator.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
