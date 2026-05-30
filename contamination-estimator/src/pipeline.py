"""Contamination Estimator — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    bam_path: Path | None,
    reference_fasta: Path | None,
    output_dir: Path,
    demo: bool = False,
    methods: list[str] | None = None,
    contamination_threshold: float = 0.03,
    no_plot: bool = False,
) -> dict:
    """Estimate contamination rate in an ancient DNA BAM file.

    Uses mitochondrial consensus deviation, X-chromosome heterozygosity (for
    male samples), and ANGSD-based likelihood approaches to estimate the fraction
    of reads derived from modern human contamination.

    Args:
        bam_path: BAM file of ancient DNA reads. None if demo.
        reference_fasta: Reference genome FASTA. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        methods: Contamination methods to use (mt_consensus, nuclear_X, ANGSD, schmutzi).
        contamination_threshold: Contamination rate above which sample fails QC.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'estimates', 'combined_estimate', 'passes_threshold', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Contamination Estimator v%s ===", __version__)
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
