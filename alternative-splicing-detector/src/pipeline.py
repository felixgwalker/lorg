"""Alternative Splicing Detector — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    bam_a: Path | None,
    bam_b: Path | None,
    annotation: Path | None,
    output_dir: Path,
    demo: bool = False,
    event_types: list[str] | None = None,
    min_delta_psi: float = 0.1,
    no_plot: bool = False,
) -> dict:
    """Detect differential alternative splicing events between two conditions.

    Quantifies percent spliced in (PSI) for exon skipping, intron retention,
    and alternative donor/acceptor events from RNA-seq BAMs, testing for
    significant changes between conditions.

    Args:
        bam_a: RNA-seq BAM for condition A. None if demo.
        bam_b: RNA-seq BAM for condition B. None if demo.
        annotation: GTF gene annotation. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        event_types: Event types to detect (exon_skipping, intron_retention, etc.).
        min_delta_psi: Minimum |ΔPSI| to report as differential.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'events', 'n_genes_with_events', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Alternative Splicing Detector v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_detector.py"
    spec = importlib.util.spec_from_file_location("run_detector", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_detector.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
