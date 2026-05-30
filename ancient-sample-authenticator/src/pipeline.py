"""Ancient Sample Authenticator — pipeline orchestrator."""

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
    sample_id: str = "sample",
    no_plot: bool = False,
) -> dict:
    """Authenticate an ancient DNA sample against multiple QC criteria.

    Evaluates fragment length, terminal deamination rate, contamination estimate,
    endogenous DNA fraction, and coverage, combining them into a composite
    authentication score and verdict.

    Args:
        bam_path: BAM file of aligned ancient DNA reads. None if demo.
        reference_fasta: Reference genome FASTA. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        sample_id: Sample identifier for report.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'result', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Ancient Sample Authenticator v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_authenticator.py"
    spec = importlib.util.spec_from_file_location("run_authenticator", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_authenticator.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
