"""Palaeogenomic Coverage Assessor — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    bam_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    min_mapq: int = 25,
    remove_duplicates: bool = True,
    no_plot: bool = False,
) -> dict:
    """Assess genome-wide coverage statistics for a palaeogenomic BAM.

    Computes mapping rate, mean/median depth, breadth of coverage at 1×/5×/10×,
    per-chromosome coverage, endogenous fraction, and duplication rate to
    characterise the quality and utility of an ancient genome dataset.

    Args:
        bam_path: BAM file of aligned ancient DNA reads. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        min_mapq: Minimum mapping quality for depth calculation.
        remove_duplicates: Exclude duplicate reads from coverage calculation.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'assessment', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Palaeogenomic Coverage Assessor v%s ===", __version__)
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
