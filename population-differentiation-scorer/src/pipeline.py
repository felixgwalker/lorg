"""Population Differentiation Scorer — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    vcf_path: Path | None,
    population_map: Path | None,
    output_dir: Path,
    demo: bool = False,
    metrics: list[str] | None = None,
    window_size: int = 50000,
    outlier_percentile: float = 99.0,
    no_plot: bool = False,
) -> dict:
    """Score population differentiation using Fst and related metrics.

    Computes Weir-Cockerham Fst, Gst, and Jost's D between all population pairs
    genome-wide and in sliding windows, flagging Fst outlier windows as candidates
    for divergent selection.

    Args:
        vcf_path: Multi-population VCF. None if demo.
        population_map: TSV mapping sample IDs to population labels. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        metrics: Differentiation metrics to compute (Fst, Gst, Jost_D, Phi_st).
        window_size: Sliding window size in base pairs.
        outlier_percentile: Percentile threshold to flag Fst outlier windows.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'genome_wide', 'windows', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Population Differentiation Scorer v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_scorer.py"
    spec = importlib.util.spec_from_file_location("run_scorer", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_scorer.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
