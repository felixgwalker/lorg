"""Genomic Diversity Index — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    vcf_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    window_size: int = 50000,
    step_size: int = 10000,
    metrics: list[str] | None = None,
    no_plot: bool = False,
) -> dict:
    """Compute a comprehensive genomic diversity index across the genome.

    Calculates per-window θW, θπ, Tajima's D, observed heterozygosity (Ho),
    expected heterozygosity (He), and inbreeding coefficient (Fis) in sliding
    windows, reporting genome-wide summary statistics.

    Args:
        vcf_path: VCF of SNPs (mono- or multi-sample). None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        window_size: Sliding window size in base pairs.
        step_size: Step size in base pairs.
        metrics: Subset of metrics to compute (theta_w, theta_pi, tajimas_d, Ho, He, Fis).
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'population_index', 'windows', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Genomic Diversity Index v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_index.py"
    spec = importlib.util.spec_from_file_location("run_index", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_index.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
