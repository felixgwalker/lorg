"""Bottleneck Detector — pipeline orchestrator."""

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
    tests: list[str] | None = None,
    simulation_reps: int = 10000,
    no_plot: bool = False,
) -> dict:
    """Detect historical population bottlenecks from SNP data.

    Tests for bottleneck signatures using Tajima's D, the heterozygosity
    excess Wilcoxon (HEW) test, SFS mode-shift detection, and the M-ratio
    (allele count range to number of alleles), each sensitive to different
    bottleneck timescales.

    Args:
        vcf_path: VCF of SNPs for one population. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        tests: List of tests to run (tajimas_d, HEW, mode_shift, m_ratio).
        simulation_reps: Number of coalescent simulations for null distribution.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'sfs_shape', 'results', 'combined_signal', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Bottleneck Detector v%s ===", __version__)
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
