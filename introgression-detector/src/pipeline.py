"""Introgression Detector — pipeline orchestrator."""

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
    p1: str = "pop1",
    p2: str = "pop2",
    p3: str = "pop3",
    outgroup: str = "outgroup",
    test: str = "D_statistic",
    window_kb: int = 50,
    no_plot: bool = False,
) -> dict:
    """Detect introgression between populations using ABBA-BABA statistics.

    Computes the D-statistic (Patterson's D) and f4-ratio across the genome
    in sliding windows to identify genomic regions showing evidence of gene
    flow between P3 and P1 or P2.

    Args:
        vcf_path: Multi-population VCF. None if demo.
        population_map: TSV mapping sample IDs to population labels. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        p1: Name of population P1 (no introgression expected).
        p2: Name of population P2 (no introgression expected).
        p3: Name of putative donor population.
        outgroup: Name of outgroup population.
        test: Statistical test to use (D_statistic, f4_ratio, RND_min, Dfoil).
        window_kb: Sliding window size in kilobases.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'genome_wide_result', 'segments', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Introgression Detector v%s ===", __version__)
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
