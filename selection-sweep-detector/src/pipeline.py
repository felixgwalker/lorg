"""Selection Sweep Detector — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    vcf_path: Path | None,
    population_map: Path | None,
    gene_annotation: Path | None,
    output_dir: Path,
    demo: bool = False,
    tests: list[str] | None = None,
    outlier_percentile: float = 99.0,
    no_plot: bool = False,
) -> dict:
    """Detect positive selection sweeps using haplotype and diversity statistics.

    Computes iHS (integrated haplotype score), XP-EHH (cross-population EHH),
    the composite likelihood ratio (CLR), and Tajima's D in sliding windows to
    identify genomic regions consistent with recent positive selection sweeps.

    Args:
        vcf_path: Phased multi-population VCF. None if demo.
        population_map: TSV mapping sample IDs to population labels. None if demo.
        gene_annotation: Gene annotation BED for candidate gene reporting. None optional.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        tests: Tests to run (iHS, XP-EHH, CLR, tajimas_d).
        outlier_percentile: Percentile threshold for outlier windows.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'windows', 'sweep_regions', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Selection Sweep Detector v%s ===", __version__)
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
