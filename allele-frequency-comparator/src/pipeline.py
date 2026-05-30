"""Allele Frequency Comparator — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    vcf_path: Path | None,
    af_table: Path | None,
    output_dir: Path,
    demo: bool = False,
    populations: list[str] | None = None,
    fold_change_threshold: float = 5.0,
    min_af: float = 1e-5,
    no_plot: bool = False,
) -> dict:
    """Compare allele frequencies across gnomAD populations.

    For each variant, retrieves per-population allele counts and frequencies,
    computes pairwise fold changes, tests for population differentiation via
    Fisher's exact test, and estimates Fst.

    Args:
        vcf_path: VCF or variant list to query. None if demo.
        af_table: gnomAD-format population AF TSV. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        populations: Subset of gnomAD populations to compare.
        fold_change_threshold: Min fold change to flag a differential variant.
        min_af: Minimum AF in at least one population to include.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'comparisons', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Allele Frequency Comparator v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_comparator.py"
    spec = importlib.util.spec_from_file_location("run_comparator", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_comparator.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
