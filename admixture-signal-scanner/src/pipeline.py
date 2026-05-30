"""Admixture Signal Scanner — pipeline orchestrator."""

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
    k_min: int = 2,
    k_max: int = 10,
    model: str = "unsupervised",
    reference_panel: Path | None = None,
    no_plot: bool = False,
) -> dict:
    """Scan for admixture signals and estimate ancestry proportions.

    Fits a parametric admixture model for K = k_min … k_max using an
    EM algorithm on LD-pruned SNPs, selects the best K by cross-validation
    error, and reports per-sample ancestry proportions.

    Args:
        vcf_path: VCF of LD-pruned SNPs. None if demo.
        output_dir: Directory for output and optional plot.
        demo: Run on synthetic data without real inputs.
        k_min: Minimum number of ancestry components to test.
        k_max: Maximum number of ancestry components to test.
        model: Admixture model type (unsupervised, supervised).
        reference_panel: Reference population panel VCF for supervised mode.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'best_k', 'runs', 'sample_ancestries', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Admixture Signal Scanner v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_scanner.py"
    spec = importlib.util.spec_from_file_location("run_scanner", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_scanner.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
