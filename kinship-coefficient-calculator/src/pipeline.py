"""Kinship Coefficient Calculator — pipeline orchestrator."""

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
    method: str = "KING",
    min_kinship: float = 0.0442,
    no_plot: bool = False,
) -> dict:
    """Calculate pairwise kinship coefficients between all samples.

    Computes the KING-robust kinship estimator (or genomic relatedness / IBD)
    for all pairs of samples in the VCF, classifying each pair by relationship
    (identical / 1st / 2nd / 3rd degree / unrelated).

    Args:
        vcf_path: Multi-sample VCF of LD-pruned SNPs. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        method: Kinship estimation method (KING, genomic_relatedness, IBD).
        min_kinship: Minimum kinship coefficient to report in output.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'kinship_pairs', 'n_samples', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Kinship Coefficient Calculator v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_calculator.py"
    spec = importlib.util.spec_from_file_location("run_calculator", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_calculator.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
