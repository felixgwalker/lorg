"""Founder Effect Estimator — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    study_vcf: Path | None,
    reference_vcf: Path | None,
    output_dir: Path,
    demo: bool = False,
    study_pop: str = "study",
    reference_pop: str = "reference",
    generation_time: float = 30.0,
    no_plot: bool = False,
) -> dict:
    """Estimate the strength of the founder effect in a study population.

    Compares nucleotide diversity, private variant fraction, and haplotype block
    length between the study and reference populations to infer whether a founder
    event occurred and estimate its magnitude and timing.

    Args:
        study_vcf: VCF of the study (potentially founded) population. None if demo.
        reference_vcf: VCF of the ancestral/reference population. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        study_pop: Label for the study population.
        reference_pop: Label for the reference population.
        generation_time: Assumed generation time in years.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'study_stats', 'reference_stats', 'estimate', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Founder Effect Estimator v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_estimator.py"
    spec = importlib.util.spec_from_file_location("run_estimator", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_estimator.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
