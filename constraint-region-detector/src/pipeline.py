"""Constraint Region Detector — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    vcf_path: Path | None,
    constraint_file: Path | None,
    output_dir: Path,
    demo: bool = False,
    loeuf_threshold: float = 0.35,
    z_score_threshold: float = 3.09,
    primary_metric: str = "LOEUF",
    no_plot: bool = False,
) -> dict:
    """Detect whether variants fall in genomically constrained regions.

    Intersects a VCF with a gnomAD-style constraint annotation table, flagging
    variants that fall in genes or regions with high constraint scores (low LOEUF,
    high pLI, or high missense Z-score).

    Args:
        vcf_path: VCF of variants to evaluate. None if demo.
        constraint_file: gnomAD-format constraint TSV. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        loeuf_threshold: LOEUF upper bound; genes below this are constrained.
        z_score_threshold: Missense Z-score above which a gene is constrained.
        primary_metric: Metric to use for primary constraint call (LOEUF, pLI, z_score).
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'overlaps', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Constraint Region Detector v%s ===", __version__)
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
