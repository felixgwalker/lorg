"""Splice Impact Predictor — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    vcf_path: Path | None,
    fasta_path: Path | None,
    gtf_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    window: int = 200,
    delta_threshold: float = 2.0,
    canonical_only: bool = False,
    no_plot: bool = False,
) -> dict:
    """Predict the impact of variants on splicing.

    Scores each variant against MaxEntScan-style position weight matrices for
    donor and acceptor sites, reporting the delta score (alt − ref) and
    classifying the effect as disruption, creation, weakening, or neutral.

    Args:
        vcf_path: VCF of variants to evaluate. None if demo.
        fasta_path: Reference genome FASTA. None if demo.
        gtf_path: Gene annotation GTF. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        window: Bases around each variant to search for splice sites.
        delta_threshold: Absolute delta score to call a significant effect.
        canonical_only: Restrict analysis to canonical (GT-AG) splice sites.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'splice_scores', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Splice Impact Predictor v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_predictor.py"
    spec = importlib.util.spec_from_file_location("run_predictor", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_predictor.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
