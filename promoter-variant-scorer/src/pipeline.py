"""Promoter Variant Scorer — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    vcf_path: Path | None,
    fasta_path: Path | None,
    pwm_db: Path | None,
    output_dir: Path,
    demo: bool = False,
    promoter_window: int = 2000,
    ic_threshold: float = 8.0,
    no_plot: bool = False,
) -> dict:
    """Score variants in promoter regions for transcription factor binding site disruption.

    Scans the promoter context of each variant against a JASPAR-format PWM database,
    computing reference and alternate allele scores to identify TFBS disruptions and
    de novo TFBS creations.

    Args:
        vcf_path: VCF of promoter variants. None if demo.
        fasta_path: Reference genome FASTA. None if demo.
        pwm_db: JASPAR-format PWM database (MEME or TRANSFAC). None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        promoter_window: Upstream window (bp) from TSS to define promoter region.
        ic_threshold: Minimum information content (bits) to include a PWM in the scan.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'disruptions', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Promoter Variant Scorer v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_scorer.py"
    spec = importlib.util.spec_from_file_location("run_scorer", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_scorer.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
