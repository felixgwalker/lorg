"""Synthetic Promoter Designer — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    tf_list: Path | None,
    pwm_db: Path | None,
    output_dir: Path,
    demo: bool = False,
    promoter_type: str = "constitutive",
    n_designs: int = 10,
    promoter_length: int = 200,
    no_plot: bool = False,
) -> dict:
    """Design synthetic promoter sequences with configurable TFBS arrangements.

    Assembles synthetic promoter sequences from core promoter elements (TATA,
    Inr) and user-specified TF binding sites, scoring predicted strength based
    on TFBS density, spacing, and information content.

    Args:
        tf_list: Text file of TF names to include as binding sites. None if demo.
        pwm_db: JASPAR/MEME PWM database for TFBS motifs. None if demo.
        output_dir: Directory for output FASTA/TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        promoter_type: Promoter class (constitutive, inducible, tissue_specific).
        n_designs: Number of distinct promoter designs to generate.
        promoter_length: Total promoter length in base pairs.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'designs', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Synthetic Promoter Designer v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_designer.py"
    spec = importlib.util.spec_from_file_location("run_designer", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_designer.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
