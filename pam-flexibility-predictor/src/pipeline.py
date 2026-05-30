"""PAM Flexibility Predictor — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    fasta_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    cas_variants: list[str] | None = None,
    no_plot: bool = False,
) -> dict:
    """Score PAM site availability for a panel of Cas variants at a target locus.

    Uses IUPAC position weight matrices for each Cas variant's PAM preference
    (e.g. NGG for SpCas9, NNGRRT for SaCas9, TTTV for AsCas12a) to enumerate
    and score all candidate PAM sites in both strands of the input sequence.
    Returns a per-variant density table and a ranked compatibility matrix.

    Args:
        fasta_path: Target locus FASTA. None if demo.
        output_dir: Directory for compatibility TSV and optional bar chart.
        demo: Use a synthetic 500 bp sequence.
        cas_variants: Subset of variants to score (default: full database).
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'pam_scores' (dict), 'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== PAM Flexibility Predictor v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_predictor.py"
    spec = importlib.util.spec_from_file_location("run_predictor", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_predictor.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
