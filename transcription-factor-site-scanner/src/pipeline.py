"""Transcription Factor Site Scanner — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    sequences: Path | None,
    pwm_db: Path | None,
    output_dir: Path,
    demo: bool = False,
    p_value_threshold: float = 1e-4,
    min_ic: float = 8.0,
    both_strands: bool = True,
    no_plot: bool = False,
) -> dict:
    """Scan sequences for transcription factor binding sites using PWMs.

    Scores each position in the input sequences against all PWMs in the
    database above the IC threshold, reporting hits with p-values calculated
    against a background nucleotide distribution.

    Args:
        sequences: FASTA of sequences to scan. None if demo.
        pwm_db: JASPAR/MEME-format PWM database. None if demo.
        output_dir: Directory for output BED/TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        p_value_threshold: P-value cutoff for reporting hits.
        min_ic: Minimum information content (bits) of PWM to include.
        both_strands: Scan both strands.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'hits', 'n_sequences_scanned', 'n_pwms_used', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Transcription Factor Site Scanner v%s ===", __version__)
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
