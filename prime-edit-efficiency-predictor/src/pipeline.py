"""Prime Edit Efficiency Predictor — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    pegrna_path: Path | None,
    target_fasta: Path | None,
    output_dir: Path,
    demo: bool = False,
    no_plot: bool = False,
) -> dict:
    """Predict prime editing efficiency from pegRNA features.

    Extracts DeepPrime-style features: PBS GC content, RT template length, RT
    GC content, nick distance, spacer melting temperature, and MFE approximation
    of the RT template secondary structure.  Scores are normalised to 0–1 using
    logistic regression weights derived from published PE efficiency data.

    Args:
        pegrna_path: TSV or JSON with pegRNA designs (spacer, PBS, RT template).
        target_fasta: Target locus FASTA for context features. None if demo.
        output_dir: Directory for scored TSV and optional feature plot.
        demo: Use synthetic pegRNA designs.
        no_plot: Skip matplotlib figure generation.

    Returns:
        dict with keys 'predictions' (list), 'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Prime Edit Efficiency Predictor v%s ===", __version__)
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
