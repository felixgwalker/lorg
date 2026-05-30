"""CRISPR Array Designer — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    targets_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    cas_system: str = "Cas12a",
    direct_repeat: str | None = None,
    no_plot: bool = False,
) -> dict:
    """Design a CRISPR array encoding multiple spacers.

    For each target in the input, identifies a valid PAM site (TTTV for Cas12a,
    NGG for Cas9), extracts a 20–24 nt spacer, checks spacer uniqueness by
    exact k-mer matching within the input set, and assembles spacers into a
    final array sequence interspersed with direct repeats.

    Args:
        targets_path: FASTA or TSV of target sequences. None if demo.
        output_dir: Directory for array FASTA and spacer summary TSV.
        demo: Use synthetic target set.
        cas_system: Cas nuclease system ('Cas12a', 'Cas9', 'Cas12b').
        direct_repeat: Override direct repeat sequence (uses system default if None).
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'spacers' (list), 'array_sequence' (str),
        'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== CRISPR Array Designer v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_designer.py"
    spec = importlib.util.spec_from_file_location("run_designer", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_designer.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
