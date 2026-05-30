"""Cas Variant Selector — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    locus_fasta: Path | None,
    output_dir: Path,
    demo: bool = False,
    editing_goal: str = "knockout",
    no_plot: bool = False,
) -> dict:
    """Select the optimal Cas variant for a given target locus and editing goal.

    Scores a panel of Cas nucleases and base/prime editors (SpCas9, SaCas9,
    Cas9-NG, SpRY, AsCas12a, Cas12b, CasX, ABE8e, CBE4max, PE2) on: PAM site
    density at the locus, compatibility with the editing goal (DSB for KO,
    editing window position for base editors, RT template accessibility for PE),
    and practical constraints (protein size for AAV delivery, immunogenicity).

    Args:
        locus_fasta: Target locus FASTA. None if demo.
        output_dir: Directory for variant ranking TSV and optional bar chart.
        demo: Use a synthetic 300 bp locus.
        editing_goal: Desired edit type ('knockout', 'base-edit', 'prime-edit',
            'activation', 'repression').
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'rankings' (list), 'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Cas Variant Selector v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_selector.py"
    spec = importlib.util.spec_from_file_location("run_selector", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_selector.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
