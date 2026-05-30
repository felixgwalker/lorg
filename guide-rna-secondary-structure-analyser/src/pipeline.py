"""Guide RNA Secondary Structure Analyser — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    guides_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    scaffold: str = "SpCas9",
    no_plot: bool = False,
) -> dict:
    """Analyse secondary structure of guide RNA sequences.

    Concatenates each 20 nt spacer with the sgRNA scaffold sequence and
    predicts RNA secondary structure using a Zuker nearest-neighbour MFE
    approximation (Turner 2004 parameters).  Reports minimum free energy,
    the fraction of seed-region nucleotides (positions 1–12) that are
    unpaired (accessibility score), and a flag for guide–scaffold duplex
    formation that would occlude the spacer.

    Args:
        guides_path: FASTA or TSV of 20 nt spacer sequences. None if demo.
        output_dir: Directory for structure TSV and optional dot-bracket plot.
        demo: Use synthetic guide sequences.
        scaffold: sgRNA scaffold sequence identifier ('SpCas9', 'SaCas9', custom).
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'structures' (list), 'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Guide RNA Secondary Structure Analyser v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_analyser.py"
    spec = importlib.util.spec_from_file_location("run_analyser", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_analyser.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
