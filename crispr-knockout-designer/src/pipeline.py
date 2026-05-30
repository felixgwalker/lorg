"""CRISPR Knockout Designer — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    gene_fasta: Path | None,
    output_dir: Path,
    demo: bool = False,
    pam: str = "NGG",
    n_guides: int = 5,
    no_plot: bool = False,
) -> dict:
    """Design sgRNAs for CRISPR-mediated gene knockout.

    Targets early coding exons (exons 1–3 when annotation is provided, or the
    5′ third of the input sequence), scans both strands for NGG PAM sites,
    and scores each 20 nt spacer using Doench 2016-style on-target features
    (GC content, seed region composition, thermodynamic accessibility, position
    relative to start codon).  Returns the top-N guides ranked by on-target
    score with predicted frameshift efficiency.

    Args:
        gene_fasta: Gene or CDS FASTA. None if demo.
        output_dir: Directory for guide TSV and optional schematic plot.
        demo: Use synthetic gene sequence.
        pam: PAM motif (default NGG).
        n_guides: Number of top guides to return.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'guides' (list), 'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== CRISPR Knockout Designer v%s ===", __version__)
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
