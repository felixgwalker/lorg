"""Safe Harbour Integration Finder — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    genome_fasta: Path | None,
    output_dir: Path,
    demo: bool = False,
    regulatory_bed: Path | None = None,
    oncogene_bed: Path | None = None,
    species: str = "human",
    no_plot: bool = False,
) -> dict:
    """Identify genomic safe harbour sites for transgene integration.

    Screens for homologs of validated safe harbour loci (AAVS1/PPP1R12C,
    H11/ROSA26, CCR5) by local alignment.  Scores candidate intergenic regions
    by: distance from the nearest oncogene (>1 Mb preferred), distance from
    regulatory elements (>50 kb preferred), open chromatin accessibility
    (proxy: low repeat density), and expression level of neighbouring genes.
    Returns a ranked BED file of integration candidates.

    Args:
        genome_fasta: Genome FASTA (or chromosomal subset). None if demo.
        output_dir: Directory for candidate BED, ranked TSV, and optional plot.
        demo: Use a synthetic 10 Mb chromosome fragment.
        regulatory_bed: BED of regulatory elements to avoid.
        oncogene_bed: BED of oncogene coordinates for distance scoring.
        species: Species identifier for known-harbour lookup ('human', 'mouse').
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'candidates' (list), 'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Safe Harbour Integration Finder v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_finder.py"
    spec = importlib.util.spec_from_file_location("run_finder", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_finder.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
