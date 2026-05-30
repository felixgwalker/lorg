"""Gene Model Validator — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    annotation: Path | None,
    fasta: Path | None,
    output_dir: Path,
    demo: bool = False,
    genetic_code: int = 1,
    allow_non_canonical: bool = False,
    no_plot: bool = False,
) -> dict:
    """Validate gene models from a GTF annotation against a reference genome FASTA.

    Checks each transcript for start codon, stop codon, internal stop codons,
    canonical splice sites, minimum intron/exon lengths, and overlapping exons,
    reporting errors and the translated protein sequence.

    Args:
        annotation: GTF gene annotation to validate. None if demo.
        fasta: Reference genome FASTA. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        genetic_code: NCBI genetic code table (default 1 = standard).
        allow_non_canonical: Allow GC/AT donor splice sites without flagging.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'validations', 'n_valid', 'n_invalid', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Gene Model Validator v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_validator.py"
    spec = importlib.util.spec_from_file_location("run_validator", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_validator.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
