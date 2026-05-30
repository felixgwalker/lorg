"""UTR Variant Analyser — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    vcf_path: Path | None,
    fasta_path: Path | None,
    annotation_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    utr_type: str = "both",
    kozak_threshold: float = 0.5,
    no_plot: bool = False,
) -> dict:
    """Analyse variants in 5'/3' UTR regions.

    Identifies uORF creation and disruption events, changes in Kozak context
    strength, and disruption of polyadenylation signals in 3' UTRs.

    Args:
        vcf_path: VCF of UTR variants. None if demo.
        fasta_path: Reference genome FASTA. None if demo.
        annotation_path: UTR annotation BED or GTF. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        utr_type: Which UTRs to analyse (5prime, 3prime, or both).
        kozak_threshold: Kozak score delta to flag as significant.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'results', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== UTR Variant Analyser v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_analyser.py"
    spec = importlib.util.spec_from_file_location("run_analyser", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_analyser.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
