"""Prime Edit Design Assistant — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    fasta_path: Path | None,
    edit_json: Path | None,
    output_dir: Path,
    demo: bool = False,
    pam: str = "NGG",
    pbs_min: int = 8,
    pbs_max: int = 15,
    rt_min: int = 10,
    rt_max: int = 16,
    no_plot: bool = False,
) -> dict:
    """Design pegRNAs for prime editing.

    Scans the target locus FASTA for PAM sites, enumerates RT template and PBS
    length combinations, and ranks designs by PBS GC content, RT length, spacer
    GC, and homopolymer penalties.  Optionally designs PE3 nicking guides.

    Args:
        fasta_path: Target locus FASTA (≥200 bp around edit site). None if demo.
        edit_json: JSON specifying {"position": int, "ref": str, "alt": str}.
        output_dir: Directory for output TSV and optional plot.
        demo: Generate synthetic data and run without real inputs.
        pam: PAM motif (default NGG for SpCas9 PE2/PE3).
        pbs_min/pbs_max: Primer binding site length search range (nt).
        rt_min/rt_max: Reverse transcriptase template length search range (nt).
        no_plot: Skip matplotlib figure generation.

    Returns:
        dict with keys 'designs' (list), 'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Prime Edit Design Assistant v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_assistant.py"
    spec = importlib.util.spec_from_file_location("run_assistant", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_assistant.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
