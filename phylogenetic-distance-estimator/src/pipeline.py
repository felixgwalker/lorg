"""Phylogenetic Distance Estimator — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

try:
    from src import __version__
    from src.alignment_reader import generate_demo_alignment, read_alignment, validate_alignment
    from src.distance_calculator import compute_distance_matrix, normalize_matrix
    from src.nj_tree import build_nj_tree
    from src.plot import plot_heatmap_dendrogram
    from src.report import write_distance_matrix, write_newick, write_ranked_comparisons
except ImportError:
    from __init__ import __version__  # type: ignore[no-redef]
    from alignment_reader import generate_demo_alignment, read_alignment, validate_alignment  # type: ignore[no-redef]
    from distance_calculator import compute_distance_matrix, normalize_matrix  # type: ignore[no-redef]
    from nj_tree import build_nj_tree  # type: ignore[no-redef]
    from plot import plot_heatmap_dendrogram  # type: ignore[no-redef]
    from report import write_distance_matrix, write_newick, write_ranked_comparisons  # type: ignore[no-redef]

logger = logging.getLogger(__name__)


def run_pipeline(
    fasta_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    distance_model: str = "K2P",
    no_plot: bool = False,
) -> dict:
    """Execute the full Phylogenetic Distance Estimator pipeline."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=== Phylogenetic Distance Estimator v%s ===", __version__)

    if demo:
        logger.info("Demo mode: generating synthetic alignment.")
        sequences = generate_demo_alignment(n_species=5, seq_length=500)
    else:
        if fasta_path is None:
            raise ValueError("FASTA path required when not in demo mode.")
        logger.info("Reading alignment: %s", fasta_path)
        sequences = read_alignment(str(fasta_path))

    validate_alignment(sequences)
    logger.info("Alignment: %d sequences, length %d bp.", len(sequences), len(sequences[0][1]))

    logger.info("Computing pairwise distances using %s model.", distance_model)
    names, matrix = compute_distance_matrix(sequences, model=distance_model)
    normalized = normalize_matrix(matrix)

    logger.info("Building Neighbor-Joining tree.")
    newick = build_nj_tree(names, matrix)

    dist_csv = write_distance_matrix(names, matrix, output_dir)
    nwk_path = write_newick(newick, output_dir)
    ranked_csv = write_ranked_comparisons(names, matrix, normalized, output_dir)
    logger.info("Written: %s", dist_csv)
    logger.info("Written: %s", nwk_path)
    logger.info("Written: %s", ranked_csv)

    output_files: dict[str, str] = {
        "distance_matrix_csv": str(dist_csv),
        "newick_tree": str(nwk_path),
        "ranked_comparisons_csv": str(ranked_csv),
    }

    if no_plot:
        logger.info("Plot generation skipped (--no-plot).")
    else:
        try:
            png, svg = plot_heatmap_dendrogram(names, matrix, output_dir)
            output_files["plot_png"] = str(png)
            output_files["plot_svg"] = str(svg)
            logger.info("Plot written: %s", png)
        except Exception as exc:
            logger.warning("Plot generation failed: %s", exc)

    logger.info("Done. Output: %s", output_dir)
    return {
        "sequences": sequences,
        "names": names,
        "distance_matrix": matrix,
        "normalized_matrix": normalized,
        "newick": newick,
        "output_files": output_files,
        "pipeline_version": __version__,
    }


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_estimator.py"
    spec = importlib.util.spec_from_file_location("run_estimator", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_estimator.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
