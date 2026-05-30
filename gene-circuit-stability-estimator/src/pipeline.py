"""Gene Circuit Stability Estimator — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    circuit_json: Path | None,
    output_dir: Path,
    demo: bool = False,
    simulation_time: float = 100.0,
    n_perturbations: int = 100,
    no_plot: bool = False,
) -> dict:
    """Estimate the stability and dynamic behaviour of a synthetic gene circuit.

    Simulates circuit ODE dynamics from multiple initial conditions, identifies
    steady states, classifies circuit behaviour (stable/oscillating/bistable/
    unstable), and computes a robustness score under parameter perturbation.

    Args:
        circuit_json: JSON of circuit nodes (proteins, rates) and edges (interactions). None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        simulation_time: Simulation duration in hours.
        n_perturbations: Number of parameter perturbations for robustness analysis.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'analysis', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Gene Circuit Stability Estimator v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_estimator.py"
    spec = importlib.util.spec_from_file_location("run_estimator", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_estimator.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
