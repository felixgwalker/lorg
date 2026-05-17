import os
import numpy as np

from .vcf_parser import parse_vcf, generate_synthetic_population
from .froh_calculator import froh_from_population_data
from .wf_simulator import run_wf_simulation, run_baseline_simulation
from .viability_scorer import compute_viability_score
from .report import write_froh_csv, write_trajectories_csv, write_viability_json
from .plot import plot_trajectories, plot_distribution


def run_pipeline(args):
    os.makedirs(args.output_dir, exist_ok=True)
    rng = np.random.default_rng(args.seed if hasattr(args, "seed") and args.seed else 42)

    if args.demo:
        population_data = generate_synthetic_population(
            n_individuals=20, n_snps=500, rng=rng
        )
        N_current = 50
        N_target = 200
        n_rescue = 5
        n_generations = 50
        n_replicates = 50
        selection_coeff = 0.1
    else:
        if args.vcf:
            population_data = parse_vcf(args.vcf)
        else:
            raise ValueError("Provide --vcf or use --demo")
        N_current = args.N_current
        N_target = args.N_target
        n_rescue = args.n_rescue
        n_generations = args.n_generations
        n_replicates = args.n_replicates
        selection_coeff = args.selection_coeff

    population_data = froh_from_population_data(population_data, rng=rng)

    initial_F = float(np.mean([d["F_initial"] for d in population_data]))

    sim_result = run_wf_simulation(
        initial_F=initial_F,
        N_current=N_current,
        N_target=N_target,
        n_rescue=n_rescue,
        n_generations=n_generations,
        n_replicates=n_replicates,
        selection_coeff=selection_coeff,
        rng=rng,
    )

    baseline_result = run_baseline_simulation(
        initial_F=initial_F,
        N_current=N_current,
        N_target=N_target,
        n_generations=n_generations,
        n_replicates=n_replicates,
        rng=np.random.default_rng(99),
    )

    score_result = compute_viability_score(sim_result, initial_F)

    froh_path = write_froh_csv(population_data, args.output_dir)
    traj_path = write_trajectories_csv(sim_result, args.output_dir)
    score_path = write_viability_json(score_result, args.output_dir)

    plot_paths = []
    if not getattr(args, "no_plot", False):
        fmt = getattr(args, "plot_format", "png")
        p1 = plot_trajectories(sim_result, baseline_result, args.output_dir, fmt=fmt)
        p2 = plot_distribution(sim_result, baseline_result, args.output_dir, fmt=fmt)
        plot_paths = [p1, p2]

    return {
        "population_data": population_data,
        "sim_result": sim_result,
        "score_result": score_result,
        "outputs": {
            "froh_csv": froh_path,
            "trajectories_csv": traj_path,
            "viability_json": score_path,
            "plots": plot_paths,
        },
    }


def main():
    import sys
    import argparse

    parser = argparse.ArgumentParser(
        description="Genetic Rescue Viability Estimator"
    )
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--vcf", type=str, default=None)
    parser.add_argument("--output-dir", dest="output_dir", type=str, default="output_estimator")
    parser.add_argument("--no-plot", dest="no_plot", action="store_true")
    parser.add_argument("--plot-format", dest="plot_format", choices=["png", "svg"], default="png")
    parser.add_argument("--N-current", dest="N_current", type=int, default=50)
    parser.add_argument("--N-target", dest="N_target", type=int, default=200)
    parser.add_argument("--n-rescue", dest="n_rescue", type=int, default=5)
    parser.add_argument("--n-generations", dest="n_generations", type=int, default=50)
    parser.add_argument("--n-replicates", dest="n_replicates", type=int, default=100)
    parser.add_argument("--selection-coeff", dest="selection_coeff", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    if not args.demo and args.vcf is None:
        parser.error("Provide --vcf or use --demo.")

    result = run_pipeline(args)
    sr = result["score_result"]
    print(f"Viability Score: {sr['score']:.1f}/100  Grade: {sr['grade']}")
    print(f"Interpretation: {sr['interpretation']}")
    print(f"Outputs written to: {args.output_dir}")
    return 0
