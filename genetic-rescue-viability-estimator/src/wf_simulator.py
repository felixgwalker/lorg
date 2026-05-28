import numpy as np


def simulate_wright_fisher(
    Ne: int,
    n_gen: int,
    initial_freq: float,
    n_reps: int = 100,
    seed: int = 42,
) -> np.ndarray:
    """Wright-Fisher allele frequency drift simulation.

    Each generation the allele count is drawn from Binomial(2*Ne, p) where p
    is the current frequency, modelling pure genetic drift.

    Parameters
    ----------
    Ne          : effective population size.
    n_gen       : number of generations to simulate.
    initial_freq: starting allele frequency (in [0, 1]).
    n_reps      : number of independent replicate trajectories.
    seed        : random seed for reproducibility.

    Returns
    -------
    np.ndarray of shape (n_reps, n_gen) containing the allele frequency at
    each generation (generation 0 / initial state is NOT included so the
    array has exactly n_gen columns, one per generation step).
    """
    rng = np.random.default_rng(seed)
    Ne = max(1, int(Ne))
    n_alleles = 2 * Ne

    # Initialise frequencies: shape (n_reps,)
    freqs = np.full(n_reps, float(initial_freq))

    trajectories = np.zeros((n_reps, n_gen), dtype=float)

    for gen in range(n_gen):
        counts = rng.binomial(n=n_alleles, p=freqs)
        freqs = counts / n_alleles
        trajectories[:, gen] = freqs

    return trajectories


# ---------------------------------------------------------------------------
# Legacy simulation functions used by the pipeline
# ---------------------------------------------------------------------------

def run_wf_simulation(
    initial_F,
    N_current,
    N_target,
    n_rescue,
    n_generations,
    n_replicates,
    selection_coeff,
    F_rescue=0.01,
    rng=None,
):
    if rng is None:
        rng = np.random.default_rng(42)

    H_trajectories = np.zeros((n_replicates, n_generations + 1))
    F_trajectories = np.zeros((n_replicates, n_generations + 1))

    for rep in range(n_replicates):
        F = float(initial_F)
        H_trajectories[rep, 0] = 1 - F
        F_trajectories[rep, 0] = F

        rescue_applied = False
        rescue_gen = 1

        for gen in range(1, n_generations + 1):
            N_eff = N_current + (N_target - N_current) * (gen / n_generations)
            N_eff = max(int(N_eff), 1)

            if gen == rescue_gen and not rescue_applied:
                F = (N_current * F + n_rescue * F_rescue) / (N_current + n_rescue)
                rescue_applied = True
                N_eff = N_current + n_rescue

            F = 1 / (2 * N_eff) + (1 - 1 / (2 * N_eff)) * F
            F = float(np.clip(F, 0.0, 1.0))

            H_trajectories[rep, gen] = 1 - F
            F_trajectories[rep, gen] = F

    mean_H = H_trajectories.mean(axis=0)
    ci_lower = np.percentile(H_trajectories, 2.5, axis=0)
    ci_upper = np.percentile(H_trajectories, 97.5, axis=0)
    mean_F = F_trajectories.mean(axis=0)

    return {
        "H_trajectories": H_trajectories,
        "F_trajectories": F_trajectories,
        "mean_H": mean_H,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "mean_F": mean_F,
        "generations": np.arange(n_generations + 1),
    }


def run_baseline_simulation(initial_F, N_current, N_target, n_generations, n_replicates, rng=None):
    if rng is None:
        rng = np.random.default_rng(99)

    H_trajectories = np.zeros((n_replicates, n_generations + 1))
    F_trajectories = np.zeros((n_replicates, n_generations + 1))

    for rep in range(n_replicates):
        F = float(initial_F)
        H_trajectories[rep, 0] = 1 - F
        F_trajectories[rep, 0] = F

        for gen in range(1, n_generations + 1):
            N_eff = N_current + (N_target - N_current) * (gen / n_generations)
            N_eff = max(int(N_eff), 1)
            F = 1 / (2 * N_eff) + (1 - 1 / (2 * N_eff)) * F
            F = float(np.clip(F, 0.0, 1.0))
            H_trajectories[rep, gen] = 1 - F
            F_trajectories[rep, gen] = F

    return {
        "H_trajectories": H_trajectories,
        "F_trajectories": F_trajectories,
        "mean_H": H_trajectories.mean(axis=0),
        "mean_F": F_trajectories.mean(axis=0),
        "ci_lower": np.percentile(H_trajectories, 2.5, axis=0),
        "ci_upper": np.percentile(H_trajectories, 97.5, axis=0),
        "generations": np.arange(n_generations + 1),
    }
