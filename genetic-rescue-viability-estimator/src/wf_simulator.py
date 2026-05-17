import numpy as np


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
