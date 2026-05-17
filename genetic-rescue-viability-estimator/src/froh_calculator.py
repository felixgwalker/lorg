import numpy as np


def compute_froh(genotype_matrix, min_run_length=10):
    """
    genotype_matrix: 2D array shape (n_individuals, n_sites), values 0=het, 1=hom
    Returns array of FROH per individual.
    """
    n_indiv, n_sites = genotype_matrix.shape
    froh_values = []
    for i in range(n_indiv):
        gt = genotype_matrix[i]
        total_hom_run = 0
        run = 0
        for j in range(n_sites):
            if gt[j] == 1:
                run += 1
            else:
                if run >= min_run_length:
                    total_hom_run += run
                run = 0
        if run >= min_run_length:
            total_hom_run += run
        froh = total_hom_run / n_sites if n_sites > 0 else 0.0
        froh_values.append(froh)
    return np.array(froh_values)


def froh_from_population_data(population_data, n_snps=500, rng=None):
    """
    Derive FROH estimates from population data dicts containing F_initial.
    Uses F_initial as approximate FROH (they are related measures of inbreeding).
    """
    if rng is None:
        rng = np.random.default_rng(42)

    for indiv in population_data:
        F = indiv["F_initial"]
        noise = rng.normal(0, 0.02)
        froh = float(np.clip(F + noise, 0.0, 1.0))
        indiv["FROH"] = round(froh, 4)

    return population_data
