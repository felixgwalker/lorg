import numpy as np
from scipy import stats


def lrt_test(dnds_results):
    """Apply a likelihood ratio test to a list of per-window dN/dS results.

    For each window result apply an LRT comparing:
    - Null model  (M0): dN = dS, omega = 1  — one fewer free parameter
    - Alternative (M1): omega unconstrained

    The test statistic is 2 * delta_lnL ~ chi-squared(1 df).

    Parameters
    ----------
    dnds_results : list of dict
        Each dict must contain keys 'window_pos' and 'omega'.  Additional keys
        ('dN', 'dS', 'n_pairs') are used when present to improve the
        log-likelihood approximation.

    Returns
    -------
    list of (window_pos, omega, lrt_stat, p_raw) tuples
    """
    output = []
    for entry in dnds_results:
        window_pos = entry.get("window_pos", 0)
        omega = entry.get("omega", 1.0)
        n_pairs = entry.get("n_pairs", 1)

        omega_clipped = min(float(omega), 10.0) if not np.isinf(omega) else 10.0

        if omega_clipped <= 1.0:
            lrt_stat = 0.0
            p_raw = 1.0
        else:
            # Approximate log-likelihood based on observed omega and n_pairs.
            # lnL0 (null: omega=1):  each pair contributes -0.5 (unit variance at omega=1)
            # lnL1 (alt):  gain is proportional to log(omega) per pair
            lnL0 = -n_pairs * 0.5
            lnL1 = lnL0 + n_pairs * np.log(omega_clipped) * 0.5
            lrt_stat = float(max(2.0 * (lnL1 - lnL0), 0.0))
            p_raw = float(stats.chi2.sf(lrt_stat, df=1))

        output.append((window_pos, omega, lrt_stat, p_raw))
    return output


def compute_lrt_pvalue(omega, n_pairs):
    if n_pairs == 0:
        return 1.0

    omega_clipped = min(omega, 10.0)

    if omega_clipped <= 1.0:
        return 1.0

    excess = omega_clipped - 1.0
    lnL0 = -n_pairs * 0.5
    lnL1 = lnL0 + n_pairs * np.log(max(omega_clipped, 1.0)) * 0.5

    lrt_stat = 2 * (lnL1 - lnL0)
    lrt_stat = max(lrt_stat, 0.0)

    pval = stats.chi2.sf(lrt_stat, df=1)
    return float(pval)


def test_all_genes(dnds_results):
    tested = []
    for gene, result in dnds_results.items():
        omega = result["omega"]
        n_pairs = result["n_pairs"]
        pval = compute_lrt_pvalue(omega, n_pairs)
        entry = dict(result)
        entry["gene"] = gene
        entry["lrt_pval"] = pval
        entry["selection_signal"] = omega > 1.5
        tested.append(entry)
    return tested
