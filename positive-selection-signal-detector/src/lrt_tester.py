import numpy as np
from scipy import stats


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
