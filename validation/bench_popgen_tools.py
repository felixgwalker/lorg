"""Benchmark population-genetics PIVOT tools against msprime / SLiM simulations.

Each test uses session fixtures from conftest.py.

TODO: implement each benchmark once the corresponding PIVOT tool is functional.
"""

import pytest


@pytest.mark.benchmark
def test_ne_estimator_rmse(msprime_tree_sequence):
    """effective-population-size-estimator: Ne RMSE < 20 % relative error on temporal samples."""
    # TODO:
    # 1. Convert msprime_tree_sequence to VCF with ancient tip-dating samples.
    # 2. Run effective-population-size-estimator in temporal/aDNA mode.
    # 3. Compare estimated Ne trajectory to msprime truth.
    raise NotImplementedError


@pytest.mark.benchmark
def test_roh_interpreter_precision_recall(msprime_tree_sequence):
    """roh-interpreter: ROH detection recall > 0.85 at low (2x) simulated coverage."""
    # TODO: downsample msprime VCF to 2x, run roh-interpreter, compare to true ROH segments.
    raise NotImplementedError


@pytest.mark.benchmark
def test_introgression_detector_dstat_power(msprime_tree_sequence):
    """introgression-detector: D-stat detects planted ancient introgression at 5 % proportion."""
    # TODO: simulate admixture event, verify D-stat p-value < 0.01.
    raise NotImplementedError


@pytest.mark.benchmark
def test_selection_detector_sweep_recall(slim_vcf):
    """positive-selection-signal-detector: planted sweeps recovered at recall > 0.80."""
    # TODO: run positive-selection-signal-detector on SLiM VCF with known sweep positions.
    raise NotImplementedError
