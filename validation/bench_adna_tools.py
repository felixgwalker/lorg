"""Benchmark aDNA damage and authentication PIVOT tools against simulated ground truth.

Simulation source: validation/post-mortem-damage-simulator (Briggs model).
Each test uses the `simulated_damage_reads` session fixture from conftest.py.

TODO: implement each benchmark once the corresponding PIVOT tool is functional.
"""

import pytest


@pytest.mark.benchmark
def test_damage_classifier_ct_recovery(simulated_damage_reads):
    """ancient-dna-damage-classifier: recovered C->T rates match simulated rates."""
    # TODO:
    # 1. Run ancient-dna-damage-classifier on simulated_damage_reads.
    # 2. Load output damage frequency table.
    # 3. Assert abs(observed_ct_rate_pos1 - simulated_ct_rate) < tolerance.
    raise NotImplementedError


@pytest.mark.benchmark
def test_contamination_estimator_fpr(simulated_damage_reads):
    """contamination-estimator: false-positive contamination rate < 5 % on clean simulations."""
    # TODO: run contamination-estimator on undamaged synthetic reads; assert FPR < 0.05.
    raise NotImplementedError


@pytest.mark.benchmark
def test_coverage_assessor_endogenous_fraction(simulated_damage_reads):
    """palaeogenomic-coverage-assessor: endogenous fraction estimate within 2 % of truth."""
    # TODO: run palaeogenomic-coverage-assessor on simulated reads with known endogenous %.
    raise NotImplementedError


@pytest.mark.benchmark
def test_damage_classifier_authenticity_auc(simulated_damage_reads):
    """ancient-dna-damage-classifier: AUROC > 0.95 on ancient vs modern-contaminated reads."""
    # TODO: mix authentic and contaminated reads; check classification AUROC.
    raise NotImplementedError
