"""Shared pytest fixtures for the validation harness."""

import pytest


@pytest.fixture(scope="session")
def simulated_damage_reads(tmp_path_factory):
    """Generate a small set of simulated ancient reads via post-mortem-damage-simulator.

    TODO: wire up actual simulator call once src/pipeline.py is implemented.
    Returns path to output directory containing simulated_reads.fastq.gz.
    """
    raise NotImplementedError("post-mortem-damage-simulator not yet implemented")


@pytest.fixture(scope="session")
def msprime_tree_sequence():
    """Return an msprime TreeSequence simulating a proxy-species bottleneck scenario.

    TODO: implement scenario parameterised by mammoth/elephant demographic estimates.
    """
    raise NotImplementedError("msprime fixture not yet implemented")


@pytest.fixture(scope="session")
def slim_vcf(tmp_path_factory):
    """Return path to a SLiM-simulated VCF with planted selection sweeps.

    TODO: implement SLiM script for proxy-species positive-selection scenario.
    """
    raise NotImplementedError("SLiM fixture not yet implemented")
