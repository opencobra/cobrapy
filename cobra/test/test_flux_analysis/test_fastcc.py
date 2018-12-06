# -*- coding: utf-8 -*-

"""Test functionalities of FASTCC."""

from __future__ import absolute_import

import pytest

from cobra.core import Model, Reaction
from cobra.flux_analysis import fastcc


@pytest.fixture(scope="module")
def fastcc_model():
    """
    Generate FASTCC model as described in [1]_

    References
    ----------
    .. [1] Vlassis N, Pacheco MP, Sauter T (2014)
           Fast Reconstruction of Compact Context-Specific Metabolic Network
           Models.
           PLoS Comput Biol 10(1): e1003424. doi:10.1371/journal.pcbi.1003424
    """
    test_model = Model("FASTCC_model")
    r1 = Reaction("r1")
    r2 = Reaction("r2")
    r3 = Reaction("r3")
    r4 = Reaction("r4")
    r5 = Reaction("r5")
    r6 = Reaction("r6")

    test_model.add_reactions([r1, r2, r3, r4, r5, r6])

    r1.reaction = "-> 2 A"
    r2.reaction = "A <-> B"
    r3.reaction = "A -> D"
    r4.reaction = "A -> C"
    r5.reaction = "C -> D"
    r6.reaction = "D ->"

    r1.bounds = (0.0, 3.0)
    r2.bounds = (-3.0, 3.0)
    r3.bounds = (0.0, 3.0)
    r4.bounds = (0.0, 3.0)
    r5.bounds = (0.0, 3.0)
    r6.bounds = (0.0, 3.0)

    test_model.objective = r6
    return test_model


def test_fastcc_benchmark(model, benchmark, all_solvers):
    """Benchmark fastcc."""
    model.solver = all_solvers
    benchmark(fastcc, model)


def test_fastcc(fastcc_model, all_solvers):
    """Test fastcc."""
    fastcc_model.solver = all_solvers
    consistent_model = fastcc(fastcc_model)
    expected_reactions = ['r1', 'r3', 'r4', 'r5', 'r6']
    assert expected_reactions == [rxn.id for rxn in consistent_model.reactions]
