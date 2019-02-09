# -*- coding: utf-8 -*-

"""Test functionalities of FASTCC."""


from __future__ import absolute_import

import pytest

import cobra
from cobra.flux_analysis import fastcc


@pytest.fixture(scope="module")
def figure1_model():
    """
    Generate a toy model as described in [1]_ figure 1.

    References
    ----------
    .. [1] Vlassis N, Pacheco MP, Sauter T (2014)
           Fast Reconstruction of Compact Context-Specific Metabolic Network
           Models.
           PLoS Comput Biol 10(1): e1003424. doi:10.1371/journal.pcbi.1003424

    """
    test_model = cobra.Model("figure 1")
    v1 = cobra.Reaction("v1")
    v2 = cobra.Reaction("v2")
    v3 = cobra.Reaction("v3")
    v4 = cobra.Reaction("v4")
    v5 = cobra.Reaction("v5")
    v6 = cobra.Reaction("v6")

    test_model.add_reactions([v1, v2, v3, v4, v5, v6])

    v1.reaction = "-> 2 A"
    v2.reaction = "A <-> B"
    v3.reaction = "A -> D"
    v4.reaction = "A -> C"
    v5.reaction = "C -> D"
    v6.reaction = "D ->"

    v1.bounds = (0.0, 3.0)
    v2.bounds = (-3.0, 3.0)
    v3.bounds = (0.0, 3.0)
    v4.bounds = (0.0, 3.0)
    v5.bounds = (0.0, 3.0)
    v6.bounds = (0.0, 3.0)

    test_model.objective = v6
    return test_model


@pytest.fixture(scope="module")
def opposing_model():
    """
    Generate a toy model with opposing reversible reactions.

    This toy model ensures that two opposing reversible reactions do not
    appear as blocked.

    """
    test_model = cobra.Model("opposing")
    v1 = cobra.Reaction("v1")
    v2 = cobra.Reaction("v2")
    v3 = cobra.Reaction("v3")
    v4 = cobra.Reaction("v4")

    test_model.add_reactions([v1, v2, v3, v4])

    v1.reaction = "-> 2 A"
    v2.reaction = "A -> C"  # Later made reversible via bounds.
    v3.reaction = "D -> C"  # Later made reversible via bounds.
    v4.reaction = "D ->"

    v1.bounds = 0.0, 3.0
    v2.bounds = -3.0, 3.0
    v3.bounds = -3.0, 3.0
    v4.bounds = 0.0, 3.0

    test_model.objective = v4
    return test_model


def test_fastcc_benchmark(model, benchmark, all_solvers):
    """Benchmark fastcc."""
    model.solver = all_solvers
    benchmark(fastcc, model)


def test_figure1(figure1_model, all_solvers):
    """Test fastcc."""
    figure1_model.solver = all_solvers
    consistent_model = fastcc(figure1_model)
    expected_reactions = {'v1', 'v3', 'v4', 'v5', 'v6'}
    assert expected_reactions == {rxn.id for rxn in consistent_model.reactions}


def test_opposing(opposing_model, all_solvers):
    """Test fastcc."""
    opposing_model.solver = all_solvers
    consistent_model = fastcc(opposing_model)
    expected_reactions = {'v1', 'v2', 'v3', 'v4'}
    assert expected_reactions == {rxn.id for rxn in consistent_model.reactions}
