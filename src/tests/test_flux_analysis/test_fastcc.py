"""Test functionalities of FASTCC."""

from typing import Callable, List

import pytest

from cobra import Model, Reaction
from cobra.flux_analysis import fastcc, flux_variability_analysis


@pytest.fixture(scope="module")
def figure1_model() -> Model:
    """Generate a toy model as described in [1]_ figure 1.

    References
    ----------
    .. [1] Vlassis N, Pacheco MP, Sauter T (2014)
           Fast Reconstruction of Compact Context-Specific Metabolic Network
           Models.
           PLoS Comput Biol 10(1): e1003424. doi:10.1371/journal.pcbi.1003424

    """
    test_model = Model("figure 1")
    v1 = Reaction("v1")
    v2 = Reaction("v2")
    v3 = Reaction("v3")
    v4 = Reaction("v4")
    v5 = Reaction("v5")
    v6 = Reaction("v6")

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
def opposing_model() -> Model:
    """Generate a toy model with opposing reversible reactions.

    This toy model ensures that two opposing reversible reactions do not
    appear as blocked.

    """
    test_model = Model("opposing")
    v1 = Reaction("v1")
    v2 = Reaction("v2")
    v3 = Reaction("v3")
    v4 = Reaction("v4")

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


def test_fastcc_benchmark(
    model: Model, benchmark: Callable, all_solvers: List[str]
) -> None:
    """Benchmark fastcc."""
    model.solver = all_solvers
    benchmark(fastcc, model)


def test_figure1(figure1_model: Model, all_solvers: List[str]) -> None:
    """Test FASTCC."""
    figure1_model.solver = all_solvers
    consistent_model = fastcc(figure1_model)
    expected_reactions = {"v1", "v3", "v4", "v5", "v6"}
    assert expected_reactions == {rxn.id for rxn in consistent_model.reactions}


def test_opposing(opposing_model: Model, all_solvers: List[str]) -> None:
    """Test FASTCC."""
    opposing_model.solver = all_solvers
    consistent_model = fastcc(opposing_model)
    expected_reactions = {"v1", "v2", "v3", "v4"}
    assert expected_reactions == {rxn.id for rxn in consistent_model.reactions}


def test_fastcc_against_fva_nonblocked_rxns(
    model: Model, all_solvers: List[str]
) -> None:
    """Test non-blocked reactions obtained by FASTCC against FVA."""
    model.solver = all_solvers
    fastcc_consistent_model = fastcc(model)
    fva = flux_variability_analysis(model, fraction_of_optimum=0.0)
    nonblocked_rxns_fva = fva.index[
        (fva.minimum.abs() > model.tolerance) | (fva.maximum.abs() > model.tolerance)
    ]
    assert all(fastcc_consistent_model.reactions) == all(nonblocked_rxns_fva.tolist())
