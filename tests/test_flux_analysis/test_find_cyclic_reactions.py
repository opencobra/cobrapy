"""Test finding cyclic reactions in metabolic model."""

from typing import Callable, List

import pytest

from cobra import Metabolite, Model, Reaction
from cobra.flux_analysis.find_cyclic_reactions import find_cyclic_reactions


@pytest.mark.parametrize("method", ["basic", "optimized"])
def test_find_cyclic_reactions_benchmark(
    large_model: Model, benchmark: Callable, all_solvers: List[str], method: str
) -> None:
    """Benchmark find_cyclic_reactions."""
    large_model.solver = all_solvers
    benchmark(
        find_cyclic_reactions,
        large_model,
        method=method,
    )


@pytest.mark.parametrize("method", ["basic", "optimized"])
def test_find_cyclic_reactions(
    model: Model, all_solvers: List[str], method: str
) -> None:
    """Test find_cyclic_reactions."""
    model.solver = all_solvers
    cyclic_reactions, cyclic_directions = find_cyclic_reactions(model, method=method)

    assert len(cyclic_reactions) == 2
    assert len(cyclic_directions) == 2
    for reaction in ["FRD7", "SUCDi"]:
        assert reaction in cyclic_reactions, f"Expected {reaction} in cyclic reactions"
    for directions in cyclic_directions:
        assert directions == (False, True)


@pytest.mark.parametrize("method", ["basic", "optimized"])
@pytest.mark.parametrize("boundary_only", [False, True])
def test_find_cyclic_reactions_no_internal(method: str, boundary_only: bool) -> None:
    """No internal reactions means no internal cyclic directions."""
    model = Model()
    model.solver = "glpk"
    if boundary_only:
        metabolite = Metabolite("A")
        uptake = Reaction("uptake")
        uptake.add_metabolites({metabolite: 1})
        drain = Reaction("drain")
        drain.add_metabolites({metabolite: -1})
        model.add_reactions([uptake, drain])

    assert find_cyclic_reactions(model, method=method) == ([], [])
