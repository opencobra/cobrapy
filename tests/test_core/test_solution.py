"""Test functions of solution.py ."""

from typing import TYPE_CHECKING, Tuple

from cobra.core import Solution


if TYPE_CHECKING:
    from cobra import Model


def test_solution_contains_only_reaction_specific_values(
    solved_model: Tuple[Solution, "Model"]
) -> None:
    """Test solution contains specific reaction values."""
    solution, model = solved_model
    reaction_ids = set([reaction.id for reaction in model.reactions])
    assert set(solution.fluxes.index) == reaction_ids
    # assert set(solution.reduced_costs.index) == reaction_ids
