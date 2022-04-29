"""Contains module level fixtures and utility functions."""

from pathlib import Path
from pickle import load
from typing import Callable

import pytest
from optlang.interface import OPTIMAL

from cobra import Model


@pytest.fixture(scope="module")
def mini_model(data_directory: Path) -> Model:
    """Fixture for mini model."""
    with (data_directory / "mini.pickle").open(mode="rb") as infile:
        return load(infile)


@pytest.fixture
def compare_models() -> Callable:
    """Compare models as a fixture.

    Returns
    -------
    _compare_models: Callable
        A function that will assert that the models are identical.

    """

    def _compare_models(model_1: Model, model_2: Model) -> None:
        """Compare two models (only for testing purposes).

        Parameters
        ----------
        model_1: Model
            First model to compare.

        model_2: Model
            Second model to compare.

        This function does not return anything, and does a series of asserts on the
        two models. If the models are not identical, the asserts, and the test that
        called them will fail.

        """
        assert len(model_1.reactions) == len(model_2.reactions)
        assert len(model_1.metabolites) == len(model_2.metabolites)
        assert len(model_1.genes) == len(model_2.genes)
        assert model_1.objective.direction == model_2.objective.direction

        # check Reaction attributes
        for attr in (
            "id",
            "name",
            "lower_bound",
            "upper_bound",
            "objective_coefficient",
            "gene_reaction_rule",
        ):
            assert getattr(model_1.reactions[0], attr) == getattr(
                model_2.reactions[0], attr
            )
            assert getattr(model_1.reactions[5], attr) == getattr(
                model_2.reactions[5], attr
            )
            assert getattr(model_1.reactions[-1], attr) == getattr(
                model_2.reactions[-1], attr
            )

        # check Metabolite attributes
        for attr in ("id", "name", "compartment", "formula", "charge"):
            assert getattr(model_1.metabolites[0], attr) == getattr(
                model_2.metabolites[0], attr
            )
            assert getattr(model_1.metabolites[5], attr) == getattr(
                model_2.metabolites[5], attr
            )
            assert getattr(model_1.metabolites[-1], attr) == getattr(
                model_2.metabolites[-1], attr
            )
            assert len(model_1.reactions[0].metabolites) == len(
                model_2.reactions[0].metabolites
            )
        # TODO: either relax gene attribute checking or fix models for testing.
        # check Gene attributes
        # for attr in ("id", "name"):
        #     assert getattr(model_1.genes[0], attr) == getattr(model_2.genes[0],
        #                                                       attr)
        #     assert getattr(model_1.genes[10], attr) == getattr(model_2.genes[10],
        #                                                        attr)
        #     assert getattr(model_1.genes[-1], attr) == getattr(model_2.genes[-1],
        #                                                        attr)

        assert len(model_1.reactions[8].metabolites) == len(
            model_2.reactions[8].metabolites
        )
        assert len(model_1.reactions[-1].metabolites) == len(
            model_2.reactions[-1].metabolites
        )
        assert len(model_1.genes) == len(model_2.genes)

        # ensure they have the same solution max
        solution_1 = model_1.optimize()
        solution_2 = model_2.optimize()
        if solution_1.status == OPTIMAL and solution_2.status == OPTIMAL:
            assert abs(
                solution_1.objective_value - solution_2.objective_value
            ) == pytest.approx(0.0)
        else:
            assert solution_1.status == solution_2.status

        # ensure the references are correct
        # metabolite -> model reference
        assert model_1.metabolites[0]._model is model_1
        assert model_2.metabolites[0]._model is model_2

        # reaction -> model reference
        assert model_1.reactions[0]._model is model_1
        assert model_2.reactions[0]._model is model_2

        # gene -> model reference
        if model_1.genes:
            assert model_1.genes[0]._model is model_1
            assert model_2.genes[0]._model is model_2

        # extra comparisons
        # assert model_1.compartments == model_2.compartments
        # assert dict(model_1.metabolites[4].annotation) == dict(
        #     model_2.metabolites[4].annotation)
        # assert dict(model_1.reactions[4].annotation) == dict(
        #     model_2.reactions[4].annotation)
        # assert dict(model_1.genes[5].annotation) == dict(
        #     model_2.genes[5].annotation)

    return _compare_models
