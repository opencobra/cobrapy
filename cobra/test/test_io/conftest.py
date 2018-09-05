# -*- coding: utf-8 -*-

"""Contains module level fixtures and utility functions."""

from __future__ import absolute_import

from os.path import join
from pickle import dump, load

import pytest


@pytest.fixture(scope="module")
def mini_model(data_directory):
    """Fixture for mini model."""
    with open(join(data_directory, "mini.pickle"), "rb") as infile:
        return load(infile)


def compare_models(model_1, model_2):
    """Compare two models (only for testing purposes)."""
    assert len(model_1.reactions) == len(model_2.reactions)
    assert len(model_1.metabolites) == len(model_2.metabolites)
    assert len(model_1.genes) == len(model_2.genes)
    assert model_1.objective.direction == model_2.objective.direction

    # check Reaction attributes
    for attr in ("id", "name", "lower_bound", "upper_bound",
                 "objective_coefficient", "gene_reaction_rule"):
        assert getattr(model_1.reactions[0], attr) == getattr(
            model_2.reactions[0], attr)
        assert getattr(model_1.reactions[5], attr) == getattr(
            model_2.reactions[5], attr)
        assert getattr(model_1.reactions[-1], attr) == getattr(
            model_2.reactions[-1], attr)

    # check Metabolite attributes
    for attr in ("id", "name", "compartment", "formula", "charge"):
        assert getattr(model_1.metabolites[0], attr) == getattr(
            model_2.metabolites[0], attr)
        assert getattr(model_1.metabolites[5], attr) == getattr(
            model_2.metabolites[5], attr)
        assert getattr(model_1.metabolites[-1], attr) == getattr(
            model_2.metabolites[-1], attr)
        assert len(model_1.reactions[0].metabolites) == len(
            model_2.reactions[0].metabolites)
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
        model_2.reactions[8].metabolites)
    assert len(model_1.reactions[-1].metabolites) == len(
        model_2.reactions[-1].metabolites)
    assert len(model_1.genes) == len(model_2.genes)

    # ensure they have the same solution max
    solution_1 = model_1.optimize()
    solution_2 = model_2.optimize()
    assert abs(solution_1.objective_value -
               solution_2.objective_value) == pytest.approx(0.0)

    # ensure the references are correct
    # metabolite -> model reference
    assert model_1.metabolites[0]._model is model_1
    assert model_2.metabolites[0]._model is model_2

    # reaction -> model reference
    assert model_1.reactions[0]._model is model_1
    assert model_2.reactions[0]._model is model_2

    # gene -> model reference
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
