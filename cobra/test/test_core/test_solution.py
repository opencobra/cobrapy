# -*- coding: utf-8 -*-

"""Test functions of solution.py"""

from __future__ import absolute_import

from cobra.core import Solution


def test_solution_contains_only_reaction_specific_values(solved_model):
    solution, model = solved_model
    reaction_ids = set([reaction.id for reaction in model.reactions])
    if isinstance(solution, Solution):
        assert set(solution.fluxes.index) == reaction_ids
        # assert set(solution.reduced_costs.index) == reaction_ids
    else:
        raise TypeError(
            "solutions of type {0:r} are untested".format(type(solution)))
