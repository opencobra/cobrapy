# -*- coding: utf-8 -*-

from __future__ import absolute_import

import copy
import os

import numpy
import optlang
import pytest
import six

import cobra
from cobra.core import Metabolite, Model, Reaction, Solution
from cobra.util.solver import SolverNotFound, set_objective, solvers

solver_trials = ['glpk',
                 pytest.mark.skipif('cplex' not in solvers,
                                    reason='no cplex')]


@pytest.fixture(scope="function", params=solver_trials)
def solved_model(request, model):
    model.solver = request.param
    solution = model.optimize()
    return solution, model


def same_ex(ex1, ex2):
    """Compare to expressions for mathematical equality."""
    return ex1.simplify() == ex2.simplify()


class TestSolution:
    def test_solution_contains_only_reaction_specific_values(self,
                                                             solved_model):
        solution, model = solved_model
        reaction_ids = set([reaction.id for reaction in model.reactions])
        if isinstance(solution, Solution):
            assert set(solution.fluxes.index) == reaction_ids
#            assert set(solution.reduced_costs.index) == reaction_ids
        else:
            raise TypeError(
                "solutions of type {0:r} are untested".format(type(solution)))


class TestReaction:
    def test_str(self, model):
        assert model.reactions[0].__str__().startswith('ACALD')

    def test_add_metabolite(self, solved_model):
        solution, model = solved_model
        pgi_reaction = model.reactions.PGI
        test_met = model.metabolites[0]
        pgi_reaction.add_metabolites({test_met: 42}, combine=False)
        assert pgi_reaction.metabolites[test_met] == 42.0
        assert model.constraints[
                   test_met.id].expression.as_coefficients_dict()[
                   pgi_reaction.forward_variable] == 42.0
        assert model.constraints[
                   test_met.id].expression.as_coefficients_dict()[
                   pgi_reaction.reverse_variable] == -42.0

        pgi_reaction.add_metabolites({test_met: -10}, combine=True)
        assert pgi_reaction.metabolites[test_met] == 32.0
        assert model.constraints[
                   test_met.id].expression.as_coefficients_dict()[
                   pgi_reaction.forward_variable] == 32.0
        assert model.constraints[
                   test_met.id].expression.as_coefficients_dict()[
                   pgi_reaction.reverse_variable] == -32.0

        pgi_reaction.add_metabolites({test_met: 0}, combine=False)
        with pytest.raises(KeyError):
            pgi_reaction.metabolites[test_met]
        assert model.constraints[
                   test_met.id].expression.as_coefficients_dict()[
                   pgi_reaction.forward_variable] == 0
        assert model.constraints[
                   test_met.id].expression.as_coefficients_dict()[
                   pgi_reaction.reverse_variable] == 0

    def test_removal_from_model_retains_bounds(self, model):
        model_cp = model.copy()
        reaction = model_cp.reactions.ACALD
        assert reaction.model == model_cp
        assert reaction.lower_bound == -1000.0
        assert reaction.upper_bound == 1000.0
        assert reaction._lower_bound == -1000.0
        assert reaction._upper_bound == 1000.0
        model_cp.remove_reactions([reaction])
        assert reaction.model is None
        assert reaction.lower_bound == -1000.0
        assert reaction.upper_bound == 1000.0
        assert reaction._lower_bound == -1000.0
        assert reaction._upper_bound == 1000.0

    def test_set_bounds_scenario_1(self, model):
        acald_reaction = model.reactions.ACALD
        assert acald_reaction.lower_bound == -1000.
        assert acald_reaction.upper_bound == 1000.
        assert acald_reaction.forward_variable.lb == 0.
        assert acald_reaction.forward_variable.ub == 1000.
        assert acald_reaction.reverse_variable.lb == 0
        assert acald_reaction.reverse_variable.ub == 1000.
        acald_reaction.upper_bound = acald_reaction.lower_bound - 100
        assert acald_reaction.lower_bound == -1100.0
        assert acald_reaction.upper_bound == -1100.0
        assert acald_reaction.forward_variable.lb == 0
        assert acald_reaction.forward_variable.ub == 0
        assert acald_reaction.reverse_variable.lb == 1100.
        assert acald_reaction.reverse_variable.ub == 1100.
        acald_reaction.upper_bound = 100
        assert acald_reaction.lower_bound == -1100.0
        assert acald_reaction.upper_bound == 100
        assert acald_reaction.forward_variable.lb == 0
        assert acald_reaction.forward_variable.ub == 100
        assert acald_reaction.reverse_variable.lb == 0
        assert acald_reaction.reverse_variable.ub == 1100.0

    def test_set_bounds_scenario_3(self, model):
        reac = model.reactions.ACALD
        reac.upper_bound = -10
        reac.lower_bound = -10
        assert reac.lower_bound == -10
        assert reac.upper_bound == -10
        reac.lower_bound = -9
        assert reac.lower_bound == -9
        assert reac.upper_bound == -9
        reac.lower_bound = 2
        assert reac.lower_bound == 2
        assert reac.upper_bound == 2
        reac.upper_bound = -10
        assert reac.lower_bound == -10
        assert reac.upper_bound == -10
        reac.upper_bound = -11
        assert reac.lower_bound == -11
        assert reac.upper_bound == -11
        reac.upper_bound = 2
        assert reac.lower_bound == -11
        assert reac.upper_bound == 2

    def test_set_bounds_scenario_4(self, model):
        reac = model.reactions.ACALD
        reac.lower_bound = reac.upper_bound = 0
        reac.lower_bound = 2
        assert reac.lower_bound == 2
        assert reac.upper_bound == 2
        assert reac.forward_variable.lb == 2
        assert reac.forward_variable.ub == 2
        reac.knock_out()
        reac.upper_bound = -2
        assert reac.lower_bound == -2
        assert reac.upper_bound == -2
        assert reac.reverse_variable.lb == 2
        assert reac.reverse_variable.ub == 2

    def test_set_upper_before_lower_bound_to_0(self, model):
        model.reactions.GAPD.upper_bound = 0
        model.reactions.GAPD.lower_bound = 0
        assert model.reactions.GAPD.lower_bound == 0
        assert model.reactions.GAPD.upper_bound == 0
        assert model.reactions.GAPD.forward_variable.lb == 0
        assert model.reactions.GAPD.forward_variable.ub == 0
        assert model.reactions.GAPD.reverse_variable.lb == 0
        assert model.reactions.GAPD.reverse_variable.ub == 0

    def test_set_bounds_scenario_2(self, model):
        acald_reaction = model.reactions.ACALD
        assert acald_reaction.lower_bound == -1000.
        assert acald_reaction.upper_bound == 1000.
        assert acald_reaction.forward_variable.lb == 0.
        assert acald_reaction.forward_variable.ub == 1000.
        assert acald_reaction.reverse_variable.lb == 0
        assert acald_reaction.reverse_variable.ub == 1000.
        acald_reaction.lower_bound = acald_reaction.upper_bound + 100
        assert acald_reaction.lower_bound == 1100.0
        assert acald_reaction.upper_bound == 1100.0
        assert acald_reaction.forward_variable.lb == 1100.0
        assert acald_reaction.forward_variable.ub == 1100.0
        assert acald_reaction.reverse_variable.lb == 0
        assert acald_reaction.reverse_variable.ub == 0
        acald_reaction.lower_bound = -100
        assert acald_reaction.lower_bound == -100.
        assert acald_reaction.upper_bound == 1100.
        assert acald_reaction.forward_variable.lb == 0
        assert acald_reaction.forward_variable.ub == 1100.
        assert acald_reaction.reverse_variable.lb == 0
        assert acald_reaction.reverse_variable.ub == 100

    def test_change_bounds(self, model):
        reac = model.reactions.ACALD
        reac.bounds = (2, 2)
        assert reac.lower_bound == 2
        assert reac.upper_bound == 2
        with model:
            reac.bounds = (5, 5)
            assert reac.lower_bound == 5
            assert reac.upper_bound == 5
        assert reac.lower_bound == 2
        assert reac.upper_bound == 2

    def test_make_irreversible(self, model):
        acald_reaction = model.reactions.ACALD
        assert acald_reaction.lower_bound == -1000.
        assert acald_reaction.upper_bound == 1000.
        assert acald_reaction.forward_variable.lb == 0.
        assert acald_reaction.forward_variable.ub == 1000.
        assert acald_reaction.reverse_variable.lb == 0
        assert acald_reaction.reverse_variable.ub == 1000.
        acald_reaction.lower_bound = 0
        assert acald_reaction.lower_bound == 0
        assert acald_reaction.upper_bound == 1000.
        assert acald_reaction.forward_variable.lb == 0
        assert acald_reaction.forward_variable.ub == 1000.0
        assert acald_reaction.reverse_variable.lb == 0
        assert acald_reaction.reverse_variable.ub == 0
        acald_reaction.lower_bound = -100
        assert acald_reaction.lower_bound == -100.
        assert acald_reaction.upper_bound == 1000.
        assert acald_reaction.forward_variable.lb == 0
        assert acald_reaction.forward_variable.ub == 1000.
        assert acald_reaction.reverse_variable.lb == 0
        assert acald_reaction.reverse_variable.ub == 100

    def test_make_reversible(self, model):
        pfk_reaction = model.reactions.PFK
        assert pfk_reaction.lower_bound == 0.
        assert pfk_reaction.upper_bound == 1000.
        assert pfk_reaction.forward_variable.lb == 0.
        assert pfk_reaction.forward_variable.ub == 1000.
        assert pfk_reaction.reverse_variable.lb == 0
        assert pfk_reaction.reverse_variable.ub == 0
        pfk_reaction.lower_bound = -100.
        assert pfk_reaction.lower_bound == -100.
        assert pfk_reaction.upper_bound == 1000.
        assert pfk_reaction.forward_variable.lb == 0
        assert pfk_reaction.forward_variable.ub == 1000.0
        assert pfk_reaction.reverse_variable.lb == 0
        assert pfk_reaction.reverse_variable.ub == 100.
        pfk_reaction.lower_bound = 0
        assert pfk_reaction.lower_bound == 0
        assert pfk_reaction.upper_bound == 1000.
        assert pfk_reaction.forward_variable.lb == 0
        assert pfk_reaction.forward_variable.ub == 1000.
        assert pfk_reaction.reverse_variable.lb == 0
        assert pfk_reaction.reverse_variable.ub == 0

    def test_make_irreversible_irreversible_to_the_other_side(self, model):
        pfk_reaction = model.reactions.PFK
        assert pfk_reaction.lower_bound == 0.
        assert pfk_reaction.upper_bound == 1000.
        assert pfk_reaction.forward_variable.lb == 0.
        assert pfk_reaction.forward_variable.ub == 1000.
        assert pfk_reaction.reverse_variable.lb == 0
        assert pfk_reaction.reverse_variable.ub == 0
        pfk_reaction.upper_bound = -100.
        assert pfk_reaction.forward_variable.lb == 0
        assert pfk_reaction.forward_variable.ub == 0
        assert pfk_reaction.reverse_variable.lb == 100
        assert pfk_reaction.reverse_variable.ub == 100
        pfk_reaction.lower_bound = -1000.
        assert pfk_reaction.lower_bound == -1000.
        assert pfk_reaction.upper_bound == -100.
        assert pfk_reaction.forward_variable.lb == 0
        assert pfk_reaction.forward_variable.ub == 0
        assert pfk_reaction.reverse_variable.lb == 100
        assert pfk_reaction.reverse_variable.ub == 1000.

    def test_make_lhs_irreversible_reversible(self, model):
        rxn = Reaction('test')
        rxn.add_metabolites(
            {model.metabolites[0]: -1., model.metabolites[1]: 1.})
        rxn.lower_bound = -1000.
        rxn.upper_bound = -100
        model.add_reaction(rxn)
        assert rxn.lower_bound == -1000.
        assert rxn.upper_bound == -100.
        assert rxn.forward_variable.lb == 0.
        assert rxn.forward_variable.ub == 0.
        assert rxn.reverse_variable.lb == 100.
        assert rxn.reverse_variable.ub == 1000.
        rxn.upper_bound = 666.
        assert rxn.lower_bound == -1000.
        assert rxn.upper_bound == 666.
        assert rxn.forward_variable.lb == 0.
        assert rxn.forward_variable.ub == 666
        assert rxn.reverse_variable.lb == 0.
        assert rxn.reverse_variable.ub == 1000.

    def test_model_less_reaction(self, model):
        model.slim_optimize()
        for reaction in model.reactions:
            assert isinstance(reaction.flux, float)
            assert isinstance(reaction.reduced_cost, float)
        for reaction in model.reactions:
            model.remove_reactions([reaction])
            with pytest.raises(RuntimeError):
                reaction.flux
            with pytest.raises(RuntimeError):
                reaction.reduced_cost

    def test_knockout(self, model):
        original_bounds = dict()
        for reaction in model.reactions:
            original_bounds[reaction.id] = (
                reaction.lower_bound, reaction.upper_bound)
            reaction.knock_out()
            assert reaction.lower_bound == 0
            assert reaction.upper_bound == 0
        for k, (lb, ub) in six.iteritems(original_bounds):
            model.reactions.get_by_id(k).lower_bound = lb
            model.reactions.get_by_id(k).upper_bound = ub
        for reaction in model.reactions:
            assert reaction.lower_bound == original_bounds[reaction.id][0]
            assert reaction.upper_bound == original_bounds[reaction.id][1]
        with model:
            for reaction in model.reactions:
                original_bounds[reaction.id] = (
                    reaction.lower_bound, reaction.upper_bound)
                reaction.knock_out()
                assert reaction.lower_bound == 0
                assert reaction.upper_bound == 0
        for reaction in model.reactions:
            assert reaction.lower_bound == original_bounds[reaction.id][0]
            assert reaction.upper_bound == original_bounds[reaction.id][1]

    def test_reaction_without_model(self):
        r = Reaction('blub')
        assert r.flux_expression is None
        assert r.forward_variable is None
        assert r.reverse_variable is None

    def test_weird_left_to_right_reaction_issue(self, tiny_toy_model):
        d1 = tiny_toy_model.reactions.get_by_id('ex1')
        assert not d1.reversibility
        assert d1.lower_bound == -1000
        assert d1._lower_bound == -1000
        assert d1.upper_bound == 0
        assert d1._upper_bound == 0
        with tiny_toy_model:
            d1.knock_out()
            assert d1.lower_bound == 0
            assert d1._lower_bound == 0
            assert d1.upper_bound == 0
            assert d1._upper_bound == 0
        assert d1.lower_bound == -1000
        assert d1._lower_bound == -1000
        assert d1.upper_bound == 0
        assert d1._upper_bound == 0

    def test_one_left_to_right_reaction_set_positive_ub(self, tiny_toy_model):
        d1 = tiny_toy_model.reactions.get_by_id('ex1')
        assert d1.reverse_variable.lb == 0
        assert d1.reverse_variable.ub == 1000
        assert d1._lower_bound == -1000
        assert d1.lower_bound == -1000
        assert d1._upper_bound == 0
        assert d1.upper_bound == 0
        assert d1.forward_variable.lb == 0
        assert d1.forward_variable.ub == 0
        d1.upper_bound = .1
        assert d1.forward_variable.lb == 0
        assert d1.forward_variable.ub == .1
        assert d1.reverse_variable.lb == 0
        assert d1.reverse_variable.ub == 1000
        assert d1._lower_bound == -1000
        assert d1.upper_bound == .1
        assert d1._lower_bound == -1000
        assert d1.upper_bound == .1

    def test_irrev_reaction_set_negative_lb(self, model):
        assert not model.reactions.PFK.reversibility
        assert model.reactions.PFK.lower_bound == 0
        assert model.reactions.PFK.upper_bound == 1000.0
        assert model.reactions.PFK.forward_variable.lb == 0
        assert model.reactions.PFK.forward_variable.ub == 1000.0
        assert model.reactions.PFK.reverse_variable.lb == 0
        assert model.reactions.PFK.reverse_variable.ub == 0
        model.reactions.PFK.lower_bound = -1000
        assert model.reactions.PFK.lower_bound == -1000
        assert model.reactions.PFK.upper_bound == 1000.0
        assert model.reactions.PFK.forward_variable.lb == 0
        assert model.reactions.PFK.forward_variable.ub == 1000.0
        assert model.reactions.PFK.reverse_variable.lb == 0
        assert model.reactions.PFK.reverse_variable.ub == 1000

    def test_twist_irrev_right_to_left_reaction_to_left_to_right(self, model):
        assert not model.reactions.PFK.reversibility
        assert model.reactions.PFK.lower_bound == 0
        assert model.reactions.PFK.upper_bound == 1000.0
        assert model.reactions.PFK.forward_variable.lb == 0
        assert model.reactions.PFK.forward_variable.ub == 1000.0
        assert model.reactions.PFK.reverse_variable.lb == 0
        assert model.reactions.PFK.reverse_variable.ub == 0
        model.reactions.PFK.lower_bound = -1000
        model.reactions.PFK.upper_bound = 0
        assert model.reactions.PFK.lower_bound == -1000
        assert model.reactions.PFK.upper_bound == 0
        assert model.reactions.PFK.forward_variable.lb == 0
        assert model.reactions.PFK.forward_variable.ub == 0
        assert model.reactions.PFK.reverse_variable.lb == 0
        assert model.reactions.PFK.reverse_variable.ub == 1000

    def test_set_lb_higher_than_ub_sets_ub_to_new_lb(self, model):
        for reaction in model.reactions:
            assert reaction.lower_bound <= reaction.upper_bound
            reaction.lower_bound = reaction.upper_bound + 100
            assert reaction.lower_bound == reaction.upper_bound

    def test_set_ub_lower_than_lb_sets_lb_to_new_ub(self, model):
        for reaction in model.reactions:
            assert reaction.lower_bound <= reaction.upper_bound
            reaction.upper_bound = reaction.lower_bound - 100
            assert reaction.lower_bound == reaction.upper_bound

    def test_add_metabolites_combine_true(self, model):
        test_metabolite = Metabolite('test')
        for reaction in model.reactions:
            reaction.add_metabolites({test_metabolite: -66}, combine=True)
            assert reaction.metabolites[test_metabolite] == -66
            assert model.constraints['test'].get_linear_coefficients(
                [reaction.forward_variable])[reaction.forward_variable] == -66
            assert model.constraints['test'].get_linear_coefficients(
                [reaction.reverse_variable])[reaction.reverse_variable] == 66
            already_included_metabolite = \
                list(reaction.metabolites.keys())[0]
            previous_coefficient = reaction.get_coefficient(
                already_included_metabolite.id)
            reaction.add_metabolites({already_included_metabolite: 10},
                                     combine=True)
            new_coefficient = previous_coefficient + 10
            assert reaction.metabolites[
                       already_included_metabolite] == new_coefficient
            assert (model.constraints[
                    already_included_metabolite.id].get_linear_coefficients(
                    [reaction.forward_variable])[reaction.forward_variable] ==
                    new_coefficient)
            assert (model.constraints[
                    already_included_metabolite.id].get_linear_coefficients(
                    [reaction.reverse_variable])[
                        reaction.reverse_variable] == -new_coefficient)

    @pytest.mark.xfail(reason='non-deterministic test')
    def test_add_metabolites_combine_false(self, model):
        test_metabolite = Metabolite('test')
        for reaction in model.reactions:
            reaction.add_metabolites({test_metabolite: -66}, combine=False)
            assert reaction.metabolites[test_metabolite] == -66
            assert model.constraints['test'].expression.has(
                -66. * reaction.forward_variable)
            assert model.constraints['test'].expression.has(
                66. * reaction.reverse_variable)
            already_included_metabolite = \
                list(reaction.metabolites.keys())[0]
            reaction.add_metabolites({already_included_metabolite: 10},
                                     combine=False)
            assert reaction.metabolites[already_included_metabolite] == 10
            assert model.constraints[
                already_included_metabolite.id].expression.has(
                10 * reaction.forward_variable)
            assert model.constraints[
                already_included_metabolite.id].expression.has(
                -10 * reaction.reverse_variable)

    def test_reaction_imul(self, model):
        with model:
            model.reactions.EX_glc__D_e *= 100
            assert model.constraints.glc__D_e.expression.coeff(
                model.variables.EX_glc__D_e) == -100
            assert model.reactions.EX_glc__D_e.reaction == \
                '100.0 glc__D_e <=> '

        assert model.constraints.glc__D_e.expression.coeff(
            model.variables.EX_glc__D_e) == -1
        assert model.reactions.EX_glc__D_e.reaction == 'glc__D_e <=> '

        with model:
            model.reactions.EX_glc__D_e *= -2
            assert model.reactions.EX_glc__D_e.bounds == (-1000.0, 10.0)
            assert model.reactions.EX_glc__D_e.reaction == ' <=> 2.0 glc__D_e'

        assert model.reactions.EX_glc__D_e.bounds == (-10, 1000.0)
        assert model.reactions.EX_glc__D_e.reaction == 'glc__D_e <=> '

    # def test_pop(self, model):
    #     pgi = model.reactions.PGI
    #     g6p = model.metabolites.get_by_id("g6p_c")
    #     f6p = model.metabolites.get_by_id("f6p_c")
    #     g6p_expr = model.solver.constraints["g6p_c"].expression
    #     g6p_coef = pgi.pop("g6p_c")
    #     assert g6p not in pgi.metabolites
    #     actual = model.solver.constraints[
    #         "g6p_c"].expression.as_coefficients_dict()
    #     expected = (g6p_expr - g6p_coef * pgi.flux_expression
    #                 ).as_coefficients_dict()
    #     assert actual == expected
    #     assert pgi.metabolites[f6p] == 1
    #
    #     f6p_expr = model.solver.constraints["f6p_c"].expression
    #     f6p_coef = pgi.pop(f6p)
    #     assert f6p not in pgi.metabolites
    #     assert model.solver.constraints[
    #                "f6p_c"].expression.as_coefficients_dict() == (
    #                f6p_expr - f6p_coef * pgi.flux_expression
    #            ).as_coefficients_dict()

    def test_remove_from_model(self, model):
        pgi = model.reactions.PGI
        g6p = model.metabolites.g6p_c
        pgi_flux = model.optimize().fluxes['PGI']
        assert abs(pgi_flux) > 1E-6

        with model:
            pgi.remove_from_model()
            assert pgi.model is None
            assert "PGI" not in model.reactions
            assert pgi.id not in model.variables
            assert pgi.reverse_id not in model.variables
            assert pgi not in g6p.reactions
            model.optimize()

        assert "PGI" in model.reactions
        assert pgi.id in model.variables
        assert pgi.reverse_id in model.variables
        assert pgi.forward_variable.problem is model.solver
        assert pgi in g6p.reactions
        assert g6p in pgi.metabolites
        assert numpy.isclose(pgi_flux, model.optimize().fluxes['PGI'])

    def test_change_id_is_reflected_in_solver(self, model):
        for i, reaction in enumerate(model.reactions):
            old_reaction_id = reaction.id
            assert model.variables[
                       old_reaction_id].name == old_reaction_id
            assert old_reaction_id in model.variables
            new_reaction_id = reaction.id + '_' + str(i)
            reaction.id = new_reaction_id
            assert reaction.id == new_reaction_id
            assert not (old_reaction_id in model.variables)
            assert reaction.id in model.variables
            assert reaction.reverse_id in model.variables
            name = model.variables[reaction.id].name
            assert name == reaction.id


class TestSolverBasedModel:
    def test_objective_coefficient_reflects_changed_objective(self, model):
        biomass_r = model.reactions.get_by_id('Biomass_Ecoli_core')
        assert biomass_r.objective_coefficient == 1
        model.objective = "PGI"
        assert biomass_r.objective_coefficient == 0
        assert model.reactions.PGI.objective_coefficient == 1

    def test_change_objective_through_objective_coefficient(self, model):
        biomass_r = model.reactions.get_by_id('Biomass_Ecoli_core')
        pgi = model.reactions.PGI
        pgi.objective_coefficient = 2
        coef_dict = model.objective.expression.as_coefficients_dict()
        # Check that objective has been updated
        assert coef_dict[pgi.forward_variable] == 2.0
        assert coef_dict[pgi.reverse_variable] == -2.0
        # Check that original objective is still in there
        assert coef_dict[biomass_r.forward_variable] == 1.0
        assert coef_dict[biomass_r.reverse_variable] == -1.0

    def test_transfer_objective(self, model):
        new_mod = Model("new model")
        new_mod.add_reactions(model.reactions)
        new_mod.objective = model.objective
        assert (set(str(x) for x in model.objective.expression.args) == set(
            str(x) for x in new_mod.objective.expression.args))
        new_mod.slim_optimize()
        assert abs(new_mod.objective.value - 0.874) < 0.001

    def test_model_from_other_model(self, model):
        model = Model(id_or_model=model)
        for reaction in model.reactions:
            assert reaction == model.reactions.get_by_id(reaction.id)

    def test_add_reactions(self, model):
        r1 = Reaction('r1')
        r1.add_metabolites({Metabolite('A'): -1, Metabolite('B'): 1})
        r1.lower_bound, r1.upper_bound = -999999., 999999.
        r2 = Reaction('r2')
        r2.add_metabolites(
            {Metabolite('A'): -1, Metabolite('C'): 1, Metabolite('D'): 1})
        r2.lower_bound, r2.upper_bound = 0., 999999.
        model.add_reactions([r1, r2])
        r2.objective_coefficient = 3.
        assert r2.objective_coefficient == 3.
        assert model.reactions[-2] == r1
        assert model.reactions[-1] == r2
        assert isinstance(model.reactions[-2].reverse_variable,
                          model.problem.Variable)
        coefficients_dict = model.objective.expression. \
            as_coefficients_dict()
        biomass_r = model.reactions.get_by_id('Biomass_Ecoli_core')
        assert coefficients_dict[biomass_r.forward_variable] == 1.
        assert coefficients_dict[biomass_r.reverse_variable] == -1.
        assert coefficients_dict[
                   model.reactions.r2.forward_variable] == 3.
        assert coefficients_dict[
                   model.reactions.r2.reverse_variable] == -3.

    def test_add_reactions_single_existing(self, model):
        rxn = model.reactions[0]
        r1 = Reaction(rxn.id)
        r1.add_metabolites({Metabolite('A'): -1, Metabolite('B'): 1})
        r1.lower_bound, r1.upper_bound = -999999., 999999.
        model.add_reactions([r1])
        assert rxn in model.reactions
        assert r1 is not model.reactions.get_by_id(rxn.id)

    def test_add_reactions_duplicate(self, model):
        rxn = model.reactions[0]
        r1 = Reaction('r1')
        r1.add_metabolites({Metabolite('A'): -1, Metabolite('B'): 1})
        r1.lower_bound, r1.upper_bound = -999999., 999999.
        r2 = Reaction(rxn.id)
        r2.add_metabolites(
            {Metabolite('A'): -1, Metabolite('C'): 1, Metabolite('D'): 1})
        model.add_reactions([r1, r2])
        assert r1 in model.reactions
        assert rxn in model.reactions
        assert r2 is not model.reactions.get_by_id(rxn.id)

    def test_add_cobra_reaction(self, model):
        r = cobra.Reaction(id="c1")
        model.add_reaction(r)
        assert isinstance(model.reactions.c1, Reaction)

    def test_all_objects_point_to_all_other_correct_objects(self, model):
        for reaction in model.reactions:
            assert reaction.model == model
            for gene in reaction.genes:
                assert gene == model.genes.get_by_id(gene.id)
                assert gene.model == model
                for reaction2 in gene.reactions:
                    assert reaction2.model == model
                    assert reaction2 == model.reactions.get_by_id(
                        reaction2.id)

            for metabolite in reaction.metabolites:
                assert metabolite.model == model
                assert metabolite == model.metabolites.get_by_id(
                    metabolite.id)
                for reaction2 in metabolite.reactions:
                    assert reaction2.model == model
                    assert reaction2 == model.reactions.get_by_id(
                        reaction2.id)

    def test_objects_point_to_correct_other_after_copy(self, model):
        for reaction in model.reactions:
            assert reaction.model == model
            for gene in reaction.genes:
                assert gene == model.genes.get_by_id(gene.id)
                assert gene.model == model
                for reaction2 in gene.reactions:
                    assert reaction2.model == model
                    assert reaction2 == model.reactions.get_by_id(
                        reaction2.id)

            for metabolite in reaction.metabolites:
                assert metabolite.model == model
                assert metabolite == model.metabolites.get_by_id(
                    metabolite.id)
                for reaction2 in metabolite.reactions:
                    assert reaction2.model == model
                    assert reaction2 == model.reactions.get_by_id(
                        reaction2.id)

    def test_remove_reactions(self, model):
        reactions_to_remove = model.reactions[10:30]
        assert all([reaction.model is model for reaction in
                    reactions_to_remove])
        assert all(
            [model.reactions.get_by_id(reaction.id) == reaction for
             reaction in reactions_to_remove])

        model.remove_reactions(reactions_to_remove)
        assert all(
            [reaction.model is None for reaction in reactions_to_remove])
        for reaction in reactions_to_remove:
            assert reaction.id not in list(
                model.variables.keys())

        model.add_reactions(reactions_to_remove)
        for reaction in reactions_to_remove:
            assert reaction in model.reactions

    def test_objective(self, model):
        obj = model.objective
        assert obj.get_linear_coefficients(obj.variables) == {
                   model.variables["Biomass_Ecoli_core_reverse_2cdba"]: -1,
                   model.variables["Biomass_Ecoli_core"]: 1}
        assert obj.direction == "max"

    def test_change_objective(self, model):
        expression = 1.0 * model.variables['ENO'] + \
                     1.0 * model.variables['PFK']
        model.objective = model.problem.Objective(
            expression)
        assert same_ex(model.objective.expression, expression)
        model.objective = "ENO"
        eno_obj = model.problem.Objective(
            model.reactions.ENO.flux_expression, direction="max")
        pfk_obj = model.problem.Objective(
            model.reactions.PFK.flux_expression, direction="max")
        assert same_ex(model.objective.expression, eno_obj.expression)

        with model:
            model.objective = "PFK"
            assert same_ex(model.objective.expression, pfk_obj.expression)
        assert same_ex(model.objective.expression, eno_obj.expression)
        expression = model.objective.expression
        atpm = model.reactions.get_by_id("ATPM")
        biomass = model.reactions.get_by_id("Biomass_Ecoli_core")
        with model:
            model.objective = atpm
        assert same_ex(model.objective.expression, expression)
        with model:
            atpm.objective_coefficient = 1
            biomass.objective_coefficient = 2
        assert same_ex(model.objective.expression, expression)

        with model:
            set_objective(model, model.problem.Objective(
                atpm.flux_expression))
            assert same_ex(model.objective.expression, atpm.flux_expression)
        assert same_ex(model.objective.expression, expression)

        expression = model.objective.expression
        with model:
            with model:  # Test to make sure nested contexts are OK
                set_objective(model, atpm.flux_expression,
                              additive=True)
                assert same_ex(model.objective.expression,
                               expression + atpm.flux_expression)
        assert same_ex(model.objective.expression, expression)

    def test_set_reaction_objective(self, model):
        model.objective = model.reactions.ACALD
        assert str(model.objective.expression) == str(
            1.0 * model.reactions.ACALD.forward_variable -
            1.0 * model.reactions.ACALD.reverse_variable)

    def test_set_reaction_objective_str(self, model):
        model.objective = model.reactions.ACALD.id
        assert str(model.objective.expression) == str(
            1.0 * model.reactions.ACALD.forward_variable -
            1.0 * model.reactions.ACALD.reverse_variable)

    def test_invalid_objective_raises(self, model):
        with pytest.raises(ValueError):
            setattr(model, 'objective', 'This is not a valid objective!')
        with pytest.raises(TypeError):
            setattr(model, 'objective', 3.)

    @pytest.mark.skipif("cplex" not in solvers, reason="need cplex")
    def test_solver_change(self, model):
        solver_id = id(model.solver)
        problem_id = id(model.solver.problem)
        solution = model.optimize().fluxes
        model.solver = "cplex"
        assert id(model.solver) != solver_id
        assert id(model.solver.problem) != problem_id
        new_solution = model.optimize().fluxes
        assert numpy.allclose(solution, new_solution, rtol=0, atol=1E-06)

    def test_no_change_for_same_solver(self, model):
        solver_id = id(model.solver)
        problem_id = id(model.solver.problem)
        model.solver = "glpk"
        assert id(model.solver) == solver_id
        assert id(model.solver.problem) == problem_id

    def test_invalid_solver_change_raises(self, model):
        with pytest.raises(SolverNotFound):
            setattr(model, 'solver', [1, 2, 3])
        with pytest.raises(SolverNotFound):
            setattr(model, 'solver',
                    'ThisIsDefinitelyNotAvalidSolver')
        with pytest.raises(SolverNotFound):
            setattr(model, 'solver', os)

    @pytest.mark.skipif('cplex' not in solvers, reason='no cplex')
    def test_change_solver_to_cplex_and_check_copy_works(self, model):
        assert round(abs(model.optimize().f - 0.8739215069684306), 7) == 0
        model_copy = model.copy()
        assert round(abs(model_copy.optimize().f - 0.8739215069684306),
                     7) == 0
        # Second, change existing glpk based model to cplex
        model.solver = 'cplex'
        assert round(abs(model.optimize().f - 0.8739215069684306),
                     7) == 0
        model_copy = copy.copy(model)
        assert round(abs(model_copy.optimize().f - 0.8739215069684306),
                     7) == 0

    def test_copy_preserves_existing_solution(self, solved_model):
        solution, model = solved_model
        model_cp = copy.copy(model)
        primals_original = [variable.primal for variable in
                            model.variables]
        primals_copy = [variable.primal for variable in
                        model_cp.variables]
        abs_diff = abs(
            numpy.array(primals_copy) - numpy.array(primals_original))
        assert not any(abs_diff > 1e-6)


class TestMetabolite:
    def test_set_id(self, solved_model):
        solution, model = solved_model
        met = Metabolite("test")
        with pytest.raises(TypeError):
            setattr(met, 'id', 1)
        model.add_metabolites([met])
        with pytest.raises(ValueError):
            setattr(met, "id", 'g6p_c')
        met.id = "test2"
        assert "test2" in model.metabolites
        assert "test" not in model.metabolites

    def test_remove_from_model(self, solved_model):
        solution, model = solved_model
        met = model.metabolites.get_by_id("g6p_c")
        met.remove_from_model()
        assert not (met.id in model.metabolites)
        assert not (met.id in model.constraints)
