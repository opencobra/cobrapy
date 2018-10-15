# -*- coding: utf-8 -*-

"""Test functions of model.py"""

from __future__ import absolute_import

import warnings
from copy import copy, deepcopy
from math import isnan
import os

import numpy as np
import pandas as pd
import pytest
from optlang.symbolics import Zero

import cobra.util.solver as su
from cobra.core import Metabolite, Model, Reaction
from cobra.util.solver import SolverNotFound, set_objective, solvers
from cobra.exceptions import OptimizationError

stable_optlang = ["glpk", "cplex", "gurobi"]
optlang_solvers = ["optlang-" + s for s in stable_optlang if s in su.solvers]


def same_ex(ex1, ex2):
    """Compare to expressions for mathematical equality."""
    return ex1.simplify() == ex2.simplify()


@pytest.mark.parametrize("solver", stable_optlang)
def test_add_remove_reaction_benchmark(model, benchmark, solver):
    metabolite_foo = Metabolite("test_foo")
    metabolite_bar = Metabolite("test_bar")
    metabolite_baz = Metabolite("test_baz")
    actual_metabolite = model.metabolites[0]
    dummy_reaction = Reaction("test_foo_reaction")
    dummy_reaction.add_metabolites({metabolite_foo: -1,
                                    metabolite_bar: 1,
                                    metabolite_baz: -2,
                                    actual_metabolite: 1})

    def benchmark_add_reaction():
        model.add_reaction(dummy_reaction)
        if not getattr(model, 'solver', None):
            solver_dict[solver].create_problem(model)
        model.remove_reactions([dummy_reaction])

    benchmark(benchmark_add_reaction)


def test_add_metabolite(model):
    new_metabolite = Metabolite('test_met')
    assert new_metabolite not in model.metabolites
    with model:
        model.add_metabolites(new_metabolite)
        assert new_metabolite._model == model
        assert new_metabolite in model.metabolites
        assert new_metabolite.id in model.solver.constraints

    assert new_metabolite._model is None
    assert new_metabolite not in model.metabolites
    assert new_metabolite.id not in model.solver.constraints


def test_remove_metabolite_subtractive(model):
    test_metabolite = model.metabolites[4]
    test_reactions = test_metabolite.reactions
    with model:
        model.remove_metabolites(test_metabolite, destructive=False)
        assert test_metabolite._model is None
        assert test_metabolite not in model.metabolites
        assert test_metabolite.id not in model.solver.constraints
        for reaction in test_reactions:
            assert reaction in model.reactions

    assert test_metabolite._model is model
    assert test_metabolite in model.metabolites
    assert test_metabolite.id in model.solver.constraints


def test_remove_metabolite_destructive(model):
    test_metabolite = model.metabolites[4]
    test_reactions = test_metabolite.reactions
    with model:
        model.remove_metabolites(test_metabolite, destructive=True)
        assert test_metabolite._model is None
        assert test_metabolite not in model.metabolites
        assert test_metabolite.id not in model.solver.constraints
        for reaction in test_reactions:
            assert reaction not in model.reactions

    assert test_metabolite._model is model
    assert test_metabolite in model.metabolites
    assert test_metabolite.id in model.solver.constraints
    for reaction in test_reactions:
        assert reaction in model.reactions


def test_compartments(model):
    assert set(model.compartments) == {"c", "e"}
    model = Model("test", "test")
    met_c = Metabolite("a_c", compartment="c")
    met_e = Metabolite("a_e", compartment="e")
    rxn = Reaction("foo")
    rxn.add_metabolites({met_e: -1, met_c: 1})
    model.add_reactions([rxn])
    assert model.compartments == {'c': '', 'e': ''}
    model.compartments = {'c': 'cytosol'}
    assert model.compartments == {'c': 'cytosol', 'e': ''}


def test_add_reaction(model):
    old_reaction_count = len(model.reactions)
    old_metabolite_count = len(model.metabolites)
    dummy_metabolite_1 = Metabolite("test_foo_1")
    dummy_metabolite_2 = Metabolite("test_foo_2")
    actual_metabolite = model.metabolites[0]
    copy_metabolite = model.metabolites[1].copy()
    dummy_reaction = Reaction("test_foo_reaction")
    dummy_reaction.add_metabolites({dummy_metabolite_1: -1,
                                    dummy_metabolite_2: 1,
                                    copy_metabolite: -2,
                                    actual_metabolite: 1})

    model.add_reaction(dummy_reaction)
    assert model.reactions.get_by_id(dummy_reaction.id) == dummy_reaction
    for x in [dummy_metabolite_1, dummy_metabolite_2]:
        assert model.metabolites.get_by_id(x.id) == x
    # Should add 1 reaction and 2 metabolites
    assert len(model.reactions) == old_reaction_count + 1
    assert len(model.metabolites) == old_metabolite_count + 2
    # Test the added reaction
    reaction_in_model = model.reactions.get_by_id(dummy_reaction.id)
    assert type(reaction_in_model) is Reaction
    assert reaction_in_model is dummy_reaction
    assert len(reaction_in_model._metabolites) == 4
    for i in reaction_in_model._metabolites:
        assert type(i) == Metabolite
    # Test the added metabolites
    met1_in_model = model.metabolites.get_by_id(dummy_metabolite_1.id)
    assert met1_in_model is dummy_metabolite_1
    copy_in_model = model.metabolites.get_by_id(copy_metabolite.id)
    assert copy_metabolite is not copy_in_model
    assert type(copy_in_model) is Metabolite
    assert dummy_reaction in actual_metabolite._reaction
    # Test addition of a different metabolite with the same name as an
    # existing one uses the metabolite in the model
    r2 = Reaction("test_foo_reaction2")
    model.add_reaction(r2)
    r2.add_metabolites({Metabolite(model.metabolites[0].id): 1})
    assert model.metabolites[0] is list(r2._metabolites)[0]


def test_add_reaction_context(model):
    old_reaction_count = len(model.reactions)
    old_metabolite_count = len(model.metabolites)
    dummy_metabolite_1 = Metabolite("test_foo_1")
    dummy_metabolite_2 = Metabolite("test_foo_2")
    actual_metabolite = model.metabolites[0]
    copy_metabolite = model.metabolites[1].copy()
    dummy_reaction = Reaction("test_foo_reaction")
    dummy_reaction.add_metabolites({dummy_metabolite_1: -1,
                                    dummy_metabolite_2: 1,
                                    copy_metabolite: -2,
                                    actual_metabolite: 1})
    dummy_reaction.gene_reaction_rule = 'dummy_gene'

    with model:
        model.add_reaction(dummy_reaction)
        assert model.reactions.get_by_id(
            dummy_reaction.id) == dummy_reaction
        assert len(model.reactions) == old_reaction_count + 1
        assert len(model.metabolites) == old_metabolite_count + 2
        assert dummy_metabolite_1._model == model
        assert 'dummy_gene' in model.genes

    assert len(model.reactions) == old_reaction_count
    assert len(model.metabolites) == old_metabolite_count
    with pytest.raises(KeyError):
        model.reactions.get_by_id(dummy_reaction.id)
    assert dummy_metabolite_1._model is None
    assert 'dummy_gene' not in model.genes


def test_add_reaction_from_other_model(model):
    other = model.copy()
    for i in other.reactions:
        i.id += "_other"
    other.repair()
    model.add_reactions(other.reactions)

    # Check if the other reaction has an error in its GPR
    m1 = model.copy()
    m2 = model.copy()
    m1.reactions.PGI.remove_from_model()
    m2.genes.b4025._reaction.clear()
    m1.add_reaction(m2.reactions.PGI)


def test_model_remove_reaction(model):
    old_reaction_count = len(model.reactions)

    with model:
        model.remove_reactions(["PGI"])
        assert len(model.reactions) == old_reaction_count - 1
        with pytest.raises(KeyError):
            model.reactions.get_by_id("PGI")
        model.remove_reactions(model.reactions[:1])
        assert len(model.reactions) == old_reaction_count - 2

    assert len(model.reactions) == old_reaction_count
    assert "PGI" in model.reactions

    tmp_metabolite = Metabolite("testing")
    model.reactions[0].add_metabolites({tmp_metabolite: 1})
    assert tmp_metabolite in model.metabolites
    model.remove_reactions(model.reactions[:1],
                           remove_orphans=True)
    assert tmp_metabolite not in model.metabolites

    with model:
        model.reactions[0].add_metabolites({tmp_metabolite: 1})
        assert tmp_metabolite in model.metabolites
    assert tmp_metabolite not in model.metabolites

    biomass_before = model.slim_optimize()
    with model:
        model.remove_reactions([model.reactions.Biomass_Ecoli_core])
        assert np.isclose(model.slim_optimize(), 0)

    assert np.isclose(model.slim_optimize(), biomass_before)


def test_reaction_remove(model):
    old_reaction_count = len(model.reactions)
    tmp_metabolite = Metabolite("testing")

    # Delete without removing orphan
    model.reactions[0].add_metabolites({tmp_metabolite: 1})
    assert len(tmp_metabolite.reactions) == 1

    # Esnsure the stoichiometry is still the same using different objects
    removed_reaction = model.reactions[0]
    original_stoich = {i.id: value for i, value
                       in removed_reaction._metabolites.items()}
    model.reactions[0].remove_from_model(remove_orphans=False)
    assert len(original_stoich) == len(removed_reaction._metabolites)
    for met in removed_reaction._metabolites:
        assert original_stoich[met.id] == removed_reaction._metabolites[
            met]
        assert met is not model.metabolites

    # Make sure it's still in the model
    assert tmp_metabolite in model.metabolites
    assert len(tmp_metabolite.reactions) == 0
    assert len(model.reactions) == old_reaction_count - 1

    # Now try with removing orphans
    model.reactions[0].add_metabolites({tmp_metabolite: 1})
    assert len(tmp_metabolite.reactions) == 1
    model.reactions[0].remove_from_model(remove_orphans=True)
    assert tmp_metabolite not in model.metabolites
    assert len(tmp_metabolite.reactions) == 0
    assert len(model.reactions) == old_reaction_count - 2

    # It shouldn't remove orphans if it's in 2 reactions however
    model.reactions[0].add_metabolites({tmp_metabolite: 1})
    model.reactions[1].add_metabolites({tmp_metabolite: 1})
    assert len(tmp_metabolite.reactions) == 2
    model.reactions[0].remove_from_model(remove_orphans=False)
    assert tmp_metabolite in model.metabolites
    assert len(tmp_metabolite.reactions) == 1
    assert len(model.reactions) == old_reaction_count - 3


def test_reaction_delete(model):
    old_reaction_count = len(model.reactions)
    tmp_metabolite = Metabolite("testing")

    # Delete without removing orphan
    model.reactions[0].add_metabolites({tmp_metabolite: 1})
    assert len(tmp_metabolite.reactions) == 1
    with pytest.warns(DeprecationWarning):
        model.reactions[0].delete(remove_orphans=False)

    # Make sure it's still in the model
    assert tmp_metabolite in model.metabolites
    assert len(tmp_metabolite.reactions) == 0
    assert len(model.reactions) == old_reaction_count - 1

    # Now try it with removing orphans
    model.reactions[0].add_metabolites({tmp_metabolite: 1})
    assert len(tmp_metabolite.reactions) == 1
    model.reactions[0].delete(remove_orphans=True)
    assert tmp_metabolite not in model.metabolites
    assert len(tmp_metabolite.reactions) == 0
    assert len(model.reactions) == old_reaction_count - 2

    # It shouldn't remove orphans if it's in 2 reactions however
    model.reactions[0].add_metabolites({tmp_metabolite: 1})
    model.reactions[1].add_metabolites({tmp_metabolite: 1})
    assert len(tmp_metabolite.reactions) == 2
    model.reactions[0].delete(remove_orphans=False)
    assert tmp_metabolite in model.metabolites
    assert len(tmp_metabolite.reactions) == 1
    assert len(model.reactions) == old_reaction_count - 3


def test_remove_gene(model):
    target_gene = model.genes[0]
    gene_reactions = list(target_gene.reactions)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        target_gene.remove_from_model()
    assert target_gene.model is None

    # Make sure the reaction was removed from the model
    assert target_gene not in model.genes

    # Ensure the old reactions no longer have a record of the gene
    for reaction in gene_reactions:
        assert target_gene not in reaction.genes


def test_exchange_reactions(model):
    assert set(model.exchanges) == set([rxn for rxn in model.reactions
                                        if rxn.id.startswith("EX")])


@pytest.mark.parametrize("metabolites, reaction_type, prefix", [
    ("exchange", "exchange", "EX_"),
    ("demand", "demand", "DM_"),
    ("sink", "sink", "SK_")
], indirect=["metabolites"])
def test_add_boundary(model, metabolites, reaction_type, prefix):
    for metabolite in metabolites:
        reaction = model.add_boundary(metabolite, reaction_type)
        assert model.reactions.get_by_id(
            reaction.id) == reaction
        assert reaction.reactants == [metabolite]
        assert model.constraints[metabolite.id].expression.has(
            model.variables[prefix + metabolite.id])


@pytest.mark.parametrize("metabolites, reaction_type, prefix", [
    ("exchange", "exchange", "EX_"),
    ("demand", "demand", "DM_"),
    ("sink", "sink", "SK_")
], indirect=["metabolites"])
def test_add_boundary_context(model, metabolites, reaction_type,
                              prefix):
    with model:
        for metabolite in metabolites:
            reaction = model.add_boundary(metabolite, reaction_type)
            assert model.reactions.get_by_id(
                reaction.id) == reaction
            assert reaction.reactants == [metabolite]
            assert -model.constraints[
                metabolite.id].expression.has(
                model.variables[prefix + metabolite.id])
    for metabolite in metabolites:
        assert prefix + metabolite.id not in model.reactions
        assert prefix + metabolite.id not in model.variables.keys()


@pytest.mark.parametrize("metabolites, reaction_type", [
    ("exchange", "exchange"),
    ("demand", "demand"),
    ("sink", "sink")
], indirect=["metabolites"])
def test_add_existing_boundary(model, metabolites, reaction_type):
    for metabolite in metabolites:
        model.add_boundary(metabolite, reaction_type)
        with pytest.raises(ValueError):
            model.add_boundary(metabolite, reaction_type)


@pytest.mark.parametrize("solver", stable_optlang)
def test_copy_benchmark(model, solver, benchmark):
    def _():
        model.copy()
        if not getattr(model, 'solver', None):
            solver_dict[solver].create_problem(model)

    benchmark(_)


@pytest.mark.parametrize("solver", stable_optlang)
def test_copy_benchmark_large_model(large_model, solver, benchmark):
    def _():
        large_model.copy()
        if not getattr(large_model, 'solver', None):
            solver_dict[solver].create_problem(large_model)

    benchmark(_)


def test_copy(model):
    # Modifying copy should not modify the original
    # Test that deleting reactions in the copy does not change the
    # number of reactions in the original model
    model_copy = model.copy()
    old_reaction_count = len(model.reactions)
    assert model_copy.notes is not model.notes
    assert model_copy.annotation is not model.annotation
    assert len(model.reactions) == len(model_copy.reactions)
    assert len(model.metabolites) == len(model_copy.metabolites)
    model_copy.remove_reactions(model_copy.reactions[0:5])
    assert old_reaction_count == len(model.reactions)
    assert len(model.reactions) != len(model_copy.reactions)
    # Copying a model should not copy its context
    with model:
        model.remove_reactions([model.reactions.ACALD])
        cp_model = model.copy()
        assert len(cp_model._contexts) == 0
    assert 'ACALD' not in cp_model.reactions


def test_deepcopy_benchmark(model, benchmark):
    benchmark(deepcopy, model)


def test_deepcopy(model):
    # Reference structures are maintained when deepcopying
    model_copy = deepcopy(model)
    for gene, gene_copy in zip(model.genes, model_copy.genes):
        assert gene.id == gene_copy.id
        reactions = sorted(i.id for i in gene.reactions)
        reactions_copy = sorted(i.id for i in gene_copy.reactions)
        assert reactions == reactions_copy
    for reaction, reaction_copy in zip(model.reactions,
                                       model_copy.reactions):
        assert reaction.id == reaction_copy.id
        metabolites = sorted(i.id for i in reaction._metabolites)
        metabolites_copy = sorted(i.id for i in reaction_copy._metabolites)
        assert metabolites == metabolites_copy


def test_add_reaction_orphans(model):
    # Test reaction addition
    # Need to verify that no orphan genes or metabolites are
    # contained in reactions after adding them to the model.
    model = model.__class__('test')
    model.add_reactions((x.copy() for x in model.reactions))
    genes = []
    metabolites = []
    for x in model.reactions:
        genes.extend(x.genes)
        metabolites.extend(x._metabolites)
    orphan_genes = [x for x in genes if x.model is not model]
    orphan_metabolites = [x for x in metabolites if x.model is not model]
    # Check for dangling genes when running Model.add_reactions
    assert len(orphan_genes) == 0
    # Check for dangling metabolites when running Model.add_reactions
    assert len(orphan_metabolites) == 0


def test_merge_models(model, tiny_toy_model):
    with model, tiny_toy_model:
        # Add some cons/vars to tiny_toy_model for testing merging
        tiny_toy_model.add_reactions([Reaction('EX_glc__D_e')])
        variable = tiny_toy_model.problem.Variable('foo')
        constraint = tiny_toy_model.problem.Constraint(
            variable, ub=0, lb=0, name='constraint')
        tiny_toy_model.add_cons_vars([variable, constraint])

        # Test merging to new model
        merged = model.merge(tiny_toy_model, inplace=False,
                             objective='sum', prefix_existing='tiny_')
        assert 'ex1' in merged.reactions
        assert 'ex1' not in model.reactions
        assert merged.reactions.ex1.objective_coefficient == 1
        assert merged.reactions.get_by_id(
            'Biomass_Ecoli_core').objective_coefficient == 1
        assert 'tiny_EX_glc__D_e' in merged.reactions
        assert 'foo' in merged.variables

        # Test reversible in-place model merging
        with model:
            model.merge(tiny_toy_model, inplace=True, objective='left',
                        prefix_existing='tiny_')
            assert 'ex1' in model.reactions
            assert 'constraint' in model.constraints
            assert 'foo' in model.variables
            assert 'tiny_EX_glc__D_e' in model.reactions
            assert (model.objective.expression.simplify() ==
                    model.reactions.get_by_id(
                        'Biomass_Ecoli_core').flux_expression.simplify())
        assert 'ex1' not in model.reactions
        assert 'constraint' not in model.constraints
        assert 'foo' not in model.variables
        assert 'tiny_EX_glc__D_e' not in model.reactions

    # Test the deprecated operator overloading
    with model:
        merged = model + tiny_toy_model
        assert 'ex1' in merged.reactions
    with model:
        model += tiny_toy_model
        assert 'ex1' in model.reactions


@pytest.mark.parametrize("solver", stable_optlang)
def test_change_objective_benchmark(model, benchmark, solver):
    atpm = model.reactions.get_by_id("ATPM")

    def benchmark_change_objective():
        model.objective = atpm.id
        if not getattr(model, 'solver', None):
            solver_dict[solver].create_problem(model)

    benchmark(benchmark_change_objective)


def test_get_objective_direction(model):
    assert model.objective_direction == "max"
    value = model.slim_optimize()
    assert np.isclose(value, 0.874, 1e-3)


def test_set_objective_direction(model):
    with model:
        model.objective_direction = "min"
        assert model.objective_direction == "min"
        value = model.slim_optimize()
        assert value == 0.0
    assert model.objective_direction == "max"


def test_slim_optimize(model):
    with model:
        assert model.slim_optimize() > 0.872
        model.reactions.Biomass_Ecoli_core.lower_bound = 10
        assert isnan(model.slim_optimize())
        with pytest.raises(OptimizationError):
            model.slim_optimize(error_value=None)


@pytest.mark.parametrize("solver", optlang_solvers)
def test_optimize(model, solver):
    model.solver = solver
    with model:
        assert model.optimize().objective_value > 0.872
        model.reactions.Biomass_Ecoli_core.lower_bound = 10
        with pytest.warns(UserWarning):
            model.optimize()
        with pytest.raises(OptimizationError):
            model.optimize(raise_error=True)


def test_change_objective(model):
    # Test for correct optimization behavior
    model.optimize()
    assert model.reactions.Biomass_Ecoli_core.x > 0.5
    with model:
        model.objective = model.reactions.EX_etoh_e
        model.optimize()
    assert model.reactions.Biomass_Ecoli_core.x < 0.5
    assert model.reactions.Biomass_Ecoli_core.objective_coefficient == 1
    model.optimize()
    assert model.reactions.Biomass_Ecoli_core.x > 0.5
    # Test changing objective
    biomass = model.reactions.get_by_id("Biomass_Ecoli_core")
    atpm = model.reactions.get_by_id("ATPM")
    model.objective = atpm.id
    assert atpm.objective_coefficient == 1.
    assert biomass.objective_coefficient == 0.
    assert su.linear_reaction_coefficients(model) == {atpm: 1.}
    # Change it back using object itself
    model.objective = biomass
    assert atpm.objective_coefficient == 0.
    assert biomass.objective_coefficient == 1.
    # Set both to 1 with a list
    model.objective = [atpm, biomass]
    assert atpm.objective_coefficient == 1.
    assert biomass.objective_coefficient == 1.
    # Set both using a dict
    model.objective = {atpm: 0.2, biomass: 0.3}
    assert abs(atpm.objective_coefficient - 0.2) < 10 ** -9
    assert abs(biomass.objective_coefficient - 0.3) < 10 ** -9
    # Test setting by index
    model.objective = model.reactions.index(atpm)
    assert su.linear_reaction_coefficients(model) == {atpm: 1.}
    # Test by setting list of indexes
    model.objective = [model.reactions.index(reaction) for
                       reaction in [atpm, biomass]]
    assert su.linear_reaction_coefficients(model) == {atpm: 1.,
                                                      biomass: 1.}


def test_problem_properties(model):
    new_variable = model.problem.Variable("test_variable")
    new_constraint = model.problem.Constraint(Zero,
                                              name="test_constraint")
    model.add_cons_vars([new_variable, new_constraint])
    assert "test_variable" in model.variables.keys()
    assert "test_constraint" in model.constraints.keys()
    model.remove_cons_vars([new_constraint, new_variable])
    assert "test_variable" not in model.variables.keys()
    assert "test_constraint" not in model.variables.keys()


def test_solution_data_frame(model):
    solution = model.optimize().to_frame()
    assert isinstance(solution, pd.DataFrame)
    assert 'fluxes' in solution
    assert 'reduced_costs' in solution


def test_context_manager(model):
    bounds0 = model.reactions[0].bounds
    bounds1 = (1, 2)
    bounds2 = (3, 4)

    # Trigger a nested model context, ensuring that bounds are
    # preserved at each level
    with model:
        model.reactions[0].bounds = bounds1
        with model:
            model.reactions[0].bounds = bounds2

            assert model.reactions[0].bounds == bounds2
        assert model.reactions[0].bounds == bounds1
    assert model.reactions[0].bounds == bounds0


def test_objective_coefficient_reflects_changed_objective(model):
    biomass_r = model.reactions.get_by_id('Biomass_Ecoli_core')
    assert biomass_r.objective_coefficient == 1
    model.objective = "PGI"
    assert biomass_r.objective_coefficient == 0
    assert model.reactions.PGI.objective_coefficient == 1


def test_change_objective_through_objective_coefficient(model):
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


def test_transfer_objective(model):
    new_mod = Model("new model")
    new_mod.add_reactions(model.reactions)
    new_mod.objective = model.objective
    assert (set(str(x) for x in model.objective.expression.args) == set(
        str(x) for x in new_mod.objective.expression.args))
    new_mod.slim_optimize()
    assert abs(new_mod.objective.value - 0.874) < 0.001


def test_model_from_other_model(model):
    model = Model(id_or_model=model)
    for reaction in model.reactions:
        assert reaction == model.reactions.get_by_id(reaction.id)


def test_add_reactions(model):
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


def test_add_reactions_single_existing(model):
    rxn = model.reactions[0]
    r1 = Reaction(rxn.id)
    r1.add_metabolites({Metabolite('A'): -1, Metabolite('B'): 1})
    r1.lower_bound, r1.upper_bound = -999999., 999999.
    model.add_reactions([r1])
    assert rxn in model.reactions
    assert r1 is not model.reactions.get_by_id(rxn.id)


def test_add_reactions_duplicate(model):
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


def test_add_cobra_reaction(model):
    r = Reaction(id="c1")
    model.add_reaction(r)
    assert isinstance(model.reactions.c1, Reaction)


def test_all_objects_point_to_all_other_correct_objects(model):
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


def test_objects_point_to_correct_other_after_copy(model):
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


def test_remove_reactions(model):
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


def test_objective(model):
    obj = model.objective
    assert obj.get_linear_coefficients(obj.variables) == {
               model.variables["Biomass_Ecoli_core_reverse_2cdba"]: -1,
               model.variables["Biomass_Ecoli_core"]: 1}
    assert obj.direction == "max"


def test_change_objective(model):
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


def test_set_reaction_objective(model):
    model.objective = model.reactions.ACALD
    assert str(model.objective.expression) == str(
        1.0 * model.reactions.ACALD.forward_variable -
        1.0 * model.reactions.ACALD.reverse_variable)


def test_set_reaction_objective_str(model):
    model.objective = model.reactions.ACALD.id
    assert str(model.objective.expression) == str(
        1.0 * model.reactions.ACALD.forward_variable -
        1.0 * model.reactions.ACALD.reverse_variable)


def test_invalid_objective_raises(model):
    with pytest.raises(ValueError):
        setattr(model, 'objective', 'This is not a valid objective!')
    with pytest.raises(TypeError):
        setattr(model, 'objective', 3.)


@pytest.mark.skipif("cplex" not in solvers, reason="need cplex")
def test_solver_change(model):
    solver_id = id(model.solver)
    problem_id = id(model.solver.problem)
    solution = model.optimize().fluxes
    model.solver = "cplex"
    assert id(model.solver) != solver_id
    assert id(model.solver.problem) != problem_id
    new_solution = model.optimize().fluxes
    assert np.allclose(solution, new_solution, rtol=0, atol=1E-06)


def test_no_change_for_same_solver(model):
    solver_id = id(model.solver)
    problem_id = id(model.solver.problem)
    model.solver = "glpk"
    assert id(model.solver) == solver_id
    assert id(model.solver.problem) == problem_id


def test_invalid_solver_change_raises(model):
    with pytest.raises(SolverNotFound):
        setattr(model, 'solver', [1, 2, 3])
    with pytest.raises(SolverNotFound):
        setattr(model, 'solver',
                'ThisIsDefinitelyNotAvalidSolver')
    with pytest.raises(SolverNotFound):
        setattr(model, 'solver', os)


@pytest.mark.skipif('cplex' not in solvers, reason='no cplex')
def test_change_solver_to_cplex_and_check_copy_works(model):
    assert round(abs(model.optimize().f - 0.8739215069684306), 7) == 0
    model_copy = model.copy()
    assert round(abs(model_copy.optimize().f - 0.8739215069684306),
                 7) == 0
    # Second, change existing glpk based model to cplex
    model.solver = 'cplex'
    assert round(abs(model.optimize().f - 0.8739215069684306),
                 7) == 0
    model_copy = copy(model)
    assert round(abs(model_copy.optimize().f - 0.8739215069684306),
                 7) == 0


def test_copy_preserves_existing_solution(solved_model):
    solution, model = solved_model
    model_cp = copy(model)
    primals_original = [variable.primal for variable in
                        model.variables]
    primals_copy = [variable.primal for variable in
                    model_cp.variables]
    abs_diff = abs(
        np.array(primals_copy) - np.array(primals_original))
    assert not any(abs_diff > 1e-6)


def test_repr_html_(model):
    assert '<table>' in model._repr_html_()
