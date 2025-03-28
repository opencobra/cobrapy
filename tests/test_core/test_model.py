"""Test functions of model.py."""

import os
import warnings
from copy import copy, deepcopy
from math import isnan
from typing import List, Tuple

import numpy as np
import pandas as pd
import pytest
from optlang.symbolics import Expr, Zero

from cobra import Solution
from cobra.core import Group, Metabolite, Model, Reaction
from cobra.exceptions import OptimizationError
from cobra.manipulation.delete import remove_genes
from cobra.util import solver as su
from cobra.util.solver import SolverNotFound, set_objective, solvers


try:
    import pytest_benchmark
except ImportError:
    pytest_benchmark = None

if pytest_benchmark:
    from pytest_benchmark.fixture import BenchmarkFixture


stable_optlang = ["glpk", "cplex", "gurobi"]
optlang_solvers = ["optlang-" + s for s in stable_optlang if s in su.solvers]


def same_ex(ex1: Expr, ex2: Expr) -> bool:
    """Compare two sympy-expressions for mathematical equality.

    Parameters
    ----------
    ex1 : optlang.symbolics.Expr
        The first sympy-expression.
    ex2 : optlang.symbolics.Expr
        The second sympy-expression.

    Returns
    -------
    bool
        Whether the two expression are mathematically equal.

    """
    return ex1.simplify() == ex2.simplify()


def test_add_metabolite(model: Model) -> None:
    """Tests adding a metabolite to a model, including with context.

    Parameters
    ----------
    model: cobra.Model
    """
    new_metabolite = Metabolite("test_met")
    assert new_metabolite not in model.metabolites
    with model:
        model.add_metabolites(new_metabolite)
        assert new_metabolite._model == model
        assert new_metabolite in model.metabolites
        assert new_metabolite.id in model.solver.constraints

    assert new_metabolite._model is None
    assert new_metabolite not in model.metabolites
    assert new_metabolite.id not in model.solver.constraints


def test_remove_metabolite_subtractive(model: Model) -> None:
    """Remove metabolite from model in a subtractive (not destructive) way.

    Checks that the changes to model are reversed when using context.

    Parameters
    ----------
    model: cobra.Model
    """
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


def test_remove_metabolite_destructive(model: Model) -> None:
    """Remove metabolite from a model in a destructive way.

    Checks that the changes to model are reversed when using context.

    Parameters
    ----------
    model: cobra.Model
    """
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


def test_compartments(model: Model) -> None:
    """Test setting and modifying model compartments.

    Parameters
    ----------
    model: cobra.Model
    """
    assert set(model.compartments) == {"c", "e"}
    model = Model("test", "test")
    met_c = Metabolite("a_c", compartment="c")
    met_e = Metabolite("a_e", compartment="e")
    rxn = Reaction("foo")
    rxn.add_metabolites({met_e: -1, met_c: 1})
    model.add_reactions([rxn])
    assert model.compartments == {"c": "", "e": ""}
    model.compartments = {"c": "cytosol"}
    assert model.compartments == {"c": "cytosol", "e": ""}


def test_model_remove_reaction(model: Model) -> None:
    """Test remove_reactions() to remove reaction(s).

    Parameters
    ----------
    model: cobra.Model

    """
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
    model.remove_reactions(model.reactions[:1], remove_orphans=True)
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


def test_reaction_remove(model: Model) -> None:
    """Test remove orphans in Reaction().remove_from_model.

    This function test that remove_orphans=True removes related metabolites when
    supposed to and doesn't remove when not supposed to (when this metabolite is in
    reactions not removed).

    Parameters
    ----------
    model: cobra.Model to use
    """
    old_reaction_count = len(model.reactions)
    tmp_metabolite = Metabolite("testing")

    # Delete without removing orphan
    model.reactions[0].add_metabolites({tmp_metabolite: 1})
    assert len(tmp_metabolite.reactions) == 1

    # Esnsure the stoichiometry is still the same using different objects
    removed_reaction = model.reactions[0]
    original_stoich = {
        i.id: value for i, value in removed_reaction._metabolites.items()
    }
    model.reactions[0].remove_from_model(remove_orphans=False)
    assert len(original_stoich) == len(removed_reaction._metabolites)
    for met in removed_reaction._metabolites:
        assert original_stoich[met.id] == removed_reaction._metabolites[met]
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


def test_reaction_delete(model: Model) -> None:
    """Test reaction removal using the Reaction.delete() function.

    This function calls Reaction().remove_from_model, since it is deprecated.

    Parameters
    ----------
    model: cobra.Model to use
    """
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


def test_remove_gene(model: Model) -> None:
    """Test remove_gene from model.

    Parameters
    ----------
    model: cobra.Model
    """
    target_gene = model.genes[0]
    gene_reactions = list(target_gene.reactions)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        remove_genes(model, [target_gene])
    assert target_gene.model is None

    # Make sure the reaction was removed from the model
    assert target_gene not in model.genes

    # Ensure the old reactions no longer have a record of the gene
    for reaction in gene_reactions:
        assert target_gene not in reaction.genes


def test_group_model_reaction_association(model: Model) -> None:
    """Test associating reactions with group in a model.

    This function will also remove the group from the model and check that
    reactions are no longer associated with the group

    Parameters
    ----------
    model: cobra.Model
    """
    num_members = 5
    reactions_for_group = model.reactions[0:num_members]
    group = Group("arbitrary_group1")
    group.add_members(reactions_for_group)
    group.kind = "collection"
    model.add_groups([group])
    # group should point to and be associated with the model
    assert group._model is model
    assert group in model.groups

    # model.get_associated_groups should find the group for each reaction
    # we added to the group
    for reaction in reactions_for_group:
        assert group in model.get_associated_groups(reaction)

    # remove the group from the model
    model.remove_groups([group])
    assert group not in model.groups
    assert group._model is not model
    for reaction in reactions_for_group:
        assert group not in model.get_associated_groups(reaction)


def test_group_members_add_to_model(model: Model) -> None:
    """Test adding a group with reactions to a model.

    This function will remove some reactions from a model, add them to a group, and
    test that the reactions aren't in the model. Later, the test will add the group
    to the model, and check the reactions were added to the model.

    Parameters
    ----------
    model: cobra.Model
    """
    # remove a few reactions from the model and add them to a new group
    num_members = 5
    reactions_for_group = model.reactions[0:num_members]
    model.remove_reactions(reactions_for_group, remove_orphans=False)
    group = Group("arbitrary_group1")
    group.add_members(reactions_for_group)
    group.kind = "collection"
    # the old reactions should not be in the model
    for reaction in reactions_for_group:
        assert reaction not in model.reactions

    # add the group to the model and check that the reactions were added
    model.add_groups([group])
    assert group in model.groups
    for reaction in reactions_for_group:
        assert reaction in model.reactions


def test_group_loss_of_elements(model: Model) -> None:
    """Test removal from model removes elements from group.

    This function will test that when a metabolite, reaction or gene is removed from
    a model, it no longer is a member of any groups.

    Parameters
    ----------
    model: cobra.Model
    """
    num_members_each = 5
    elements_for_group = model.reactions[0:num_members_each]
    elements_for_group.extend(model.metabolites[0:num_members_each])
    elements_for_group.extend(model.genes[0:num_members_each])
    group = Group("arbitrary_group1")
    group.add_members(elements_for_group)
    group.kind = "collection"
    model.add_groups([group])

    remove_met = model.metabolites[0]
    model.remove_metabolites([remove_met])
    remove_rxn = model.reactions[0]
    model.remove_reactions([remove_rxn])
    remove_gene = model.genes[0]
    remove_genes(model, [remove_gene])
    assert remove_met not in group.members
    assert remove_rxn not in group.members
    assert remove_gene not in group.members


def test_exchange_reactions(model: Model) -> None:
    """Test model.exchanges works as intended.

    Parameters
    ----------
    model: cobra.Model
    """
    assert set(model.exchanges) == {
        rxn for rxn in model.reactions if rxn.id.startswith("EX")
    }


@pytest.mark.parametrize(
    "metabolites, reaction_type, prefix",
    [
        ("exchange", "exchange", "EX_"),
        ("demand", "demand", "DM_"),
        ("sink", "sink", "SK_"),
    ],
    indirect=["metabolites"],
)
def test_add_boundary(
    model: Model, metabolites: List[Metabolite], reaction_type: str, prefix: str
) -> None:
    """Test add_boundary() function for model.

    Parameters
    ----------
    model: cobra.Model
    metabolites: List[Metabolites]
        This list is generated by the pytest.fixture metabolites(), see conftest.py
        in the root test directory.
    reaction_type: {"exchange", "demand", "sink"}
        The allowed types for boundary, see add_boundary() for types.
    prefix: str
    """
    for metabolite in metabolites:
        reaction = model.add_boundary(metabolite, reaction_type)
        assert model.reactions.get_by_id(reaction.id) == reaction
        assert reaction.reactants == [metabolite]
        assert model.constraints[metabolite.id].expression.has(
            model.variables[prefix + metabolite.id]
        )


@pytest.mark.parametrize(
    "metabolites, reaction_type, prefix",
    [
        ("exchange", "exchange", "EX_"),
        ("demand", "demand", "DM_"),
        ("sink", "sink", "SK_"),
    ],
    indirect=["metabolites"],
)
def test_add_boundary_context(
    model: Model, metabolites: List[Metabolite], reaction_type: str, prefix: str
) -> None:
    """Test add_boundary() function for model with context.

    Parameters
    ----------
    model: cobra.Model
    metabolites: List[Metabolites]
        This list is generated by the pytest.fixture metabolites(), see conftest.py
        in the root test directory.
    reaction_type: {"exchange", "demand", "sink"}
        The allowed types for boundary, see add_boundary() for types.
    prefix: str
    """
    with model:
        for metabolite in metabolites:
            reaction = model.add_boundary(metabolite, reaction_type)
            assert model.reactions.get_by_id(reaction.id) == reaction
            assert reaction.reactants == [metabolite]
            assert -model.constraints[metabolite.id].expression.has(
                model.variables[prefix + metabolite.id]
            )
    for metabolite in metabolites:
        assert prefix + metabolite.id not in model.reactions
        assert prefix + metabolite.id not in model.variables.keys()


@pytest.mark.parametrize(
    "metabolites, reaction_type",
    [("exchange", "exchange"), ("demand", "demand"), ("sink", "sink")],
    indirect=["metabolites"],
)
def test_add_existing_boundary(
    model: Model, metabolites: List[Metabolite], reaction_type: str
) -> None:
    """Test add_boundary() function for model with existing boundary/metabolite.

    Parameters
    ----------
    model: cobra.Model
    metabolites: List[Metabolites]
        This list is generated by the pytest.fixture metabolites(), see conftest.py
        in the root test directory.
    reaction_type: {"exchange", "demand", "sink"}
        The allowed types for boundary, see add_boundary() for types.
    """
    for metabolite in metabolites:
        rxn_added = model.add_boundary(metabolite, reaction_type)
        rxn_dup = model.add_boundary(metabolite, reaction_type)
        assert rxn_dup is rxn_added


@pytest.mark.parametrize("solver", optlang_solvers)
def test_copy_benchmark(model: Model, solver: str, benchmark: BenchmarkFixture) -> None:
    """Test copying a model with benchmark.

    Parameters
    ----------
    model: cobra.Model
    solver: str
        It is a string representing which solver to use. Parametized using
        'optlang_solvers' defined above.
    benchmark: BenchmarkFixture

    """

    def _() -> None:
        """Copy a model.

        If the model has no solver, it creates the problem using the solver given as
        parameter to the external function.
        """
        model.solver = solver
        model.copy()

    benchmark(_)


@pytest.mark.parametrize("solver", optlang_solvers)
def test_copy_benchmark_large_model(
    large_model: Model,
    solver: str,
    benchmark: BenchmarkFixture,
) -> None:
    """Test copying a large model with benchmark.

    Parameters
    ----------
    large_model: cobra.Model
    solver: str
        It is a string representing which solver to use. Parametized using
        'optlang_solvers' defined above.
    benchmark: BenchmarkFixture
    """

    def _() -> None:
        """Copy a large model.

        If the model has no solver, it creates the problem using the solver given as
        parameter to the external function.
        """
        large_model.solver = solver
        large_model.copy()

    benchmark(_)


def test_copy(model: Model) -> None:
    """Test copying a model and modifying the copy.

    This function tests that modifying the copy should not modifying the original, by
    deleting reactions in the copy (# of reactions in the original should not change).
    This function also tests that GPRs are copied by content, not by reference, and
    that the model copy does not copy the context.

    Parameters
    ----------
    model: cobra.Model

    """
    # Deleting reactions in copy does not change number of reactions in the original
    model_copy = model.copy()
    old_reaction_count = len(model.reactions)
    assert model_copy.notes is not model.notes
    assert model_copy.annotation is not model.annotation
    assert len(model.reactions) == len(model_copy.reactions)
    assert len(model.metabolites) == len(model_copy.metabolites)
    assert len(model.groups) == len(model_copy.groups)
    assert len(model.genes) == len(model_copy.genes)
    # test if GPRs are copied by content but not by reference
    assert model.reactions[0].gpr == model_copy.reactions[0].gpr
    assert id(model.reactions[0].gpr.body) != id(model_copy.reactions[0].gpr.body)
    model_copy.remove_reactions(model_copy.reactions[0:5])
    assert old_reaction_count == len(model.reactions)
    assert len(model.reactions) != len(model_copy.reactions)
    # Copying a model should not copy its context
    with model:
        model.remove_reactions([model.reactions.ACALD])
        cp_model = model.copy()
        assert len(cp_model._contexts) == 0
    assert "ACALD" not in cp_model.reactions


def test_copy_with_groups(model: Model) -> None:
    """Copy model with groups and check that groups are copied correctly.

    Parameters
    ----------
    model: cobra.Model

    """
    sub = Group("pathway", members=[model.reactions.PFK, model.reactions.FBA])
    model.add_groups([sub])
    model_copy = model.copy()
    assert len(model_copy.groups) == len(model.groups)
    assert len(model_copy.groups.get_by_id("pathway")) == len(
        model.groups.get_by_id("pathway")
    )


def test_deepcopy_benchmark(model: Model, benchmark: BenchmarkFixture) -> None:
    """Benchmark deepcopying a model.

    Parameters
    ----------
    model: cobra.Model
    benchmark: BenchmarkFixture
    """
    benchmark(deepcopy, model)


def test_deepcopy(model: Model) -> None:
    """Test deepcopying works, and maintains reference structures.

    Parameters
    ----------
    model: cobra.Model
    """
    # Reference structures are maintained when deepcopying
    model_copy = deepcopy(model)
    for gene, gene_copy in zip(model.genes, model_copy.genes):
        assert gene.id == gene_copy.id
        reactions = sorted(i.id for i in gene.reactions)
        reactions_copy = sorted(i.id for i in gene_copy.reactions)
        assert reactions == reactions_copy
    for reaction, reaction_copy in zip(model.reactions, model_copy.reactions):
        assert reaction.id == reaction_copy.id
        metabolites = sorted(i.id for i in reaction._metabolites)
        metabolites_copy = sorted(i.id for i in reaction_copy._metabolites)
        assert metabolites == metabolites_copy


def test_add_reaction_orphans(model: Model) -> None:
    """Test orphan behavior when adding reactions.

    Need to verify that no orphan genes or metabolites are contained in reactions
    after adding them to the model.

    Parameters
    ---------
    model: cobra.Model
    """
    model = model.__class__("test")
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


def test_merge_models(model: Model, tiny_toy_model: Model) -> None:
    """Test merging models.

    Parameters
    ----------
    model: cobra.Model
    tiny_toy_model: cobra.Model
    """
    with model, tiny_toy_model:
        # Add some cons/vars to tiny_toy_model for testing merging
        tiny_toy_model.add_reactions([Reaction("EX_glc__D_e")])
        variable = tiny_toy_model.problem.Variable("foo")
        constraint = tiny_toy_model.problem.Constraint(
            variable, ub=0, lb=0, name="constraint"
        )
        tiny_toy_model.add_cons_vars([variable, constraint])

        merged = model.merge(
            tiny_toy_model, inplace=False, objective="sum", prefix_existing="tiny_"
        )
        assert "ex1" in merged.reactions
        assert "ex1" not in model.reactions
        assert merged.reactions.ex1.objective_coefficient == 1
        assert (
            merged.reactions.get_by_id("Biomass_Ecoli_core").objective_coefficient == 1
        )
        assert "tiny_EX_glc__D_e" in merged.reactions
        assert "foo" in merged.variables

        # Test reversible in-place model merging
        with model:
            model.merge(
                tiny_toy_model, inplace=True, objective="left", prefix_existing="tiny_"
            )
            assert "ex1" in model.reactions
            assert "constraint" in model.constraints
            assert "foo" in model.variables
            assert "tiny_EX_glc__D_e" in model.reactions
            assert (
                model.objective.expression.simplify()
                == model.reactions.get_by_id(
                    "Biomass_Ecoli_core"
                ).flux_expression.simplify()
            )
        assert "ex1" not in model.reactions
        assert "constraint" not in model.constraints
        assert "foo" not in model.variables
        assert "tiny_EX_glc__D_e" not in model.reactions


@pytest.mark.parametrize("solver", optlang_solvers)
def test_change_objective_benchmark(
    model: Model, benchmark: BenchmarkFixture, solver: str
) -> None:
    """Benchmark changing objective in model.

    Parameters
    ----------
    model: cobra.Model
    benchmark: BenchmarkFixture
    solver: str
        Solver to use. Parametized using 'optlang_solvers' defined above.
    """
    atpm = model.reactions.get_by_id("ATPM")

    def benchmark_change_objective():
        model.objective = atpm.id
        model.solver = solver

    benchmark(benchmark_change_objective)


def test_get_objective_direction(model: Model) -> None:
    """Test getting objective.

    Parameters
    ----------
    model: cobra.Model
    """
    assert model.objective_direction == "max"
    value = model.slim_optimize()
    assert np.isclose(value, 0.874, 1e-3)


def test_set_objective_direction(model: Model) -> None:
    """Test setting objective.

    Parameters
    ----------
    model: cobra.Model
    """
    with model:
        model.objective_direction = "min"
        assert model.objective_direction == "min"
        value = model.slim_optimize()
        assert value == 0.0
    assert model.objective_direction == "max"


def test_slim_optimize(model: Model) -> None:
    """Test slim_optimize with context.

    Parameters
    ----------
    model: cobra.Model

    """
    with model:
        assert model.slim_optimize() > 0.872
        model.reactions.Biomass_Ecoli_core.lower_bound = 10
        assert isnan(model.slim_optimize())
        with pytest.raises(OptimizationError):
            model.slim_optimize(error_value=None)


@pytest.mark.parametrize("solver", optlang_solvers)
def test_optimize(model: Model, solver: str) -> None:
    """Test optimizing a model.

    Parameters
    ----------
    model: cobra.Model
    solver: str
        Solver to use. Parametized using 'optlang_solvers' defined above.
    """
    model.solver = solver
    with model:
        assert model.optimize().objective_value > 0.872
        model.reactions.Biomass_Ecoli_core.lower_bound = 10
        with pytest.warns(UserWarning):
            model.optimize()
        with pytest.raises(OptimizationError):
            model.optimize(raise_error=True)


def test_change_objective(model: Model) -> None:
    """Test changing objective.

    Parameters
    ----------
    model: cobra.Model
    """
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
    assert atpm.objective_coefficient == 1.0
    assert biomass.objective_coefficient == 0.0
    assert su.linear_reaction_coefficients(model) == {atpm: 1.0}
    # Change it back using object itself
    model.objective = biomass
    assert atpm.objective_coefficient == 0.0
    assert biomass.objective_coefficient == 1.0
    # Set both to 1 with a list
    model.objective = [atpm, biomass]
    assert atpm.objective_coefficient == 1.0
    assert biomass.objective_coefficient == 1.0
    # Set both using a dict
    model.objective = {atpm: 0.2, biomass: 0.3}
    assert abs(atpm.objective_coefficient - 0.2) < 10**-9
    assert abs(biomass.objective_coefficient - 0.3) < 10**-9
    # Test setting by index
    model.objective = model.reactions.index(atpm)
    assert su.linear_reaction_coefficients(model) == {atpm: 1.0}
    # Test by setting list of indexes
    model.objective = [model.reactions.index(reaction) for reaction in [atpm, biomass]]
    assert su.linear_reaction_coefficients(model) == {atpm: 1.0, biomass: 1.0}


def test_problem_properties(model: Model) -> None:
    """Test model problem properties.

    Parameters
    ----------
    model: cobra.Model
    """
    new_variable = model.problem.Variable("test_variable")
    new_constraint = model.problem.Constraint(Zero, name="test_constraint", lb=0)
    model.add_cons_vars([new_variable, new_constraint])
    assert "test_variable" in model.variables
    assert "test_constraint" in model.constraints
    model.remove_cons_vars([new_constraint, new_variable])
    assert "test_variable" not in model.variables
    assert "test_constraint" not in model.variables


def test_solution_data_frame(model: Model) -> None:
    """Test that solution is transformed correctly to a Pandas data frame.

    Parameters
    ----------
    model: cobra. Model

    """
    solution = model.optimize().to_frame()
    assert isinstance(solution, pd.DataFrame)
    assert "fluxes" in solution
    assert "reduced_costs" in solution


def test_context_manager(model: Model) -> None:
    """Test that the context manager works.

    Parameters
    ----------
    model: cobra.Model

    """
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


def test_objective_coefficient_reflects_changed_objective(model: Model) -> None:
    """Test that changing objectives is reflected in the objectives changing.

    Parameters
    ----------
    model: cobra.Model
    """
    biomass_r = model.reactions.get_by_id("Biomass_Ecoli_core")
    assert biomass_r.objective_coefficient == 1
    model.objective = "PGI"
    assert biomass_r.objective_coefficient == 0
    assert model.reactions.PGI.objective_coefficient == 1


def test_change_objective_through_objective_coefficient(model: Model) -> None:
    """Test that changing the objective coefficients will change the objective.

    Parameters
    ----------
    model: cobra.Model
    """
    biomass_r = model.reactions.get_by_id("Biomass_Ecoli_core")
    pgi = model.reactions.PGI
    pgi.objective_coefficient = 2
    coef_dict = model.objective.expression.as_coefficients_dict()
    # Check that objective has been updated
    assert coef_dict[pgi.forward_variable] == 2.0
    assert coef_dict[pgi.reverse_variable] == -2.0
    # Check that original objective is still in there
    assert coef_dict[biomass_r.forward_variable] == 1.0
    assert coef_dict[biomass_r.reverse_variable] == -1.0


def test_transfer_objective(model: Model) -> None:
    """Test assigning objective from a different mdoel objective.

    Parameters
    ----------
    model: cobra.Model
    """
    new_mod = Model("new model")
    new_mod.add_reactions(model.reactions)
    new_mod.objective = model.objective
    assert {str(x) for x in model.objective.expression.args} == {
        str(x) for x in new_mod.objective.expression.args
    }
    new_mod.slim_optimize()
    assert abs(new_mod.objective.value - 0.874) < 0.001


def test_model_from_other_model(model: Model) -> None:
    """Test creating model from other model.

    Parameters
    ----------
    model: cobra.Model
    """
    model = Model(id_or_model=model)
    for reaction in model.reactions:
        assert reaction == model.reactions.get_by_id(reaction.id)


def test_add_reactions(model: Model) -> None:
    """Test add_reactions() function to add reactions to model.

    Parameters
    ----------
    model: cobra.Model
    """
    r1 = Reaction("r1")
    r1.add_metabolites({Metabolite("A"): -1, Metabolite("B"): 1})
    r1.lower_bound, r1.upper_bound = -999999.0, 999999.0
    r2 = Reaction("r2")
    r2.add_metabolites({Metabolite("A"): -1, Metabolite("C"): 1, Metabolite("D"): 1})
    r2.lower_bound, r2.upper_bound = 0.0, 999999.0
    model.add_reactions([r1, r2])
    r2.objective_coefficient = 3.0
    assert r2.objective_coefficient == 3.0
    assert model.reactions[-2] == r1
    assert model.reactions[-1] == r2
    assert isinstance(model.reactions[-2].reverse_variable, model.problem.Variable)
    coefficients_dict = model.objective.expression.as_coefficients_dict()
    biomass_r = model.reactions.get_by_id("Biomass_Ecoli_core")
    assert coefficients_dict[biomass_r.forward_variable] == 1.0
    assert coefficients_dict[biomass_r.reverse_variable] == -1.0
    assert coefficients_dict[model.reactions.r2.forward_variable] == 3.0
    assert coefficients_dict[model.reactions.r2.reverse_variable] == -3.0


def test_add_reactions_single_existing(model: Model) -> None:
    """Test adding a reaction already present to a model.

    Parameters
    ----------
    model: cobra.Model
    """
    rxn = model.reactions[0]
    r1 = Reaction(rxn.id)
    r1.add_metabolites({Metabolite("A"): -1, Metabolite("B"): 1})
    r1.lower_bound, r1.upper_bound = -999999.0, 999999.0
    model.add_reactions([r1])
    assert rxn in model.reactions
    assert r1 is not model.reactions.get_by_id(rxn.id)


def test_add_reactions_duplicate(model: Model) -> None:
    """Test adding duplicate reactions to a model.

    Parameters
    ----------
    model: cobra.Model
    """
    rxn = model.reactions[0]
    r1 = Reaction("r1")
    r1.add_metabolites({Metabolite("A"): -1, Metabolite("B"): 1})
    r1.lower_bound, r1.upper_bound = -999999.0, 999999.0
    r2 = Reaction(rxn.id)
    r2.add_metabolites({Metabolite("A"): -1, Metabolite("C"): 1, Metabolite("D"): 1})
    model.add_reactions([r1, r2])
    assert r1 in model.reactions
    assert rxn in model.reactions
    assert r2 is not model.reactions.get_by_id(rxn.id)


def test_all_objects_point_to_all_other_correct_objects(model: Model) -> None:
    """Test that objects point to needed other objects.

    This will test that the reaction.genes, reaction.metabolites point to the correct
    genes and metabolites in the model.

    Parameters
    ----------
    model: cobra.Model
    """
    for reaction in model.reactions:
        assert reaction.model == model
        for gene in reaction.genes:
            assert gene == model.genes.get_by_id(gene.id)
            assert gene.model == model
            for reaction2 in gene.reactions:
                assert reaction2.model == model
                assert reaction2 == model.reactions.get_by_id(reaction2.id)

        for metabolite in reaction.metabolites:
            assert metabolite.model == model
            assert metabolite == model.metabolites.get_by_id(metabolite.id)
            for reaction2 in metabolite.reactions:
                assert reaction2.model == model
                assert reaction2 == model.reactions.get_by_id(reaction2.id)


def test_objects_point_to_correct_other_after_copy(model: Model) -> None:
    """Test that objects point to correct other objects after copying a model.

    Parameters
    ----------
    model: cobra.Model
    """
    for reaction in model.reactions:
        assert reaction.model == model
        for gene in reaction.genes:
            assert gene == model.genes.get_by_id(gene.id)
            assert gene.model == model
            for reaction2 in gene.reactions:
                assert reaction2.model == model
                assert reaction2 == model.reactions.get_by_id(reaction2.id)

        for metabolite in reaction.metabolites:
            assert metabolite.model == model
            assert metabolite == model.metabolites.get_by_id(metabolite.id)
            for reaction2 in metabolite.reactions:
                assert reaction2.model == model
                assert reaction2 == model.reactions.get_by_id(reaction2.id)


def test_remove_reactions(model: Model) -> None:
    """Test remove_reactions() from Model.

    Parameters
    ----------
    model: cobra.Model
    """
    reactions_to_remove = model.reactions[10:30]
    assert all([reaction.model is model for reaction in reactions_to_remove])
    assert all(
        [
            model.reactions.get_by_id(reaction.id) == reaction
            for reaction in reactions_to_remove
        ]
    )

    model.remove_reactions(reactions_to_remove)
    assert all([reaction.model is None for reaction in reactions_to_remove])
    for reaction in reactions_to_remove:
        assert reaction.id not in list(model.variables.keys())

    model.add_reactions(reactions_to_remove)
    for reaction in reactions_to_remove:
        assert reaction in model.reactions


def test_objective(model: Model) -> None:
    """Test that objective contains the correct coefficients.

    Parameters
    ----------
    model: cobra.Model
    """
    obj = model.objective
    assert obj.get_linear_coefficients(obj.variables) == {
        model.variables["Biomass_Ecoli_core_reverse_2cdba"]: -1,
        model.variables["Biomass_Ecoli_core"]: 1,
    }
    assert obj.direction == "max"


def test_change_objective_with_context(model: Model) -> None:
    """Test changing objective is reversed with context.

    Parameters
    ----------
    model: cobra.Model
    """
    expression = 1.0 * model.variables["ENO"] + 1.0 * model.variables["PFK"]
    model.objective = model.problem.Objective(expression)
    assert same_ex(model.objective.expression, expression)
    model.objective = "ENO"
    eno_obj = model.problem.Objective(
        model.reactions.ENO.flux_expression, direction="max"
    )
    pfk_obj = model.problem.Objective(
        model.reactions.PFK.flux_expression, direction="max"
    )
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
        set_objective(model, model.problem.Objective(atpm.flux_expression))
        assert same_ex(model.objective.expression, atpm.flux_expression)
    assert same_ex(model.objective.expression, expression)

    expression = model.objective.expression
    with model:
        with model:  # Test to make sure nested contexts are OK
            set_objective(model, atpm.flux_expression, additive=True)
            assert same_ex(
                model.objective.expression, expression + atpm.flux_expression
            )
    assert same_ex(model.objective.expression, expression)


def test_set_reaction_objective(model: Model) -> None:
    """Test setting reaction objective.

    Parameters
    ----------
    model: cobra.Model
    """
    model.objective = model.reactions.ACALD
    assert same_ex(
        model.objective.expression,
        1.0 * model.reactions.ACALD.forward_variable
        - 1.0 * model.reactions.ACALD.reverse_variable,
    )


def test_set_reaction_objective_str(model: Model) -> None:
    """Test setting reaction objective using string.

    Parameters
    ----------
    model: cobra.Model
    """
    model.objective = model.reactions.ACALD.id
    assert same_ex(
        model.objective.expression,
        1.0 * model.reactions.ACALD.forward_variable
        - 1.0 * model.reactions.ACALD.reverse_variable,
    )


def test_invalid_objective_raises(model: Model) -> None:
    """Test that an invalid objective will raise appropriate errors.

    Parameters
    ----------
    model: cobra.Model
    """
    with pytest.raises(ValueError):
        model.objective = "This is not a valid objective!"
    with pytest.raises(TypeError):
        model.objective = 3.0


@pytest.mark.skipif("cplex" not in solvers, reason="need cplex")
def test_solver_change(model: Model) -> None:
    """Test changing the solver to cplex.

    This test will be skipped if cplex is not installed and available to python.

    Parameters
    ----------
    model: cobra.Model
    """
    model.solver = "glpk"
    solver_id = id(model.solver)
    problem_id = id(model.solver.problem)
    solution = model.optimize().fluxes
    model.solver = "cplex"
    assert id(model.solver) != solver_id
    assert id(model.solver.problem) != problem_id
    new_solution = model.optimize().fluxes
    assert np.allclose(solution, new_solution, rtol=0, atol=1e-06)


def test_no_change_for_same_solver(model: Model) -> None:
    """Test no change in variables if changing to the same solver as before.

    Parameters
    ----------
    model: cobra.Model
    """
    model.solver = "glpk"
    solver_id = id(model.solver)
    problem_id = id(model.solver.problem)
    model.solver = "glpk"
    assert id(model.solver) == solver_id
    assert id(model.solver.problem) == problem_id


def test_invalid_solver_change_raises(model: Model) -> None:
    """Test changing to an invalid solver will raise SovlerNotFound.

    Parameters
    ----------
    model: cobra.Model
    """
    with pytest.raises(SolverNotFound):
        model.solver = [1, 2, 3]
    with pytest.raises(SolverNotFound):
        model.solver = "ThisIsDefinitelyNotAvalidSolver"
    with pytest.raises(SolverNotFound):
        model.solver = os


@pytest.mark.skipif("cplex" not in solvers, reason="no cplex")
def test_change_solver_to_cplex_and_check_copy_works(model: Model) -> None:
    """Test changing solver and copying model work.

    This test will be skipped if cplex is not installed and available for python.

    Parameters
    ----------
    model: cobra.Model
    """
    assert (model.slim_optimize() - 0.8739215069684306) == pytest.approx(0.0)
    model_copy = model.copy()
    assert (model_copy.slim_optimize() - 0.8739215069684306) == pytest.approx(0.0)
    # Second, change existing glpk based model to cplex
    model.solver = "cplex"
    assert (model.slim_optimize() - 0.8739215069684306) == pytest.approx(0.0)
    model_copy = copy(model)
    assert (model_copy.slim_optimize() - 0.8739215069684306) == pytest.approx(0.0)


def test_copy_preserves_existing_solution(solved_model: Tuple[Solution, Model]) -> None:
    """Test copy keeps the existing solution.

    Primal values are the same when copying a solved model.

    Parameters
    ----------
    solved_model: Tuple
        A Tuple that contains a Solution and a Model
    """
    solution, model = solved_model
    model_cp = copy(model)
    primals_original = [variable.primal for variable in model.variables]
    primals_copy = [variable.primal for variable in model_cp.variables]
    abs_diff = abs(np.array(primals_copy) - np.array(primals_original))
    assert not any(abs_diff > 1e-6)


def test_repr_html_(model: Model) -> None:
    """Test HTML representation of model.

    Parameters
    ----------
    model: cobra.Model
    """
    assert "<table>" in model._repr_html_()
