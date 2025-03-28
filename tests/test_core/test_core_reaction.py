"""Test functions of cobra.core.reaction ."""

import warnings
from typing import Iterable

import numpy as np
import pytest

from cobra import Gene
from cobra.core import GPR, Configuration, Metabolite, Model, Reaction


config = Configuration()
stable_optlang = ["glpk", "cplex", "gurobi"]


def test_gpr() -> None:
    """Test GPR evaluation."""
    model = Model()
    reaction = Reaction("test")

    # Set GPR to an empty string
    reaction.gene_reaction_rule = ""
    # Empty gene_reaction_rule leads to an empty GPR
    assert reaction.gpr.body is None
    assert reaction.gpr.to_string() == ""
    assert reaction.gpr.to_string(names={"goo": "blah"}) == ""
    # Set GPR directly (shouldn't really be used, but just a test)
    reaction.gpr = GPR()
    assert reaction.gene_reaction_rule == ""
    assert reaction.gpr.eval()
    # Set GPR to a reaction not in a model
    reaction.gene_reaction_rule = "(g1 or g2) and g3"
    assert reaction.gene_reaction_rule == "(g1 or g2) and g3"
    assert len(reaction.genes) == 3

    # Adding reaction with a GPR propagates to the model
    model.add_reactions([reaction])
    assert len(model.genes) == 3

    # Ensure the gene objects are the same in the model and reaction
    reaction_gene = list(reaction.genes)[0]
    model_gene = model.genes.get_by_id(reaction_gene.id)
    assert reaction_gene is model_gene


def test_gpr_uppercase() -> None:
    """Test ability to handle uppercase AND/OR."""
    reaction = Reaction("test")
    with pytest.warns(SyntaxWarning):
        reaction.gene_reaction_rule = "(b1 AND b2) OR (b3 and b4)"
        assert reaction.gene_reaction_rule == "(b1 and b2) or (b3 and b4)"
        assert len(reaction.genes) == 4


@pytest.mark.parametrize("input_gpr", ["(a1 or a2", "(forT or "])
def test_gpr_malformed(input_gpr: str) -> None:
    """Test ability to deal with malformed GPR.

    Malformed GPR strings will lead to empty GPRs with no genes.

    Parameters
    ----------
    input_gpr: str
        String representing a malformed GPR string.
    """
    reaction = Reaction("test")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        reaction.gene_reaction_rule = input_gpr
        assert len(reaction.genes) == 0


def test_gpr_modification(model: Model) -> None:
    """Test GPR manipulations."""
    reaction = model.reactions.get_by_id("PGI")
    old_gene = list(reaction.genes)[0]
    new_gene = model.genes.get_by_id("s0001")

    # Add an existing 'gene' to the GPR
    reaction.gene_reaction_rule = "s0001"
    assert new_gene in reaction.genes
    assert reaction in new_gene.reactions

    # Remove old gene correctly, keep it in model
    assert old_gene not in reaction.genes
    assert reaction not in old_gene.reactions
    assert old_gene in model.genes

    # Add a new 'gene' to the GPR
    reaction.gene_reaction_rule = "fake_gene"
    assert model.genes.has_id("fake_gene")
    fake_gene = model.genes.get_by_id("fake_gene")
    assert fake_gene in reaction.genes
    assert reaction in fake_gene.reactions
    fake_gene.name = "foo_gene"
    assert reaction.gene_name_reaction_rule == fake_gene.name


def test_gpr_modification_with_context(model: Model) -> None:
    """Test GPR manipulations are reversed in context."""
    empty_model = Model()
    reaction = Reaction("test")
    reaction.gene_reaction_rule = "(g1 or g2) and g3"
    assert reaction.gene_reaction_rule == "(g1 or g2) and g3"
    assert len(reaction.genes) == 3

    with empty_model:
        # Adding reaction with a GPR propagates to the model
        empty_model.add_reactions([reaction])
        assert len(empty_model.genes) == 3
    assert len(empty_model.reactions) == 0
    assert len(empty_model.genes) == 0
    assert reaction._model is None

    reaction = model.reactions.get_by_id("PGI")
    old_reaction_rule = reaction.gene_reaction_rule
    old_gene = list(reaction.genes)[0]
    new_gene = model.genes.get_by_id("s0001")

    with model:
        # Add an existing 'gene' to the GPR
        reaction.gene_reaction_rule = "s0001"
        assert new_gene in reaction.genes
        assert reaction in new_gene.reactions

        # Remove old gene correctly, keep it in model
        assert old_gene not in reaction.genes
        assert reaction not in old_gene.reactions
        assert old_gene in model.genes

    assert reaction.gene_reaction_rule == old_reaction_rule
    assert new_gene not in reaction.genes
    assert reaction not in new_gene.reactions

    # Remove old gene correctly, keep it in model
    assert old_gene in reaction.genes
    assert reaction in old_gene.reactions
    assert old_gene in model.genes

    with model:
        # Add a new 'gene' to the GPR
        reaction.gene_reaction_rule = "fake_gene"
        assert model.genes.has_id("fake_gene")
        fake_gene = model.genes.get_by_id("fake_gene")
        assert fake_gene in reaction.genes
        assert reaction in fake_gene.reactions
        fake_gene.name = "foo_gene"
        assert reaction.gene_name_reaction_rule == fake_gene.name
    assert not model.genes.has_id("fake_gene")
    fake_gene = Gene("fake_gene")
    assert fake_gene not in reaction.genes
    assert reaction not in fake_gene.reactions


def test_gene_knock_out(model: Model) -> None:
    """Test gene knockout effect on reaction."""
    rxn = Reaction("rxn")
    rxn.add_metabolites({Metabolite("A"): -1, Metabolite("B"): 1})
    rxn.gene_reaction_rule = "A2B1 or A2B2 and A2B3"
    assert hasattr(list(rxn.genes)[0], "knock_out")
    model.add_reactions([rxn])
    with model:
        model.genes.A2B1.knock_out()
        assert not model.genes.A2B1.functional
        model.genes.A2B3.knock_out()
        assert not rxn.functional
    assert model.genes.A2B3.functional
    assert rxn.functional
    model.genes.A2B1.knock_out()
    assert not model.genes.A2B1.functional
    assert model.reactions.rxn.functional
    model.genes.A2B3.knock_out()
    assert not model.reactions.rxn.functional


def test_str() -> None:
    """Test `str` output for a reaction."""
    rxn = Reaction("rxn")
    rxn.add_metabolites({Metabolite("A"): -1, Metabolite("B"): 1})
    assert str(rxn) == "rxn: A --> B"


def test_str_from_model(model: Model) -> None:
    """Test `str` output for a reaction associated with a model."""
    assert model.reactions[0].__str__().startswith("ACALD")


def test_add_metabolite_from_solved_model(solved_model: Model) -> None:
    """Test metabolite addition to a reaction from a solved model."""
    solution, model = solved_model
    pgi_reaction = model.reactions.PGI
    test_met = model.metabolites[0]
    pgi_reaction.add_metabolites({test_met: 42}, combine=False)
    assert pgi_reaction.metabolites[test_met] == 42.0
    assert (
        model.constraints[test_met.id].expression.as_coefficients_dict()[
            pgi_reaction.forward_variable
        ]
        == 42.0
    )
    assert (
        model.constraints[test_met.id].expression.as_coefficients_dict()[
            pgi_reaction.reverse_variable
        ]
        == -42.0
    )

    pgi_reaction.add_metabolites({test_met: -10}, combine=True)
    assert pgi_reaction.metabolites[test_met] == 32.0
    assert (
        model.constraints[test_met.id].expression.as_coefficients_dict()[
            pgi_reaction.forward_variable
        ]
        == 32.0
    )
    assert (
        model.constraints[test_met.id].expression.as_coefficients_dict()[
            pgi_reaction.reverse_variable
        ]
        == -32.0
    )

    pgi_reaction.add_metabolites({test_met: 0}, combine=False)
    with pytest.raises(KeyError):
        pgi_reaction.metabolites[test_met]
    assert (
        model.constraints[test_met.id].expression.as_coefficients_dict()[
            pgi_reaction.forward_variable
        ]
        == 0
    )
    assert (
        model.constraints[test_met.id].expression.as_coefficients_dict()[
            pgi_reaction.reverse_variable
        ]
        == 0
    )


@pytest.mark.parametrize("solver", stable_optlang)
def test_add_metabolite_benchmark(model: Model, benchmark, solver: Iterable) -> None:
    """Benchmark metabolite addition to a reaction associated with a model."""
    reaction = model.reactions.get_by_id("PGI")
    many_metabolites = dict((m, 1) for m in model.metabolites[0:50])

    def add_remove_metabolite():
        reaction.add_metabolites(many_metabolites)
        if not getattr(model, "solver", None):
            stable_optlang[solver].create_problem(model)
        for met in many_metabolites:
            try:
                reaction.subtract_metabolites({met: reaction.get_coefficient(met)})
            except KeyError:
                pass

    benchmark(add_remove_metabolite)


def test_add_metabolite(model: Model) -> None:
    """Test metabolite addition to a reaction from an unsolved model."""
    with pytest.raises(ValueError):
        model.add_metabolites(Metabolite())
    with model:
        with model:
            reaction = model.reactions.get_by_id("PGI")
            reaction.add_metabolites({model.metabolites[0]: 1})
            assert model.metabolites[0] in reaction._metabolites
            fake_metabolite = Metabolite("fake")
            reaction.add_metabolites({fake_metabolite: 1})
            assert fake_metabolite in reaction._metabolites
            assert model.metabolites.has_id("fake")
            assert model.metabolites.get_by_id("fake") is fake_metabolite
            assert len(model._contexts[0]._history) == 0

    assert fake_metabolite._model is None
    assert fake_metabolite not in reaction._metabolites
    assert "fake" not in model.metabolites

    # Test adding by string
    with model:
        reaction.add_metabolites({"g6p_c": -1})  # already in reaction
        assert reaction._metabolites[model.metabolites.get_by_id("g6p_c")] == -2
        reaction.add_metabolites({"h_c": 1})
        assert reaction._metabolites[model.metabolites.get_by_id("h_c")] == 1
        with pytest.raises(KeyError):
            reaction.add_metabolites({"missing": 1})

    assert reaction._metabolites[model.metabolites.get_by_id("g6p_c")] == -1
    assert model.metabolites.h_c not in reaction._metabolites

    # Test combine=False
    reaction = model.reactions.get_by_id("ATPM")
    old_stoich = reaction._metabolites[model.metabolites.get_by_id("h2o_c")]
    with model:
        reaction.add_metabolites({"h2o_c": 2.5}, combine=False)
        assert reaction._metabolites[model.metabolites.get_by_id("h2o_c")] == 2.5

    assert reaction._metabolites[model.metabolites.get_by_id("h2o_c")] == old_stoich

    # Test adding to a new Reaction
    reaction = Reaction("test")
    assert len(reaction._metabolites) == 0
    reaction.add_metabolites({Metabolite("test_met"): -1})
    assert len(reaction._metabolites) == 1


@pytest.mark.parametrize("solver", stable_optlang)
def test_subtract_metabolite_benchmark(
    model: Model, benchmark, solver: Iterable
) -> None:
    """Benchmark metabolite deletion from a reaction."""
    benchmark(test_subtract_metabolite, model, solver)


@pytest.mark.parametrize("solver", stable_optlang)
def test_subtract_metabolite(model: Model, solver: Iterable) -> None:
    """Test metabolite deletion from a reaction associated with an unsolved model."""
    reaction = model.reactions.get_by_id("PGI")
    reaction.subtract_metabolites(reaction.metabolites)
    if not getattr(model, "solver", None):
        stable_optlang[solver].create_problem(model)
        assert len(reaction.metabolites) == 0


def test_mass_balance(model: Model) -> None:
    """Test mass balance of metabolites of a reaction."""
    reaction = model.reactions.get_by_id("PGI")
    # Should be balanced now
    assert len(reaction.check_mass_balance()) == 0
    # Should not be balanced after adding a hydrogen
    reaction.add_metabolites({model.metabolites.get_by_id("h_c"): 1})
    imbalance = reaction.check_mass_balance()
    assert "H" in imbalance
    assert imbalance["H"] == 1


def test_build_from_string(model: Model) -> None:
    """Test reaction building from string evaluation."""
    m = len(model.metabolites)
    pgi = model.reactions.get_by_id("PGI")
    old_bounds = pgi.bounds

    with model:
        pgi.reaction = "g6p_c --> f6p_c"
        assert pgi.lower_bound == 0

    assert pgi.bounds == old_bounds

    pgi.bounds = (0, 1000)
    assert pgi.bounds == (0, 1000)
    assert not pgi.reversibility
    pgi.reaction = "g6p_c <== f6p_c"
    assert pgi.upper_bound == 0
    assert pgi.reaction.strip() == "g6p_c <-- f6p_c"
    pgi.reaction = "g6p_c --> f6p_c + h2o_c"
    assert model.metabolites.h2o_c, pgi._metabolites

    with model:
        pgi.build_reaction_from_string("g6p_c --> f6p_c + foo", verbose=False)
        assert model.metabolites.h2o_c not in pgi._metabolites
        assert "foo" in model.metabolites
        assert model.metabolites.foo in pgi._metabolites
        assert len(model.metabolites) == m + 1

    assert model.metabolites.h2o_c in pgi._metabolites
    assert "foo" not in model.metabolites
    with pytest.raises(AttributeError):
        assert model.metabolites.foo
        assert len(model.metabolites) == m

    with model:
        old_bounds = config.bounds
        assert old_bounds == (-1000, 1000)
        config.bounds = (-5, 5)
        pgi.build_reaction_from_string("g6p_c <--> f6p_c + new", verbose=False)
        assert pgi.bounds == (-5, 5)
        config.bounds = old_bounds
        pgi.build_reaction_from_string("g6p_c --> f6p_c + new", verbose=False)
        assert pgi.bounds == (0, 1000)


def test_build_from_string_creating_metabolites() -> None:
    """Test that metabolites are created in the correct compartment."""
    # https://github.com/opencobra/cobrapy/issues/1418
    model = Model()
    reaction = Reaction("R1")
    model.add_reactions([reaction])
    reaction.build_reaction_from_string("[c]: a --> b")
    assert len(model.metabolites) == 2
    assert model.metabolites.get_by_id("a[c]").compartment == "c"
    assert model.metabolites.get_by_id("b[c]").compartment == "c"
    assert model.reactions.R1.compartments == set(["c"])


def test_bounds_setter(model: Model) -> None:
    """Test reaction bounds setter."""
    rxn = model.reactions.get_by_id("PGI")
    with pytest.raises(ValueError):
        rxn.bounds = (1, 0)


def test_copy(model: Model) -> None:
    """Test reaction copying."""
    PGI = model.reactions.PGI
    copied = PGI.copy()
    assert PGI is not copied
    assert PGI._model is model
    assert copied._model is not model
    # The copy should refer to different metabolites and genes
    for met in copied.metabolites:
        assert met is not model.metabolites.get_by_id(met.id)
        assert met.model is not model
    for gene in copied.genes:
        assert gene is not model.genes.get_by_id(gene.id)
        assert gene.model is not model
    assert len(model.get_associated_groups(copied.id)) == 0


def test_iadd(model: Model) -> None:
    """Test in-place addition of reaction."""
    PGI = model.reactions.PGI
    EX_h2o = model.reactions.EX_h2o_e
    original_PGI_gpr = PGI.gene_reaction_rule
    PGI += EX_h2o
    assert PGI.gene_reaction_rule == original_PGI_gpr
    assert PGI.metabolites[model.metabolites.h2o_e] == -1.0
    # Original should not change
    assert EX_h2o.gene_reaction_rule == ""
    assert EX_h2o.metabolites[model.metabolites.h2o_e] == -1.0
    # Add a reaction not in the model
    new_reaction = Reaction("test")
    new_reaction.add_metabolites({Metabolite("A"): -1, Metabolite("B"): 1})
    PGI += new_reaction
    assert PGI.gene_reaction_rule == original_PGI_gpr
    assert len(PGI.gene_reaction_rule) == 5
    # And vice versa
    new_reaction += PGI
    assert len(new_reaction.metabolites) == 5  # not
    assert len(new_reaction.genes) == 1
    assert new_reaction.gene_reaction_rule == original_PGI_gpr
    # Combine two GPRs
    model.reactions.ACKr += model.reactions.ACONTa
    expected_rule = "(b2296 or b3115 or b1849) and (b0118 or b1276)"
    assert model.reactions.ACKr.gene_reaction_rule == expected_rule
    assert len(model.reactions.ACKr.genes) == 5


def test_iadd_with_context(model: Model) -> None:
    """Test in-place addition of reaction is reversed with context."""
    PGI = model.reactions.PGI
    EX_h2o = model.reactions.EX_h2o_e
    original_PGI_gene_reaction_rule = PGI.gene_reaction_rule
    with model:
        PGI += EX_h2o
        assert PGI.gene_reaction_rule == original_PGI_gene_reaction_rule
        assert PGI.metabolites[model.metabolites.h2o_e] == -1.0
    assert PGI.gene_reaction_rule == original_PGI_gene_reaction_rule
    assert model.metabolites.h2o_e not in PGI.metabolites.keys()
    # Add a reaction not in the model
    new_reaction = Reaction("test")
    new_reaction.add_metabolites({Metabolite("A"): -1, Metabolite("B"): 1})
    with model:
        PGI += new_reaction
    assert PGI.gene_reaction_rule == original_PGI_gene_reaction_rule
    assert len(PGI.gene_reaction_rule) == 5
    # Combine two GPRs
    expected_rule = "(b2296 or b3115 or b1849) and (b0118 or b1276)"
    old_rule = model.reactions.ACKr.gene_reaction_rule
    with model:
        model.reactions.ACKr += model.reactions.ACONTa
        assert model.reactions.ACKr.gene_reaction_rule == expected_rule
        assert len(model.reactions.ACKr.genes) == 5
    assert model.reactions.ACKr.gene_reaction_rule == old_rule
    assert old_rule != expected_rule
    assert len(model.reactions.ACKr.genes) == 3


def test_add(model: Model) -> None:
    """Test reaction addition to model."""
    # Not in place addition should work on a copy
    new = model.reactions.PGI + model.reactions.EX_h2o_e
    assert new._model is not model
    assert len(new.metabolites) == 3
    # The copy should refer to different metabolites and genes
    # This currently fails because add_metabolites does not copy.
    # Should that be changed?
    # for met in new.metabolites:
    #    assert met is not model.metabolites.get_by_id(met.id)
    #    assert met.model is not model
    for gene in new.genes:
        assert gene is not model.genes.get_by_id(gene.id)
        assert gene.model is not model


def test_radd(model: Model) -> None:
    """Test __radd__ for a reaction."""
    new = sum([model.reactions.PGI, model.reactions.EX_h2o_e])
    assert new._model is not model
    assert len(new.metabolites) == 3


def test_mul(model: Model) -> None:
    """Test scalar multiplication of factors with a reaction."""
    new = model.reactions.PGI * 2
    assert set(new.metabolites.values()) == {-2, 2}


def test_sub(model: Model) -> None:
    """Test reaction subtraction."""
    new = model.reactions.PGI - model.reactions.EX_h2o_e
    assert new._model is not model
    assert len(new.metabolites) == 3


def test_removal_from_model_retains_bounds(model: Model) -> None:
    """Test reaction removal from a model, retains its bounds."""
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


def test_set_bounds_scenario_1(model: Model) -> None:
    """Test reaction bounds setting for a scenario."""
    acald_reaction = model.reactions.ACALD
    assert acald_reaction.lower_bound == -1000.0
    assert acald_reaction.upper_bound == 1000.0
    assert acald_reaction.forward_variable.lb == 0.0
    assert acald_reaction.forward_variable.ub == 1000.0
    assert acald_reaction.reverse_variable.lb == 0
    assert acald_reaction.reverse_variable.ub == 1000.0
    acald_reaction.bounds = (
        acald_reaction.lower_bound - 100,
        acald_reaction.lower_bound - 100,
    )
    assert acald_reaction.lower_bound == -1100.0
    assert acald_reaction.upper_bound == -1100.0
    assert acald_reaction.forward_variable.lb == 0
    assert acald_reaction.forward_variable.ub == 0
    assert acald_reaction.reverse_variable.lb == 1100.0
    assert acald_reaction.reverse_variable.ub == 1100.0
    acald_reaction.upper_bound = 100
    assert acald_reaction.lower_bound == -1100.0
    assert acald_reaction.upper_bound == 100
    assert acald_reaction.forward_variable.lb == 0
    assert acald_reaction.forward_variable.ub == 100
    assert acald_reaction.reverse_variable.lb == 0
    assert acald_reaction.reverse_variable.ub == 1100.0


def test_set_bounds_scenario_2(model: Model) -> None:
    """Test reaction bounds setting for a scenario."""
    acald_reaction = model.reactions.ACALD
    assert acald_reaction.lower_bound == -1000.0
    assert acald_reaction.upper_bound == 1000.0
    assert acald_reaction.forward_variable.lb == 0.0
    assert acald_reaction.forward_variable.ub == 1000.0
    assert acald_reaction.reverse_variable.lb == 0
    assert acald_reaction.reverse_variable.ub == 1000.0
    acald_reaction.bounds = (
        acald_reaction.upper_bound + 100,
        acald_reaction.upper_bound + 100,
    )
    assert acald_reaction.lower_bound == 1100.0
    assert acald_reaction.upper_bound == 1100.0
    assert acald_reaction.forward_variable.lb == 1100.0
    assert acald_reaction.forward_variable.ub == 1100.0
    assert acald_reaction.reverse_variable.lb == 0
    assert acald_reaction.reverse_variable.ub == 0
    acald_reaction.lower_bound = -100
    assert acald_reaction.lower_bound == -100.0
    assert acald_reaction.upper_bound == 1100.0
    assert acald_reaction.forward_variable.lb == 0
    assert acald_reaction.forward_variable.ub == 1100.0
    assert acald_reaction.reverse_variable.lb == 0
    assert acald_reaction.reverse_variable.ub == 100


def test_set_bounds_scenario_3(model: Model) -> None:
    """Test reaction bounds setting for a scenario."""
    reac = model.reactions.ACALD
    reac.bounds = (-10, -10)
    assert reac.lower_bound == -10
    assert reac.upper_bound == -10
    reac.bounds = (-9, -9)
    assert reac.lower_bound == -9
    assert reac.upper_bound == -9
    reac.bounds = (2, 2)
    assert reac.lower_bound == 2
    assert reac.upper_bound == 2
    reac.bounds = (-10, -10)
    assert reac.lower_bound == -10
    assert reac.upper_bound == -10
    reac.bounds = (-11, -11)
    assert reac.lower_bound == -11
    assert reac.upper_bound == -11
    reac.upper_bound = 2
    assert reac.lower_bound == -11
    assert reac.upper_bound == 2


def test_set_bounds_scenario_4(model: Model) -> None:
    """Test reaction bounds setting for a scenario."""
    reac = model.reactions.ACALD
    reac.bounds = (2, 2)
    assert reac.lower_bound == 2
    assert reac.upper_bound == 2
    assert reac.forward_variable.lb == 2
    assert reac.forward_variable.ub == 2
    reac.knock_out()
    reac.bounds = (-2, -2)
    assert reac.lower_bound == -2
    assert reac.upper_bound == -2
    assert reac.reverse_variable.lb == 2
    assert reac.reverse_variable.ub == 2


def test_set_upper_before_lower_bound_to_0(model: Model) -> None:
    """Test reaction bounds setting to zero."""
    model.reactions.GAPD.bounds = (0, 0)
    assert model.reactions.GAPD.lower_bound == 0
    assert model.reactions.GAPD.upper_bound == 0
    assert model.reactions.GAPD.forward_variable.lb == 0
    assert model.reactions.GAPD.forward_variable.ub == 0
    assert model.reactions.GAPD.reverse_variable.lb == 0
    assert model.reactions.GAPD.reverse_variable.ub == 0


def test_change_bounds(model: Model) -> None:
    """Test reaction bounds change."""
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


def test_make_irreversible(model: Model) -> None:
    """Test reaction irreversibility."""
    acald_reaction = model.reactions.ACALD
    assert acald_reaction.lower_bound == -1000.0
    assert acald_reaction.upper_bound == 1000.0
    assert acald_reaction.forward_variable.lb == 0.0
    assert acald_reaction.forward_variable.ub == 1000.0
    assert acald_reaction.reverse_variable.lb == 0
    assert acald_reaction.reverse_variable.ub == 1000.0
    acald_reaction.lower_bound = 0
    assert acald_reaction.lower_bound == 0
    assert acald_reaction.upper_bound == 1000.0
    assert acald_reaction.forward_variable.lb == 0
    assert acald_reaction.forward_variable.ub == 1000.0
    assert acald_reaction.reverse_variable.lb == 0
    assert acald_reaction.reverse_variable.ub == 0
    acald_reaction.lower_bound = -100
    assert acald_reaction.lower_bound == -100.0
    assert acald_reaction.upper_bound == 1000.0
    assert acald_reaction.forward_variable.lb == 0
    assert acald_reaction.forward_variable.ub == 1000.0
    assert acald_reaction.reverse_variable.lb == 0
    assert acald_reaction.reverse_variable.ub == 100


def test_make_reversible(model: Model) -> None:
    """Test reaction reversibility."""
    pfk_reaction = model.reactions.PFK
    assert pfk_reaction.lower_bound == 0.0
    assert pfk_reaction.upper_bound == 1000.0
    assert pfk_reaction.forward_variable.lb == 0.0
    assert pfk_reaction.forward_variable.ub == 1000.0
    assert pfk_reaction.reverse_variable.lb == 0
    assert pfk_reaction.reverse_variable.ub == 0
    pfk_reaction.lower_bound = -100.0
    assert pfk_reaction.lower_bound == -100.0
    assert pfk_reaction.upper_bound == 1000.0
    assert pfk_reaction.forward_variable.lb == 0
    assert pfk_reaction.forward_variable.ub == 1000.0
    assert pfk_reaction.reverse_variable.lb == 0
    assert pfk_reaction.reverse_variable.ub == 100.0
    pfk_reaction.lower_bound = 0
    assert pfk_reaction.lower_bound == 0
    assert pfk_reaction.upper_bound == 1000.0
    assert pfk_reaction.forward_variable.lb == 0
    assert pfk_reaction.forward_variable.ub == 1000.0
    assert pfk_reaction.reverse_variable.lb == 0
    assert pfk_reaction.reverse_variable.ub == 0


def test_make_irreversible_irreversible_to_the_other_side(model: Model) -> None:
    """Test reaction irreversibility to irreversibility."""
    pfk_reaction = model.reactions.PFK
    assert pfk_reaction.lower_bound == 0.0
    assert pfk_reaction.upper_bound == 1000.0
    assert pfk_reaction.forward_variable.lb == 0.0
    assert pfk_reaction.forward_variable.ub == 1000.0
    assert pfk_reaction.reverse_variable.lb == 0
    assert pfk_reaction.reverse_variable.ub == 0
    pfk_reaction.bounds = (-100.0, -100.0)
    assert pfk_reaction.forward_variable.lb == 0
    assert pfk_reaction.forward_variable.ub == 0
    assert pfk_reaction.reverse_variable.lb == 100
    assert pfk_reaction.reverse_variable.ub == 100
    pfk_reaction.lower_bound = -1000.0
    assert pfk_reaction.lower_bound == -1000.0
    assert pfk_reaction.upper_bound == -100.0
    assert pfk_reaction.forward_variable.lb == 0
    assert pfk_reaction.forward_variable.ub == 0
    assert pfk_reaction.reverse_variable.lb == 100
    assert pfk_reaction.reverse_variable.ub == 1000.0


def test_make_lhs_irreversible_reversible(model: Model) -> None:
    """Test reaction LHS irreversibility to reversibility."""
    rxn = Reaction("test")
    rxn.add_metabolites({model.metabolites[0]: -1.0, model.metabolites[1]: 1.0})
    rxn.bounds = (-1000.0, -100)
    model.add_reactions([rxn])
    assert rxn.lower_bound == -1000.0
    assert rxn.upper_bound == -100.0
    assert rxn.forward_variable.lb == 0.0
    assert rxn.forward_variable.ub == 0.0
    assert rxn.reverse_variable.lb == 100.0
    assert rxn.reverse_variable.ub == 1000.0
    rxn.upper_bound = 666.0
    assert rxn.lower_bound == -1000.0
    assert rxn.upper_bound == 666.0
    assert rxn.forward_variable.lb == 0.0
    assert rxn.forward_variable.ub == 666
    assert rxn.reverse_variable.lb == 0.0
    assert rxn.reverse_variable.ub == 1000.0


def test_model_less_reaction(model: Model) -> None:
    """Test model without reactions."""
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


def test_knockout(model: Model) -> None:
    """Test reaction knockouts."""
    original_bounds = dict()
    for reaction in model.reactions:
        original_bounds[reaction.id] = (reaction.lower_bound, reaction.upper_bound)
        reaction.knock_out()
        assert reaction.lower_bound == 0
        assert reaction.upper_bound == 0
    for k, (lb, ub) in original_bounds.items():
        model.reactions.get_by_id(k).bounds = (lb, ub)
    for reaction in model.reactions:
        assert reaction.lower_bound == original_bounds[reaction.id][0]
        assert reaction.upper_bound == original_bounds[reaction.id][1]
    with model:
        for reaction in model.reactions:
            original_bounds[reaction.id] = (reaction.lower_bound, reaction.upper_bound)
            reaction.knock_out()
            assert reaction.lower_bound == 0
            assert reaction.upper_bound == 0
    for reaction in model.reactions:
        assert reaction.lower_bound == original_bounds[reaction.id][0]
        assert reaction.upper_bound == original_bounds[reaction.id][1]


def test_reaction_without_model() -> None:
    """Test reaction without model association."""
    r = Reaction("blub")
    assert r.flux_expression is None
    assert r.forward_variable is None
    assert r.reverse_variable is None


def test_weird_left_to_right_reaction_issue(tiny_toy_model: Model) -> None:
    """Test absurd left to right reaction."""
    d1 = tiny_toy_model.reactions.get_by_id("ex1")
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


def test_one_left_to_right_reaction_set_positive_ub(tiny_toy_model: Model) -> None:
    """Test left to right reaction with positive upper bound."""
    d1 = tiny_toy_model.reactions.get_by_id("ex1")
    assert d1.reverse_variable.lb == 0
    assert d1.reverse_variable.ub == 1000
    assert d1._lower_bound == -1000
    assert d1.lower_bound == -1000
    assert d1._upper_bound == 0
    assert d1.upper_bound == 0
    assert d1.forward_variable.lb == 0
    assert d1.forward_variable.ub == 0
    d1.upper_bound = 0.1
    assert d1.forward_variable.lb == 0
    assert d1.forward_variable.ub == 0.1
    assert d1.reverse_variable.lb == 0
    assert d1.reverse_variable.ub == 1000
    assert d1._lower_bound == -1000
    assert d1.upper_bound == 0.1
    assert d1._lower_bound == -1000
    assert d1.upper_bound == 0.1


def test_irrev_reaction_set_negative_lb(model: Model) -> None:
    """Test reaction irreversibility with negative lower bound."""
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


def test_twist_irrev_right_to_left_reaction_to_left_to_right(model: Model) -> None:
    """Test irreversibility reversal from right to left to left to right."""
    assert not model.reactions.PFK.reversibility
    assert model.reactions.PFK.lower_bound == 0
    assert model.reactions.PFK.upper_bound == 1000.0
    assert model.reactions.PFK.forward_variable.lb == 0
    assert model.reactions.PFK.forward_variable.ub == 1000.0
    assert model.reactions.PFK.reverse_variable.lb == 0
    assert model.reactions.PFK.reverse_variable.ub == 0
    model.reactions.PFK.bounds = (-1000, 0)
    assert model.reactions.PFK.lower_bound == -1000
    assert model.reactions.PFK.upper_bound == 0
    assert model.reactions.PFK.forward_variable.lb == 0
    assert model.reactions.PFK.forward_variable.ub == 0
    assert model.reactions.PFK.reverse_variable.lb == 0
    assert model.reactions.PFK.reverse_variable.ub == 1000


def test_set_lb_higher_than_ub_sets_ub_to_new_lb(model: Model) -> None:
    """Test lower bound > upper bound makes upper bound to new lower bound."""
    for reaction in model.reactions:
        assert reaction.lower_bound <= reaction.upper_bound
        reaction.bounds = (reaction.upper_bound + 100, reaction.upper_bound + 100)
        assert reaction.lower_bound == reaction.upper_bound


def test_set_ub_lower_than_lb_sets_lb_to_new_ub(model: Model) -> None:
    """Test upper bound < lower bound makes lower bound to new upper bound."""
    for reaction in model.reactions:
        assert reaction.lower_bound <= reaction.upper_bound
        reaction.bounds = (reaction.lower_bound - 100, reaction.lower_bound - 100)
        assert reaction.lower_bound == reaction.upper_bound


def test_add_metabolites_combine_true(model: Model) -> None:
    """Test metabolite addition to reaction (with combine = True)."""
    test_metabolite = Metabolite("test")
    for reaction in model.reactions:
        reaction.add_metabolites({test_metabolite: -66}, combine=True)
        assert reaction.metabolites[test_metabolite] == -66
        assert (
            model.constraints["test"].get_linear_coefficients(
                [reaction.forward_variable]
            )[reaction.forward_variable]
            == -66
        )
        assert (
            model.constraints["test"].get_linear_coefficients(
                [reaction.reverse_variable]
            )[reaction.reverse_variable]
            == 66
        )
        already_included_metabolite = list(reaction.metabolites.keys())[0]
        previous_coefficient = reaction.get_coefficient(already_included_metabolite.id)
        reaction.add_metabolites({already_included_metabolite: 10}, combine=True)
        new_coefficient = previous_coefficient + 10
        assert reaction.metabolites[already_included_metabolite] == new_coefficient
        assert (
            model.constraints[already_included_metabolite.id].get_linear_coefficients(
                [reaction.forward_variable]
            )[reaction.forward_variable]
            == new_coefficient
        )
        assert (
            model.constraints[already_included_metabolite.id].get_linear_coefficients(
                [reaction.reverse_variable]
            )[reaction.reverse_variable]
            == -new_coefficient
        )


@pytest.mark.xfail(reason="non-deterministic test")
def test_add_metabolites_combine_false(model: Model) -> None:
    """Test metabolite addition to reaction (with combine = False)."""
    test_metabolite = Metabolite("test")
    for reaction in model.reactions:
        reaction.add_metabolites({test_metabolite: -66}, combine=False)
        assert reaction.metabolites[test_metabolite] == -66
        assert model.constraints["test"].expression.has(
            -66.0 * reaction.forward_variable
        )
        assert model.constraints["test"].expression.has(
            66.0 * reaction.reverse_variable
        )
        already_included_metabolite = list(reaction.metabolites.keys())[0]
        reaction.add_metabolites({already_included_metabolite: 10}, combine=False)
        assert reaction.metabolites[already_included_metabolite] == 10
        assert model.constraints[already_included_metabolite.id].expression.has(
            10 * reaction.forward_variable
        )
        assert model.constraints[already_included_metabolite.id].expression.has(
            -10 * reaction.reverse_variable
        )


def test_reaction_imul(model: Model) -> None:
    """Test in-place scalar factor multiplication to reaction."""
    with model:
        model.reactions.EX_glc__D_e *= 100
        assert (
            model.constraints.glc__D_e.expression.coeff(model.variables.EX_glc__D_e)
            == -100.0
        )
        assert model.reactions.EX_glc__D_e.reaction == "100.0 glc__D_e <=> "

    assert (
        model.constraints.glc__D_e.expression.coeff(model.variables.EX_glc__D_e) == -1.0
    )
    assert model.reactions.EX_glc__D_e.reaction == "glc__D_e <=> "

    with model:
        model.reactions.EX_glc__D_e *= -2
        assert model.reactions.EX_glc__D_e.bounds == (-1000.0, 10.0)
        assert model.reactions.EX_glc__D_e.reaction == " <=> 2.0 glc__D_e"

    assert model.reactions.EX_glc__D_e.bounds == (-10, 1000.0)
    assert model.reactions.EX_glc__D_e.reaction == "glc__D_e <=> "


# def test_pop(model):
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


def test_remove_from_model(model: Model) -> None:
    """Test reaction removal from model."""
    pgi = model.reactions.PGI
    g6p = model.metabolites.g6p_c
    pgi_flux = model.optimize().fluxes["PGI"]
    assert abs(pgi_flux) > 1e-6

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
    assert np.isclose(pgi_flux, model.optimize().fluxes["PGI"])


def test_change_id_is_reflected_in_solver(model: Model) -> None:
    """Test reaction ID change reflection in solver."""
    for i, reaction in enumerate(model.reactions):
        old_reaction_id = reaction.id
        assert model.variables[old_reaction_id].name == old_reaction_id
        assert old_reaction_id in model.variables
        new_reaction_id = reaction.id + "_" + str(i)
        reaction.id = new_reaction_id
        assert reaction.id == new_reaction_id
        assert not (old_reaction_id in model.variables)
        assert reaction.id in model.variables
        assert reaction.reverse_id in model.variables
        name = model.variables[reaction.id].name
        assert name == reaction.id


def test_repr_html_(model: Model) -> None:
    """Test __repr_html__ functionality."""
    assert "<table>" in model.reactions[0]._repr_html_()


def test_compartment_changes(model: Model) -> None:
    """Test reaction compartment change."""
    rxn = model.reactions.EX_ac_e
    assert rxn.reactants[0].compartment in rxn.compartments
    rxn.reactants[0].compartment = "blub"
    assert rxn.reactants[0].compartment in rxn.compartments


def test_gpr_serialization(model: Model) -> None:
    """Verify that reactions GPRs are serialized compactly as str."""
    state = model.reactions[0].__getstate__()
    assert isinstance(state["_gpr"], str)


def test_gpr_serialization_backwards_compatibility(model: Model) -> None:
    """Verify that GPR serialization is backwards compatible."""
    state = model.reactions[0].__getstate__()
    print(state)
    state["gene_reaction_rule"] = state["_gpr"]  # old format
    del state["_gpr"]
    rxn = Reaction("test")
    rxn.__setstate__(state)
    assert isinstance(rxn.gpr, GPR)
