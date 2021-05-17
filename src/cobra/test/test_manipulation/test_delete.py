"""Test functionalities of model component pruning functions."""

from itertools import chain
from typing import List, Set

from cobra.core import Gene, Metabolite, Model, Reaction
from cobra.manipulation import (
    delete_model_genes,
    find_gene_knockout_reactions,
    get_compiled_gene_reaction_rules,
    prune_unused_metabolites,
    prune_unused_reactions,
    remove_genes,
    undelete_model_genes,
)


def test_prune_unused_metabolites_output_type(model: Model) -> None:
    """Test the output type of unused metabolites pruning."""
    [model.reactions.get_by_id(x).remove_from_model() for x in ["RPI", "RPE", "GND"]]
    model_pruned, unused = prune_unused_metabolites(model)
    assert isinstance(model_pruned, Model)
    # test that the output contains metabolite objects
    assert isinstance(unused[0], Metabolite)


def test_prune_unused_metabolites_sanity(model: Model) -> None:
    """Test the sanity of unused metabolites pruning."""
    metabolite1 = model.metabolites.ru5p__D_c
    metabolite2 = model.metabolites.akg_e
    metabolite3 = model.metabolites.akg_c
    reactions = set(
        chain(metabolite1.reactions, metabolite2.reactions, metabolite3.reactions)
    )
    model.remove_reactions(reactions)
    model_pruned, unused = prune_unused_metabolites(model)
    assert metabolite1 in model.metabolites
    assert metabolite2 in model.metabolites
    assert metabolite3 in model.metabolites
    # test that the unused metabolites are not used in the model
    assert metabolite1 not in model_pruned.metabolites
    assert metabolite2 not in model_pruned.metabolites
    assert metabolite3 not in model_pruned.metabolites


def test_prune_unused_reactions_output_type(model: Model) -> None:
    """Test the output type of unused reactions pruning."""
    reaction = Reaction("foo")
    model.add_reactions([reaction])
    model_pruned, unused = prune_unused_reactions(model)
    assert isinstance(model_pruned, Model)
    # test that the output contains reaction objects
    assert isinstance(unused[0], Reaction)


def test_prune_unused_rxns_functionality(model: Model) -> None:
    """Test the sanity of unused reactions pruning."""
    for x in ["foo1", "foo2", "foo3"]:
        model.add_reactions([Reaction(x)])

    model_pruned, unused = prune_unused_reactions(model)
    assert "foo1" in model.reactions
    assert "foo2" in model.reactions
    assert "foo3" in model.reactions
    # test that the unused reactions are not used in the model
    assert "foo1" not in model_pruned.reactions
    assert "foo2" not in model_pruned.reactions
    assert "foo3" not in model_pruned.reactions


def _find_gene_knockout_reactions_fast(
    m: Model, gene_list: List[Gene]
) -> List[Reaction]:
    """Quickly find gene knockout reactions."""
    compiled_rules = get_compiled_gene_reaction_rules(m)
    return find_gene_knockout_reactions(
        m, gene_list, compiled_gene_reaction_rules=compiled_rules
    )


def _get_removed(m: Model) -> Set[str]:
    """Get trimmed reactions."""
    return {x.id for x in m._trimmed_reactions}


def _gene_knockout_computation(
    m: Model, gene_ids: List[str], expected_reaction_ids: List[str]
) -> None:
    """Compute gene knockout."""
    genes = [m.genes.get_by_id(i) for i in gene_ids]
    expected_reactions = {m.reactions.get_by_id(i) for i in expected_reaction_ids}
    removed1 = set(find_gene_knockout_reactions(m, genes))
    removed2 = set(_find_gene_knockout_reactions_fast(m, genes))
    assert removed1 == expected_reactions
    assert removed2 == expected_reactions
    delete_model_genes(m, gene_ids, cumulative_deletions=False)
    assert _get_removed(m) == expected_reaction_ids
    undelete_model_genes(m)


def test_gene_knockout(salmonella: Model) -> None:
    """Test gene knockout."""
    gene_list = ["STM1067", "STM0227"]
    dependent_reactions = {
        "3HAD121",
        "3HAD160",
        "3HAD80",
        "3HAD140",
        "3HAD180",
        "3HAD100",
        "3HAD181",
        "3HAD120",
        "3HAD60",
        "3HAD141",
        "3HAD161",
        "T2DECAI",
        "3HAD40",
    }
    _gene_knockout_computation(salmonella, gene_list, dependent_reactions)
    _gene_knockout_computation(salmonella, ["STM4221"], {"PGI"})
    _gene_knockout_computation(salmonella, ["STM1746.S"], {"4PEPTabcpp"})
    # test cumulative behavior
    delete_model_genes(salmonella, gene_list[:1])
    delete_model_genes(salmonella, gene_list[1:], cumulative_deletions=True)
    delete_model_genes(salmonella, ["STM4221"], cumulative_deletions=True)
    dependent_reactions.add("PGI")
    assert _get_removed(salmonella) == dependent_reactions
    # non-cumulative following cumulative
    delete_model_genes(salmonella, ["STM4221"], cumulative_deletions=False)
    assert _get_removed(salmonella) == {"PGI"}
    # make sure on reset that the bounds are correct
    reset_bound = salmonella.reactions.get_by_id("T2DECAI").upper_bound
    assert reset_bound == 1000.0
    # test computation when gene name is a subset of another
    test_model = Model()
    test_reaction_1 = Reaction("test1")
    test_reaction_1.gene_reaction_rule = "eggs or (spam and eggspam)"
    test_model.add_reactions([test_reaction_1])
    _gene_knockout_computation(test_model, ["eggs"], set())
    _gene_knockout_computation(test_model, ["eggs", "spam"], {"test1"})
    # test computation with nested boolean expression
    test_reaction_1.gene_reaction_rule = "g1 and g2 and (g3 or g4 or (g5 and g6))"
    _gene_knockout_computation(test_model, ["g3"], set())
    _gene_knockout_computation(test_model, ["g1"], {"test1"})
    _gene_knockout_computation(test_model, ["g5"], set())
    _gene_knockout_computation(test_model, ["g3", "g4", "g5"], {"test1"})
    # test computation when gene names are python expressions
    test_reaction_1.gene_reaction_rule = "g1 and (for or in)"
    _gene_knockout_computation(test_model, ["for", "in"], {"test1"})
    _gene_knockout_computation(test_model, ["for"], set())
    test_reaction_1.gene_reaction_rule = "g1 and g2 and g2.conjugate"
    _gene_knockout_computation(test_model, ["g2"], {"test1"})
    _gene_knockout_computation(test_model, ["g2.conjugate"], {"test1"})
    test_reaction_1.gene_reaction_rule = "g1 and (try:' or 'except:1)"
    _gene_knockout_computation(test_model, ["try:'"], set())
    _gene_knockout_computation(test_model, ["try:'", "'except:1"], {"test1"})


def test_remove_genes() -> None:
    """Test gene removal."""
    m = Model("test")
    m.add_reactions([Reaction("r" + str(i + 1)) for i in range(8)])
    assert len(m.reactions) == 8
    rxns = m.reactions
    rxns.r1.gene_reaction_rule = "(a and b) or (c and a)"
    rxns.r2.gene_reaction_rule = "(a and b and d and e)"
    rxns.r3.gene_reaction_rule = "(a and b) or (b and c)"
    rxns.r4.gene_reaction_rule = "(f and b) or (b and c)"
    rxns.r5.gene_reaction_rule = "x"
    rxns.r6.gene_reaction_rule = "y"
    rxns.r7.gene_reaction_rule = "x or     z"
    rxns.r8.gene_reaction_rule = ""
    assert "a" in m.genes
    assert "x" in m.genes
    remove_genes(m, ["a"], remove_reactions=False)
    assert "a" not in m.genes
    assert "x" in m.genes
    assert rxns.r1.gene_reaction_rule == ""
    assert rxns.r2.gene_reaction_rule == ""
    assert rxns.r3.gene_reaction_rule == "b and c"
    assert rxns.r4.gene_reaction_rule == "(f and b) or (b and c)"
    assert rxns.r5.gene_reaction_rule == "x"
    assert rxns.r6.gene_reaction_rule == "y"
    assert rxns.r7.genes == {m.genes.x, m.genes.z}
    assert rxns.r8.gene_reaction_rule == ""
    remove_genes(m, ["x"], remove_reactions=True)
    assert len(m.reactions) == 7
    assert "r5" not in m.reactions
    assert "x" not in m.genes
    assert rxns.r1.gene_reaction_rule == ""
    assert rxns.r2.gene_reaction_rule == ""
    assert rxns.r3.gene_reaction_rule == "b and c"
    assert rxns.r4.gene_reaction_rule == "(f and b) or (b and c)"
    assert rxns.r6.gene_reaction_rule == "y"
    assert rxns.r7.gene_reaction_rule == "z"
    assert rxns.r7.genes == {m.genes.z}
    assert rxns.r8.gene_reaction_rule == ""
