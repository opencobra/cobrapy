"""Test functionalities of model component pruning functions."""

from itertools import chain

from cobra.core import GPR, Metabolite, Model, Reaction
from cobra.manipulation import (
    knock_out_model_genes,
    prune_unused_metabolites,
    prune_unused_reactions,
    remove_genes,
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


def test_gene_knockout(salmonella: Model) -> None:
    """Test gene knockout with and without context."""
    gene_list = ["STM1067", "STM0227"]
    dependent_reactions = [
        "3HAD100",
        "3HAD120",
        "3HAD121",
        "3HAD140",
        "3HAD141",
        "3HAD160",
        "3HAD161",
        "3HAD180",
        "3HAD181",
        "3HAD40",
        "3HAD60",
        "3HAD80",
        "T2DECAI",
    ]
    orig_gene_len = len(salmonella.genes)
    orig_bounds = salmonella.reactions.list_attr("bounds")
    with salmonella:
        expected_reactions = [
            salmonella.reactions.get_by_id(r) for r in dependent_reactions
        ]
        knocked_out_reactions = knock_out_model_genes(salmonella, gene_list)
        assert set(expected_reactions) == set(knocked_out_reactions)
    assert len(salmonella.genes) == orig_gene_len
    assert salmonella.reactions.list_attr("bounds") == orig_bounds
    with salmonella:
        expected_reactions = [salmonella.reactions.get_by_id("PGI")]
        knocked_out_reactions = knock_out_model_genes(salmonella, ["STM4221"])
        assert expected_reactions == knocked_out_reactions
    with salmonella:
        expected_reactions = [salmonella.reactions.get_by_id("PGI")]
        knocked_out_reactions = knock_out_model_genes(
            salmonella, [salmonella.genes.get_by_id("STM4221")]
        )
        assert expected_reactions == knocked_out_reactions
    with salmonella:
        expected_reactions = [salmonella.reactions.get_by_id("PGI")]
        knocked_out_reactions = knock_out_model_genes(
            salmonella, [salmonella.genes.index("STM4221")]
        )
        assert expected_reactions == knocked_out_reactions
    with salmonella:
        expected_reactions = [salmonella.reactions.get_by_id("4PEPTabcpp")]
        knocked_out_reactions = knock_out_model_genes(salmonella, ["STM1746.S"])
        assert expected_reactions == knocked_out_reactions
    knocked_out_reactions = knock_out_model_genes(salmonella, gene_list)
    assert len(knocked_out_reactions) == 13
    expected_reactions = [
        salmonella.reactions.get_by_id(r) for r in dependent_reactions
    ]
    assert set(knocked_out_reactions) == set(expected_reactions)
    knocked_out_reactions.extend(knock_out_model_genes(salmonella, ["STM4221"]))
    expected_reactions.append(salmonella.reactions.get_by_id("PGI"))
    assert set(knocked_out_reactions) == set(expected_reactions)
    # test computation when gene name is a subset of another
    test_model = Model()
    test_reaction_1 = Reaction("test1")
    test_reaction_1.gene_reaction_rule = "eggs or (spam and eggspam)"
    test_model.add_reactions([test_reaction_1])
    with test_model:
        knocked_out_reactions = knock_out_model_genes(test_model, ["eggs"])
        assert knocked_out_reactions == list()
    with test_model:
        knocked_out_reactions = knock_out_model_genes(test_model, ["eggs", "spam"])
        expected_reactions = [test_model.reactions.get_by_id("test1")]
        assert set(knocked_out_reactions) == set(expected_reactions)
    # test computation with nested boolean expression
    test_reaction_1.gene_reaction_rule = "g1 and g2 and (g3 or g4 or (g5 and g6))"
    with test_model:
        knocked_out_reactions = knock_out_model_genes(test_model, ["g3"])
        assert knocked_out_reactions == list()
    with test_model:
        knocked_out_reactions = knock_out_model_genes(test_model, ["g1"])
        assert knocked_out_reactions == [test_reaction_1]
    with test_model:
        knocked_out_reactions = knock_out_model_genes(test_model, ["g5"])
        assert knocked_out_reactions == list()
    with test_model:
        knocked_out_reactions = knock_out_model_genes(test_model, ["g3", "g4", "g5"])
        assert knocked_out_reactions == [test_reaction_1]
    # test computation when gene names are python expressions
    test_reaction_1.gene_reaction_rule = "g1 and (for or in)"
    with test_model:
        knocked_out_reactions = knock_out_model_genes(test_model, ["for", "in"])
        assert knocked_out_reactions == [test_reaction_1]
    with test_model:
        knocked_out_reactions = knock_out_model_genes(test_model, ["for"])
        assert knocked_out_reactions == list()
    test_reaction_1.gene_reaction_rule = "g1 and g2 and g2.conjugate"
    with test_model:
        knocked_out_reactions = knock_out_model_genes(test_model, ["g2"])
        assert knocked_out_reactions == [test_reaction_1]
    with test_model:
        knocked_out_reactions = knock_out_model_genes(test_model, ["g2.conjugate"])
        assert knocked_out_reactions == [test_reaction_1]
    test_reaction_1.gene_reaction_rule = "g1 and (try:' or 'except:1)"
    with test_model:
        knocked_out_reactions = knock_out_model_genes(test_model, ["try:'"])
        assert knocked_out_reactions == []
    with test_model:
        knocked_out_reactions = knock_out_model_genes(
            test_model, ["try:'", "'except:1"]
        )
        assert knocked_out_reactions == [test_reaction_1]


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


def test_remove_genes_with_context() -> None:
    """Test gene removal is reversed in context."""
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
    with m:
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
    assert "a" in m.genes
    assert "x" in m.genes
    assert rxns.r1.gene_reaction_rule == "(a and b) or (c and a)"
    assert rxns.r2.gpr == GPR.from_string("(a and b and d and e)")
    assert rxns.r3.gene_reaction_rule == "(a and b) or (b and c)"
    assert rxns.r4.gene_reaction_rule == "(f and b) or (b and c)"
    assert rxns.r5.gene_reaction_rule == "x"
    assert rxns.r6.gene_reaction_rule == "y"
    assert rxns.r7.genes == {m.genes.x, m.genes.z}
    assert rxns.r8.gene_reaction_rule == ""
    with m:
        remove_genes(m, ["a"], remove_reactions=False)
        assert "a" not in m.genes
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
    assert "a" in m.genes
    assert "x" in m.genes
    assert len(m.reactions) == 8
    assert rxns.r1.gene_reaction_rule == "(a and b) or (c and a)"
    assert rxns.r2.gpr == GPR.from_string("(a and b and d and e)")
    assert rxns.r3.gene_reaction_rule == "(a and b) or (b and c)"
    assert rxns.r4.gene_reaction_rule == "(f and b) or (b and c)"
    assert rxns.r5.gene_reaction_rule == "x"
    assert rxns.r6.gene_reaction_rule == "y"
    assert rxns.r7.genes == {m.genes.x, m.genes.z}
    assert rxns.r8.gene_reaction_rule == ""
