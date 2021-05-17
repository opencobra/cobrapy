"""Test functionalities of model component modifcations."""

from cobra.core import Model
from cobra.manipulation import escape_ID, rename_genes


def test_escape_ids(model: Model) -> None:
    """Test model component IDs' SBML compliance."""
    model.reactions.PGI.gene_reaction_rule = "a.b or c"
    assert "a.b" in model.genes
    escape_ID(model)
    assert "a.b" not in model.genes


def test_rename_genes(model: Model) -> None:
    """Test gene renaming functionality."""
    original_name = model.genes.b1241.name
    rename_dict = {
        "b1241": "foo",
        "hello": "world",
        "b3115": "b3115",
        "b2465": "b3919",
        "bar": "2935",
    }
    rename_genes(model, rename_dict)
    for i in rename_dict.keys():
        if i not in rename_dict.values():
            assert i not in model.genes

    assert "b3115" in model.genes
    assert "foo" in model.genes
    assert "world" not in model.genes
    # make sure the object name was preserved
    assert model.genes.foo.name == original_name
    # make sure the reactions are correct
    assert len(model.genes.foo.reactions) == 2
    assert model.reactions.ACALD.gene_reaction_rule == "b0351 or foo"
    assert model.reactions.TPI.gene_reaction_rule == "b3919"
    assert model.reactions.TPI.genes == {model.genes.b3919}
    assert model.reactions.TKT1.gene_reaction_rule == "b2935 or b3919"
    assert model.reactions.TKT1.genes == {model.genes.b2935, model.genes.b3919}
    assert model.genes.b3919.reactions == {
        model.reactions.get_by_id(i) for i in ("TKT1", "TKT2", "TPI")
    }
