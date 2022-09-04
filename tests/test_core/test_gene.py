"""Test functions of gene.py ."""


from cobra.core import Model


def test_repr_html_(model: Model) -> None:
    """Test HTML represenation is correct for a gene.

    Parameters
    ----------
    model : cobra.Model
        The textbook model.

    """
    assert "<table>" in model.genes[0]._repr_html_()
