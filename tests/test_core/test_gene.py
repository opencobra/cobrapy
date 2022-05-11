"""Test functions of gene.py"""


def test_repr_html_(model) -> None:
    assert "<table>" in model.genes[0]._repr_html_()
