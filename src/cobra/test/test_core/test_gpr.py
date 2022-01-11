"""Test functions of cobra.core.gene.GPR ."""
import warnings
from ast import parse as ast_parse

import pytest

from cobra.core.gene import GPR, ast2str, eval_gpr, parse_gpr


def test_gpr():
    gpr1 = GPR().from_string(string_gpr="")
    assert isinstance(gpr1.body, list)
    assert len(gpr1.body) == 0


def test_empty_gpr() -> None:
    gpr1 = GPR()
    assert len(gpr1._genes) == 0
    assert len(gpr1.genes) == 0
    gpr1.update_genes()
    assert len(gpr1.genes) == 0
    assert gpr1.to_string() == ""
    assert gpr1.eval()
    gpr2 = gpr1.copy()
    assert not hasattr(gpr1, "body")
    gpr1 = GPR(None)
    assert len(gpr1._genes) == 0
    assert len(gpr1.genes) == 0
    gpr1.update_genes()
    assert len(gpr1.genes) == 0
    assert gpr1.to_string() == ""
    assert gpr1.eval()
    assert not hasattr(gpr1, "body")
    gpr1 = gpr1.from_string("a")
    assert len(gpr1._genes) == 1
    assert len(gpr1.genes) == 1
    gpr1.update_genes()
    assert len(gpr1.genes) == 1
    assert gpr1.genes == {"a"}
    assert gpr1.to_string() == "a"
    assert gpr1.eval()
    assert hasattr(gpr1, "body")


def test_and_gpr():
    gpr1 = GPR.from_string("a & b")
    assert len(gpr1._genes) == 2
    assert len(gpr1.genes) == 2  # This will test update_genes()
    gpr1.update_genes()
    assert len(gpr1.genes) == 2
    assert gpr1.genes == {"a", "b"}
    assert gpr1.to_string() == "a and b"
    assert gpr1.eval()
    assert hasattr(gpr1, "body")
    gpr1 = gpr1.from_string("a and b")
    assert len(gpr1._genes) == 2
    assert len(gpr1.genes) == 2  # This will test update_genes()
    gpr1.update_genes()
    assert len(gpr1.genes) == 2
    assert gpr1.to_string() == "a and b"
    assert gpr1.eval()
    assert hasattr(gpr1, "body")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gpr1.from_string("a AND b")
        assert len(gpr1._genes) == 2
        assert len(gpr1.genes) == 2  # This will test update_genes()
        gpr1.update_genes()
        assert len(gpr1.genes) == 2
        assert gpr1.to_string() == "a and b"
        assert gpr1.eval()
        assert hasattr(gpr1, "body")
    gpr1 = gpr1.from_string("a & b & c")
    assert len(gpr1._genes) == 3
    assert len(gpr1.genes) == 3  # This will test update_genes()
    gpr1.update_genes()
    assert len(gpr1.genes) == 3
    assert gpr1.genes == {"a", "b", "c"}
    assert gpr1.to_string() == "(a and b) and c"
    assert gpr1.eval()
    assert not gpr1.eval("a")
    assert not gpr1.eval("b")
    assert not gpr1.eval("c")
    assert not gpr1.eval(["a", "b"])
    assert hasattr(gpr1, "body")
    gpr1 = gpr1.from_string("a")
    assert len(gpr1.genes) == 1
    gpr1.update_genes()
    assert len(gpr1.genes) == 1
    assert gpr1.to_string() == "a"
    assert gpr1.eval()
    assert hasattr(gpr1, "body")


def test_or_gpr():
    gpr1 = GPR.from_string("a | b")
    assert len(gpr1._genes) == 2
    assert len(gpr1.genes) == 2
    gpr1.update_genes()
    assert len(gpr1.genes) == 2
    assert gpr1.genes == {"a", "b"}
    assert gpr1.to_string() == "a or b"
    assert gpr1.eval()
    assert hasattr(gpr1, "body")
    gpr1 = gpr1.from_string("a or b")
    assert len(gpr1._genes) == 2
    assert len(gpr1.genes) == 2  # This will test update_genes()
    gpr1.update_genes()
    assert len(gpr1.genes) == 2
    assert gpr1.to_string() == "a or b"
    assert gpr1.eval()
    assert hasattr(gpr1, "body")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gpr1 = gpr1.from_string("a OR b")
        assert len(gpr1._genes) == 2
        assert len(gpr1.genes) == 2  # This will test update_genes()
        gpr1.update_genes()
        assert len(gpr1.genes) == 2
        assert gpr1.to_string() == "a or b"
        assert gpr1.eval()
        assert gpr1.eval(["a"])
        assert gpr1.eval(["b"])
        assert not gpr1.eval(["a", "b"])
        assert hasattr(gpr1, "body")
    gpr1 = gpr1.from_string("a | b | c")
    assert len(gpr1._genes) == 3
    assert len(gpr1.genes) == 3  # This will test update_genes()
    gpr1.update_genes()
    assert len(gpr1.genes) == 3
    assert gpr1.genes == {"a", "b", "c"}
    assert gpr1.to_string() == "(a or b) or c"
    assert gpr1.eval()
    assert gpr1.eval("a")
    assert gpr1.eval("b")
    assert gpr1.eval("c")
    assert gpr1.eval(["a", "b"])
    assert not gpr1.eval(["a", "b", "c"])
    assert hasattr(gpr1, "body")
    gpr1 = gpr1.from_string("a")
    assert len(gpr1.genes) == 1
    gpr1.update_genes()
    assert len(gpr1.genes) == 1
    assert gpr1.to_string() == "a"
    assert gpr1.eval()
    assert hasattr(gpr1, "body")


def test_complicated_gpr():
    gpr1 = GPR.from_string("(a | b) & c")
    assert len(gpr1._genes) == 3
    assert len(gpr1.genes) == 3
    gpr1.update_genes()
    assert len(gpr1.genes) == 3
    assert gpr1.to_string() == "(a or b) and c"
    assert gpr1.eval()
    assert hasattr(gpr1, "body")
    gpr1 = gpr1.from_string("(a or b) and c")
    assert len(gpr1._genes) == 3
    assert len(gpr1.genes) == 3  # This will test update_genes()
    gpr1.update_genes()
    assert len(gpr1.genes) == 3
    assert gpr1.to_string() == "(a or b) and c"
    assert gpr1.eval()
    assert hasattr(gpr1, "body")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gpr1 = gpr1.from_string("(a OR b) AND c")
        assert len(gpr1._genes) == 3
        assert len(gpr1.genes) == 3  # This will test update_genes()
        gpr1.update_genes()
        assert len(gpr1.genes) == 3
        assert gpr1.genes == {"a", "b", "c"}
        assert gpr1.to_string() == "(a or b) and c"
        assert gpr1.eval()
        assert hasattr(gpr1, "body")
        assert gpr1.eval()
        assert gpr1.eval("a")
        assert gpr1.eval("b")
        assert not gpr1.eval("c")
        assert not gpr1.eval(["a", "b"])
        assert not gpr1.eval(["a", "b", "c"])
    gpr1 = GPR.from_string("a & c | b & c")
    assert len(gpr1._genes) == 3
    assert len(gpr1.genes) == 3  # This will test update_genes()
    gpr1.update_genes()
    assert len(gpr1.genes) == 3
    assert gpr1.to_string() == "(a and c) or (b and c)"
    assert gpr1.eval()
    assert hasattr(gpr1, "body")
    assert gpr1.eval()
    assert gpr1.eval("a")
    assert gpr1.eval("b")
    assert not gpr1.eval("c")
    assert not gpr1.eval(["a", "b"])
    assert not gpr1.eval(["a", "b", "c"])


def test_gpr_from_ast() -> None:
    string_to_ast = "a"
    ast_tree = ast_parse(string_to_ast, "<string>", "eval")
    gpr1 = GPR(ast_tree)
    assert len(gpr1._genes) == 1
    assert len(gpr1.genes) == 1
    gpr1.update_genes()
    assert len(gpr1.genes) == 1
    assert gpr1.genes == {"a"}
    assert gpr1.to_string() == "a"
    assert gpr1.eval()
    assert not gpr1.eval("a")
    assert hasattr(gpr1, "body")
    string_to_ast = "a | b"
    ast_tree = ast_parse(string_to_ast, "<string>", "eval")
    gpr1 = GPR(ast_tree)
    assert len(gpr1._genes) == 2
    assert len(gpr1.genes) == 2
    gpr1.update_genes()
    assert len(gpr1.genes) == 2
    assert gpr1.genes == {"a", "b"}
    assert gpr1.to_string() == "a or b"
    assert gpr1.eval()
    assert gpr1.eval(["a"])
    assert gpr1.eval(["b"])
    assert not gpr1.eval(["a", "b"])
    assert hasattr(gpr1, "body")
    string_to_ast = "a and b"
    ast_tree = ast_parse(string_to_ast, "<string>", "eval")
    gpr1 = GPR(ast_tree)
    assert len(gpr1._genes) == 2
    assert len(gpr1.genes) == 2
    gpr1.update_genes()
    assert gpr1.genes == {"a", "b"}
    assert len(gpr1.genes) == 2
    assert gpr1.to_string() == "a and b"
    assert gpr1.eval()
    assert not gpr1.eval(["a"])
    assert not gpr1.eval(["b"])
    assert not gpr1.eval(["a", "b"])
    assert hasattr(gpr1, "body")


def test_wrong_input_gpr():
    with pytest.raises(TypeError):
        GPR.from_string(["a", "b"])
    with pytest.raises(TypeError):
        GPR(["a", "b"])
    with pytest.raises(TypeError):
        GPR({"a", "b"})
    with pytest.warns(SyntaxWarning):
        GPR.from_string("a |")
    with pytest.warns(SyntaxWarning):
        gpr1 = GPR.from_string("a & ")
        assert not hasattr(gpr1, "body")
        assert len(gpr1.genes) == 0


def test_deprecated_gpr():
    gpr1 = GPR.from_string("(a | b) & c")
    with pytest.deprecated_call():
        assert ast2str(gpr1) == "(a or b) and c"
    with pytest.deprecated_call():
        assert eval_gpr(gpr1, {})
    with pytest.deprecated_call():
        assert eval_gpr(gpr1, {"a"})
    with pytest.deprecated_call():
        assert eval_gpr(gpr1, {"b"})
    with pytest.deprecated_call():
        assert not eval_gpr(gpr1, {"c"})
    with pytest.deprecated_call():
        gpr1, genes = parse_gpr("(a | b) & c")
        assert genes == {"a", "b", "c"}
    with pytest.deprecated_call():
        assert ast2str(gpr1) == "(a or b) and c"
    with pytest.deprecated_call():
        assert eval_gpr(gpr1, {})
    with pytest.deprecated_call():
        assert eval_gpr(gpr1, {"a"})
    with pytest.deprecated_call():
        assert eval_gpr(gpr1, {"b"})
    with pytest.deprecated_call():
        assert not eval_gpr(gpr1, {"c"})
