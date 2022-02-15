"""Test functions of cobra.core.gene.GPR ."""
import itertools
import warnings
from ast import parse as ast_parse

import pytest

from cobra.core.gene import GPR, ast2str, eval_gpr, parse_gpr


def test_gpr():
    gpr1 = GPR()
    assert len(gpr1.genes) == 0
    gpr1.update_genes()
    assert len(gpr1.genes) == 0
    assert gpr1.to_string() == ""
    assert gpr1.eval()
    gpr2 = gpr1.copy()
    assert len(gpr1.genes) == 0


@pytest.mark.parametrize("test_input", ["", "", None])
def test_empty_gpr(test_input) -> None:
    gpr1 = GPR(test_input)
    assert not gpr1.body
    assert len(gpr1.genes) == 0
    gpr1.update_genes()
    assert len(gpr1.genes) == 0
    assert gpr1.to_string() == ""
    assert gpr1.eval()
    assert len(gpr1.genes) == 0


def test_one_gene_gpr():
    gpr1 = GPR.from_string("a")
    assert len(gpr1.genes) == 1
    gpr1.update_genes()
    assert len(gpr1.genes) == 1
    assert gpr1.genes == {"a"}
    assert gpr1.to_string() == "a"
    assert gpr1.eval()
    assert len(gpr1.genes) == 1


# Gets an iterable of all combinations of genes except the empty list. Used to
# evaluate AND gprs
def powerset_ne(iterable):
    "powerset_ne([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(
        itertools.combinations(s, r) for r in range(1, len(s) + 1)
    )


@pytest.mark.parametrize(
    "gpr_input, num_genes, gpr_genes, gpr_output_string",
    [
        ("a & b", 2, {"a", "b"}, "a and b"),
        ("a and b", 2, {"a", "b"}, "a and b"),
        pytest.param(
            "a AND b", 2, {"a", "b"}, "a and b", marks=pytest.mark.filterwarnings
        ),
        ("a & b & c", 3, {"a", "b", "c"}, "(a and b) and c"),
        ("a and b and c", 3, {"a", "b", "c"}, "a and b and c"),
    ],
)
def test_and_gpr(gpr_input, num_genes, gpr_genes, gpr_output_string):
    gpr1 = GPR.from_string(gpr_input)
    assert len(gpr1.genes) == num_genes
    gpr1.update_genes()
    assert len(gpr1.genes) == num_genes
    assert gpr1.genes == gpr_genes
    assert gpr1.to_string() == gpr_output_string
    assert gpr1.eval()
    for ko_genes in powerset_ne(gpr_genes):
        assert not gpr1.eval(ko_genes)
    assert gpr1.body
    gpr1.copy()


# Gets an iterable of all combinations of genes except a single gene and the empty
# list. Used to evaluate OR gprs
def all_except_one(iterable):
    "all_except_one([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(
        itertools.combinations(s, r) for r in range(1, len(s))
    )


@pytest.mark.parametrize(
    "gpr_input, num_genes, gpr_genes, gpr_output_string",
    [
        ("a | b", 2, {"a", "b"}, "a or b"),
        ("a or b", 2, {"a", "b"}, "a or b"),
        pytest.param(
            "a OR b", 2, {"a", "b"}, "a or b", marks=pytest.mark.filterwarnings
        ),
        ("a | b | c", 3, {"a", "b", "c"}, "(a or b) or c"),
        ("a or b or c", 3, {"a", "b", "c"}, "a or b or c"),
    ],
)
def test_or_gpr(gpr_input, num_genes, gpr_genes, gpr_output_string):
    gpr1 = GPR.from_string(gpr_input)
    assert len(gpr1.genes) == num_genes
    gpr1.update_genes()
    assert len(gpr1.genes) == num_genes
    assert gpr1.genes == gpr_genes
    assert gpr1.to_string() == gpr_output_string
    assert gpr1.eval()
    for ko_genes in all_except_one(gpr_genes):
        assert gpr1.eval(ko_genes)
    assert not gpr1.eval(gpr_genes)
    assert gpr1.body
    gpr1.copy()


@pytest.mark.parametrize(
    "gpr_input",
    [
        "(a | b) & c",
        "(a or b) and c",
        pytest.param("(a OR b) AND c", marks=pytest.mark.filterwarnings),
    ],
)
def test_complicated_gpr(gpr_input):
    gpr1 = GPR.from_string(gpr_input)
    assert len(gpr1.genes) == 3
    gpr1.update_genes()
    assert len(gpr1.genes) == 3
    assert gpr1.to_string() == "(a or b) and c"
    assert gpr1.eval()
    assert gpr1.body
    assert gpr1.eval()
    assert gpr1.eval("a")
    assert gpr1.eval("b")
    assert not gpr1.eval("c")
    assert not gpr1.eval(["a", "b"])
    assert not gpr1.eval(["a", "b", "c"])
    gpr1.copy()


@pytest.mark.parametrize(
    "string_to_ast, num_genes, gpr_genes, gpr_output_string",
    [
        ("a", 1, {"a"}, "a"),
        ("a | b", 2, {"a", "b"}, "a or b"),
        ("a or b", 2, {"a", "b"}, "a or b"),
        pytest.param("a OR b", 2, {"a", "b"}, "a or b", marks=pytest.mark.xfail),
    ],
)
def test_gpr_from_ast_or(
    string_to_ast, num_genes, gpr_genes, gpr_output_string
) -> None:
    ast_tree = ast_parse(string_to_ast, "<string>", "eval")
    gpr1 = GPR(ast_tree)
    assert len(gpr1.genes) == num_genes
    gpr1.update_genes()
    assert len(gpr1.genes) == num_genes
    assert gpr1.genes == gpr_genes
    assert gpr1.to_string() == gpr_output_string
    assert gpr1.eval()
    for ko_genes in all_except_one(gpr_genes):
        assert gpr1.eval(ko_genes)
    assert not gpr1.eval(gpr_genes)
    gpr1.copy()


@pytest.mark.parametrize(
    "string_to_ast, num_genes, gpr_genes, gpr_output_string",
    [
        ("a & b", 2, {"a", "b"}, "a and b"),
        ("a and b", 2, {"a", "b"}, "a and b"),
        pytest.param("a AND b", 2, {"a", "b"}, "a and b", marks=pytest.mark.xfail),
    ],
)
def test_gpr_from_ast_and(
    string_to_ast, num_genes, gpr_genes, gpr_output_string
) -> None:
    ast_tree = ast_parse(string_to_ast, "<string>", "eval")
    gpr1 = GPR(ast_tree)
    assert len(gpr1.genes) == num_genes
    gpr1.update_genes()
    assert len(gpr1.genes) == num_genes
    assert gpr1.genes == gpr_genes
    assert gpr1.to_string() == gpr_output_string
    assert gpr1.eval()
    for ko_genes in powerset_ne(gpr_genes):
        assert not gpr1.eval(ko_genes)
    gpr1.copy()


@pytest.mark.parametrize("test_input", [["a", "b"], {"a", "b"}])
def test_wrong_input_gpr_error(test_input):
    with pytest.raises(TypeError):
        GPR.from_string(test_input)
    with pytest.raises(TypeError):
        GPR(test_input)


@pytest.mark.parametrize("test_input", ["a |", "a &"])
def test_wrong_input_gpr_warning(test_input):
    with pytest.warns(SyntaxWarning):
        gpr1 = GPR.from_string(test_input)
        assert gpr1.body is None
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
