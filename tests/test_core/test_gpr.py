"""Test functions of cobra.core.gene.GPR ."""

import itertools
from ast import parse as ast_parse
from typing import Dict, Iterable, Iterator, List, Set, Tuple, Union

import pytest
from sympy.core.expr import Expr
from sympy.core.symbol import Symbol
from sympy.logic import And, Or
from sympy.logic.boolalg import BooleanFunction

from cobra.core.gene import GPR, ast2str, eval_gpr, parse_gpr
from cobra.core.model import Model


def test_gpr() -> None:
    """Test GPR instance creation and basic usage."""
    gpr1 = GPR()
    assert len(gpr1.genes) == 0
    gpr1.update_genes()
    assert len(gpr1.genes) == 0
    assert gpr1.to_string() == ""
    assert gpr1.eval()
    gpr2 = gpr1.copy()
    assert isinstance(gpr2, GPR)
    assert len(gpr1.genes) == 0


@pytest.mark.parametrize("test_input", ["", "", None])
def test_empty_gpr(test_input) -> None:
    """Test empty GPR."""
    gpr1 = GPR(test_input)
    assert not gpr1.body
    assert len(gpr1.genes) == 0
    gpr1.update_genes()
    assert len(gpr1.genes) == 0
    assert gpr1.to_string() == ""
    assert gpr1.eval()
    assert len(gpr1.genes) == 0


def test_one_gene_gpr() -> None:
    """Test single gene GPR."""
    gpr1 = GPR.from_string("a")
    assert len(gpr1.genes) == 1
    gpr1.update_genes()
    assert len(gpr1.genes) == 1
    assert gpr1.genes == {"a"}
    assert gpr1.to_string() == "a"
    assert gpr1.eval()
    assert len(gpr1.genes) == 1


def powerset_ne(iterable: Iterable[str]) -> Iterator[Tuple[str, ...]]:
    """Get all combinations of an iterable except the empty list.

    Gets an iterable of all combinations of genes except the empty list.
    Used to evaluate AND gprs
    powerset_ne([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
    """
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
def test_and_gpr(gpr_input, num_genes, gpr_genes, gpr_output_string) -> None:
    """Test 'and' GPR."""
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


def all_except_one(iterable: Iterable[str]) -> Iterator[Tuple[str, ...]]:
    """Generate all combinations from an iterable, while leaving one out.

    Gets an iterable of all combinations of genes except the complete list and the
    empty list. Used to evaluate OR gprs
    all_except_one([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3)
    """
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
def test_or_gpr(
    gpr_input: str, num_genes: int, gpr_genes: Set, gpr_output_string: str
) -> None:
    """Test 'or' GPR."""
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
def test_complicated_gpr(gpr_input: str) -> None:
    """Test complicated GPR."""
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
    string_to_ast: str, num_genes: int, gpr_genes: set, gpr_output_string: str
) -> None:
    """Test GPR from AST 'or'."""
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
    string_to_ast: str, num_genes: int, gpr_genes: set, gpr_output_string: str
) -> None:
    """Test GPR from AST 'and'."""
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
def test_wrong_input_gpr_error(test_input: Union[list, set]) -> None:
    """Test error for incorrect GPR input."""
    with pytest.raises(TypeError):
        GPR.from_string(test_input)
    with pytest.raises(TypeError):
        GPR(test_input)


@pytest.mark.parametrize("test_input", ["a |", "a &", "a and ()", "a or ()"])
def test_wrong_input_gpr_warning(test_input: str) -> None:
    """Test warning for incorrect GPR input."""
    with pytest.warns(SyntaxWarning):
        gpr1 = GPR.from_string(test_input)
        assert gpr1.body is None
        assert len(gpr1.genes) == 0


def test_gpr_that_needs_two_replacements() -> None:
    """Test GPR with multi replacements."""
    gpr1 = GPR.from_string(
        "(591001.3.peg.1891 AND 591001.3.peg.1892 AND 591001.3.peg.1893)"
    )
    assert len(gpr1.genes) == 3
    assert "591001.3.peg.1891" in gpr1.genes
    assert "591001.3.peg.1892" in gpr1.genes
    assert "591001.3.peg.1893" in gpr1.genes


def test_deprecated_gpr() -> None:
    """Test deprecated GPR."""
    gpr1 = GPR.from_string("(a | b) & c")
    with pytest.deprecated_call():
        assert ast2str(gpr1) == "(a or b) and c"
    with pytest.deprecated_call():
        assert eval_gpr(gpr1, set())
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
        assert eval_gpr(gpr1, set())
    with pytest.deprecated_call():
        assert eval_gpr(gpr1, {"a"})
    with pytest.deprecated_call():
        assert eval_gpr(gpr1, {"b"})
    with pytest.deprecated_call():
        assert not eval_gpr(gpr1, {"c"})


def test_gpr_as_symbolic() -> None:
    """Test GPR as symbolic expression."""
    gpr1 = GPR()
    assert gpr1.as_symbolic() == Symbol("")
    gpr1 = GPR.from_string("")
    assert gpr1.as_symbolic() == Symbol("")
    gpr1 = GPR.from_string("a")
    assert isinstance(gpr1.as_symbolic(), Symbol)
    assert gpr1.as_symbolic() == Symbol("a")
    gpr2 = GPR.from_string("a & b")
    assert isinstance(gpr2.as_symbolic(), BooleanFunction)
    assert gpr2.as_symbolic() == And(Symbol("a"), Symbol("b"))
    assert gpr1 != gpr2


@pytest.mark.parametrize(
    "gpr_input, symbolic_gpr",
    [
        ("a | b", Or(Symbol("a"), Symbol("b"))),
        ("a or b", Or(Symbol("a"), Symbol("b"))),
        pytest.param(
            "a OR b", Or(Symbol("a"), Symbol("b")), marks=pytest.mark.filterwarnings
        ),
        ("a | b", Or(Symbol("b"), Symbol("a"))),
        ("a | b | c", Or(Symbol("a"), Symbol("b"), Symbol("c"))),
        ("a or b or c", Or(Symbol("a"), Symbol("b"), Symbol("c"))),
        ("(a OR b) AND c", And(Symbol("c"), Or(Symbol("a"), Symbol("b")))),
    ],
)
def test_gpr_as_symbolic_boolean(gpr_input: str, symbolic_gpr: Expr) -> None:
    """Test GPR as symbolic boolean expression."""
    gpr1 = GPR().from_string(gpr_input)
    assert isinstance(gpr1.as_symbolic(), BooleanFunction)
    assert gpr1.as_symbolic() == symbolic_gpr
    if "OR" not in gpr_input and "AND" not in gpr_input:
        ast_tree = ast_parse(gpr_input, "<string>", "eval")
        gpr1 = GPR(ast_tree)
        assert isinstance(gpr1.as_symbolic(), BooleanFunction)
        assert gpr1.as_symbolic() == symbolic_gpr


def test_gpr_equality() -> None:
    """Test GPR equality."""
    assert GPR() == GPR()
    assert GPR() == GPR.from_string("")
    assert GPR() != GPR.from_string("a")
    assert GPR.from_string("a") == GPR.from_string("a")


gpr_str_lists = {
    "a_and_b_strs": ["a & b", "a and b", "a AND b", "b & a", "b and a", "b AND a"],
    "a_or_b_strs": ["a | b", "a or b", "a OR b", "b | a", "b or a", "b OR a"],
    "a_b_c_or_strs": [
        "a | b | c",
        "a or b or c",
        "a OR b or c",
        "(a or b) or c",
        "b | a | c",
        "b or a or c",
        "b | c | a",
        "b or c | a",
        "b or c or a",
        "c or a or b",
        "c | a or b",
        "c or a | b",
        "c or b or a",
        "c | b or a",
        "c or b | a",
        "c | b | a",
    ],
    "a_b_c_and_strs": [
        "a & b & c",
        "a and b and c",
        "a AND b and c",
        "(a and b) and c",
        "b & a & c",
        "b and a and c",
        "b & c & a",
        "b and c & a",
        "b and c and a",
        "c and a and b",
        "c & a and b",
        "c and a & b",
        "c and b and a",
        "c & b and a",
        "c and b & a",
        "c & b & a",
    ],
    "a_b_c_or_and_strs": [
        "(a OR b) AND c",
        "(a | b) & c",
        "(a & c) | b & c",
        "c & (a | b)",
    ],
}


@pytest.fixture(params=list(gpr_str_lists.keys()))
def gpr_list(request: pytest.FixtureRequest) -> List[str]:
    """Provide fixture for GPR list."""
    gpr_str_list = gpr_str_lists[request.param]
    return gpr_str_list


def test_gpr_equality_with_bolean_logic(gpr_list: List[str]) -> None:
    """Test GPR equality with boolean logic."""
    for i in range(len(gpr_list)):
        for j in range(i + 1, len(gpr_list)):
            assert GPR().from_string(gpr_list[i]) == GPR.from_string(gpr_list[j])


@pytest.fixture(params=list(itertools.combinations(gpr_str_lists.keys(), 2)))
def gpr_lists(request: pytest.FixtureRequest) -> Dict[str, List[str]]:
    """Provide fixture for GPR dictionary."""
    gpr_dict = dict()
    gpr_dict["gpr1"] = gpr_str_lists[request.param[0]]
    gpr_dict["gpr2"] = gpr_str_lists[request.param[1]]
    return gpr_dict


def test_gpr_inequality_boolean(gpr_lists) -> None:
    """Test GPR inequality for boolean expression."""
    gpr_list1 = gpr_lists["gpr1"]
    gpr_list2 = gpr_lists["gpr2"]
    for i in range(len(gpr_list1)):
        for j in range(len(gpr_list2)):
            assert GPR.from_string(gpr_list1[i]) != GPR.from_string(gpr_list2[j])


def test_gpr_symbolism_benchmark(large_model: Model, benchmark) -> None:
    """Benchmark as symbolic time."""
    model = large_model.copy()

    def gpr_symbolic():
        for i in range(len(model.reactions)):
            rxn1 = model.reactions[i]
            gpr1 = rxn1.gpr
            gpr1.as_symbolic()

    benchmark(gpr_symbolic)


def test_gpr_equality_benchmark(model: Model, benchmark) -> None:
    """Benchmark equality of GPR using the mini model."""
    model2 = model.copy()

    def gpr_equality_all_reactions():
        for i in range(len(model.reactions)):
            rxn1 = model.reactions[i]
            rxn2 = model2.reactions[i]
            assert rxn1.gpr == rxn2.gpr

    benchmark(gpr_equality_all_reactions)


@pytest.mark.parametrize(
    "gpr_input, symbolic_gpr",
    [
        ("", Symbol("")),
        ("a", Symbol("a")),
        ("a & b", And(Symbol("a"), Symbol("b"))),
        ("a | b", Or(Symbol("a"), Symbol("b"))),
        ("a or b", Or(Symbol("a"), Symbol("b"))),
        pytest.param(
            "a OR b", Or(Symbol("a"), Symbol("b")), marks=pytest.mark.filterwarnings
        ),
        ("a | b", Or(Symbol("b"), Symbol("a"))),
        ("a | b | c", Or(Symbol("a"), Symbol("b"), Symbol("c"))),
        ("a or b or c", Or(Symbol("a"), Symbol("b"), Symbol("c"))),
        ("(a OR b) AND c", And(Symbol("c"), Or(Symbol("a"), Symbol("b")))),
    ],
)
def test_gpr_from_symbolic(gpr_input: str, symbolic_gpr: Expr) -> None:
    """Test GPR creation from symbolic and from string is equal."""
    gpr1 = GPR().from_symbolic(symbolic_gpr)
    gpr2 = GPR().from_string(gpr_input)
    assert gpr1 == gpr2


def test_gpr_from_as_symbolic_equality(large_model) -> None:
    """Test as_symbolic followed by from_symbolic gives a GPR equivalent to original."""
    model = large_model.copy()

    for i in range(len(model.reactions)):
        rxn1 = model.reactions[i]
        gpr1 = rxn1.gpr
        gpr2 = GPR().from_symbolic(gpr1.as_symbolic())
        assert gpr1 == gpr2
