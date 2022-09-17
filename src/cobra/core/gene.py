"""Provide functions for dealing with genes and gene product rules (GPR)."""

import logging
import re
from ast import (
    AST,
    And,
    BitAnd,
    BitOr,
    BoolOp,
    Expression,
    Module,
    Name,
    NodeTransformer,
    NodeVisitor,
    Or,
)
from ast import parse as ast_parse
from copy import deepcopy
from keyword import kwlist
from typing import FrozenSet, Iterable, Optional, Set, Tuple, Union
from warnings import warn

import sympy.logic.boolalg as spl
from sympy import Symbol

from ..util import resettable
from ..util.util import format_long_string
from .dictlist import DictList
from .species import Species


# TODO - When https://github.com/symengine/symengine.py/issues/334 is resolved,
#  change sympy.Symbol (above in imports) to optlang.symbolics.Symbol

logger = logging.getLogger(__name__)

keywords = list(kwlist)
keywords.remove("and")
keywords.remove("or")
keywords.extend(("True", "False"))
keyword_re = re.compile(rf"(?=\b({'|'.join(keywords)})\b)")
number_start_re = re.compile(r"(?=\b[0-9])")

replacements = (
    (".", "__COBRA_DOT__"),
    ("'", "__COBRA_SQUOTE__"),
    ('"', "__COBRA_DQUOTE__"),
    (":", "__COBRA_COLON__"),
    ("/", "__COBRA_FSLASH__"),
    ("\\", "__COBRA_BSLASH"),
    ("-", "__COBRA_DASH__"),
    ("=", "__COBRA_EQ__"),
)


class GPRWalker(NodeVisitor):
    """Identifies genes in an AST/GPR tree.

    Walks over the tree, and identifies the id of each Name node
    """

    def __init__(self, **kwargs) -> None:
        """Initialize a new object.

        Other Parameters
        ----------------
        **kwargs:
            Further keyword arguments are passed on to the parent class.

        """
        super().__init__(**kwargs)
        self.gene_set = set()

    def visit_Name(self, node: Name) -> None:
        """Visit a Gene node and add the id to the gene_set.

        Parameters
        ----------
        node: ast.Name
            The node to visit

        """
        self.gene_set.add(node.id)

    def visit_BoolOp(self, node: BoolOp) -> None:
        """Visit a BoolOp node (AND/OR) and visit the children to add them to gene_set.

        Parameters
        ----------
        node: ast.BoolOp
            The node to visit

        """
        self.generic_visit(node)
        for val in node.values:
            self.visit(val)


class GPRCleaner(NodeTransformer):
    """Parses compiled ast of a gene_reaction_rule and identifies genes.

    Parts of the tree are rewritten to allow periods in gene ID's and
    bitwise boolean operations
    """

    def __init__(self, **kwargs) -> None:
        """Initialize a new object.

        Other Parameters
        ----------------
        **kwargs:
            Further keyword arguments are passed on to the parent class.

        """
        super().__init__(**kwargs)
        self.gene_set = set()

    def visit_Name(self, node: Name) -> Name:
        """Visit a Gene node and add the id to the gene_set.

        The gene id will be cleaned used __cobra_escape__ and replacements
        dictionary (see above).

        Parameters
        ----------
        node: ast.Name
            The node to visit

        Returns
        -------
        node: ast.Name
            The transformed node (with the id changed).

        """
        if node.id.startswith("__cobra_escape__"):
            node.id = node.id[16:]
        for char, escaped in replacements:
            if escaped in node.id:
                node.id = node.id.replace(escaped, char)
        self.gene_set.add(node.id)
        return node

    def visit_BinOp(self, node: BoolOp) -> None:
        """Visit a BoolOp node (AND/OR) and visit the children (genes) to process them.

        Parameters
        ----------
        node: ast.BoolOp
            The node to visit. Nodes other than And() and Or() will cause an error.

        Returns
        -------
        node: ast.BoolOp
            The node with the children transformed.
        """
        self.generic_visit(node)
        if isinstance(node.op, BitAnd):
            return BoolOp(And(), (node.left, node.right))
        elif isinstance(node.op, BitOr):
            return BoolOp(Or(), (node.left, node.right))
        else:
            raise TypeError(f"unsupported operation '{node.op.__class__.__name__}'")


def parse_gpr(str_expr: str) -> Tuple:
    """Parse GPR into AST.

    Parameters
    ----------
    str_expr : string
        string with the gene reaction rule to parse

    Returns
    -------
    tuple
        elements ast_tree and gene_ids as a set

    .. deprecated ::
    Use GPR(string_gpr=str_expr) in the future. Because of the GPR() class,
    this function will be removed.
    """
    warn(
        "parse_gpr() will be removed soon."
        "Use GPR(string_gpr=str_expr) in the future",
        DeprecationWarning,
    )
    gpr_tree = GPR.from_string(str_expr)
    return gpr_tree, gpr_tree.genes


class Gene(Species):
    """A Gene in a cobra model.

    Parameters
    ----------
    id : string
        The identifier to associate the gene with
    name: string
        A longer human readable name for the gene
    functional: bool
        Indicates whether the gene is functional.  If it is not functional
        then it cannot be used in an enzyme complex nor can its products be
        used.
    """

    # noinspection PyShadowingBuiltins
    def __init__(self, id: str = None, name: str = "", functional: bool = True) -> None:
        """Initialize a gene.

        Parameters
        ----------
        id: str
            A string that will identify the gene.
        name: str
            A (longer) string that will identify the gene. Can have more special
            characters.
        functional: bool
            A flag whether or not the gene is functional
        """
        super().__init__(id=id, name=name)
        self._functional = functional

    @property
    def functional(self) -> bool:
        """Flag indicating if the gene is functional.

        Changing the flag is reverted upon exit if executed within the model
        as context.
        """
        return self._functional

    @functional.setter
    @resettable
    def functional(self, value: bool) -> None:
        if not isinstance(value, bool):
            raise ValueError("expected boolean")
        self._functional = value

    def knock_out(self) -> None:
        """Knockout gene by marking it as non-functional.

        Knockout gene by marking it as non-functional and setting all
        associated reactions bounds to zero.
        The change is reverted upon exit if executed within the model as
        context.
        """
        self.functional = False
        for reaction in self.reactions:
            if not reaction.functional:
                reaction.bounds = (0, 0)

    def _repr_html_(self):
        return f"""
        <table>
            <tr>
                <td><strong>Gene identifier</strong></td><td>{self.id}</td>
            </tr><tr>
                <td><strong>Name</strong></td><td>{self.name}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{id(self):#x}</td>
            </tr><tr>
                <td><strong>Functional</strong></td><td>{self.functional}</td>
            </tr><tr>
                <td><strong>In {len(self.reactions)} reaction(s)</strong></td><td>
                    {format_long_string(", ".join(r.id for r in self.reactions), 200)}
                    </td>
            </tr>
        </table>"""


class GPR(Module):
    """A Gene Reaction rule in a cobra model, using AST as base class.

    Parameters
    ----------
    gpr_from : Expression or Module or AST
        A GPR in AST format
    """

    def __init__(self, gpr_from: Union[Expression, Module, AST] = None, **kwargs):
        """Initialize a gene.

        Parameters
        ----------
        gpr_from: Expression, Module, AST
            An AST expression that will be parsed to GPR.
        **kwargs:
            Further keyword arguments are passed on to the parent class.
        """
        super().__init__(**kwargs)
        self._genes = set()
        self.body: Optional[list] = None
        if gpr_from:
            if isinstance(gpr_from, str):
                self.from_string(gpr_from)
                raise TypeError(
                    f"GPR accepts AST, not string. "
                    f'Next time, use GPR().from_string("{gpr_from}")'
                )
            elif isinstance(gpr_from, (Expression, Module)):
                cleaner = GPRCleaner()
                cleaner.visit(gpr_from)
                self._genes = deepcopy(cleaner.gene_set)
                self.body = deepcopy(gpr_from.body)
                self.eval()
            else:
                raise TypeError("GPR requires AST Expression or Module")

    @classmethod
    def from_string(cls, string_gpr: str) -> "GPR":
        """Construct a GPR from a string.

        Parameters
        ----------
        string_gpr: str
            a string that describes the gene rules, in a format like
            A & B

        Returns
        -------
        GPR:
            returns a new GPR while setting  self.body as
            Parsed AST tree that has the gene rules
            This function also sets self._genes with the gene ids in the AST

        """
        if not isinstance(string_gpr, str):
            raise TypeError(
                f"{cls.__name__}.from_string "
                f"requires a str argument, not {type(string_gpr)}."
            )
        gpr = cls()
        uppercase_AND = re.compile(r"\bAND\b")
        uppercase_OR = re.compile(r"\bOR\b")
        str_expr = string_gpr.strip()
        if len(str_expr) == 0:
            gpr.body = None
            return gpr
        for char, escaped in replacements:
            if char in str_expr:
                str_expr = str_expr.replace(char, escaped)
        escaped_str = keyword_re.sub("__cobra_escape__", str_expr)
        escaped_str = number_start_re.sub("__cobra_escape__", escaped_str)
        # Some mat models have () in gr_rules which leads to a complicated error later
        escaped_str = escaped_str.replace("()", "")

        try:
            tree = ast_parse(escaped_str, "<string>", "eval")
        except (SyntaxError, TypeError) as e:
            if "AND" in escaped_str or "OR" in escaped_str:
                # noinspection PyTypeChecker
                logger.warning(
                    f"Uppercase AND/OR found in rule '{string_gpr}'.",
                )
                logger.warning(e.msg)
                warn(
                    "Uppercase AND/OR found in rule '{}'.".format(string_gpr),
                    SyntaxWarning,
                )
                escaped_str = uppercase_AND.sub("and", escaped_str)
                escaped_str = uppercase_OR.sub("or", escaped_str)
            try:
                tree = ast_parse(escaped_str, "<string>", "eval")
            except SyntaxError:
                # noinspection PyTypeChecker
                logger.warning(
                    f"Malformed gene_reaction_rule '{escaped_str}' for {string_gpr}",
                    exc_info=1,
                )
                logger.warning("GPR will be empty")
                warn(
                    "Malformed gene_reaction_rule '{}'".format(escaped_str),
                    SyntaxWarning,
                )
                return gpr
        gpr = cls(tree)
        gpr.update_genes()
        return gpr

    @property
    def genes(self) -> FrozenSet:
        """To check the genes.

        This property updates the genes before returning them, in case the GPR was
        changed and the genes weren't.

        Returns
        -------
        genes: frozenset
            All the genes in a frozen set. Do not try to change them with this property.
        """
        self.update_genes()
        return frozenset(self._genes)

    def update_genes(self) -> None:
        """Update genes, used after changes in GPR.

        Walks along the AST tree of the GPR class, and modifies self._genes

        """
        if self.body:
            walker = GPRWalker()
            walker.visit(self)
            self._genes = deepcopy(walker.gene_set)

    def _eval_gpr(
        self,
        expr: Union[Expression, list, BoolOp, Name],
        knockouts: Union[DictList, set],
    ) -> bool:
        """Evaluate compiled ast of gene_reaction_rule with knockouts.

        Parameters
        ----------
        expr : Expression or GPR or list or BoolOp or Name
            The ast of the gene reaction rule
        knockouts : DictList, set
            Set of genes that are knocked out

        Returns
        -------
        bool
            True if the gene reaction rule is true with the given knockouts
            otherwise false
        """
        # just always call the recursions as self._eval_gpr(a, b)
        if isinstance(expr, (Expression, GPR)):
            if not expr.body:
                return True
            return self._eval_gpr(expr.body, knockouts)
        elif isinstance(expr, Name):
            return expr.id not in knockouts
        elif isinstance(expr, BoolOp):
            op = expr.op
            if isinstance(op, Or):
                # noinspection PyTypeChecker
                return any(self._eval_gpr(i, knockouts) for i in expr.values)
            elif isinstance(op, And):
                # noinspection PyTypeChecker
                return all(self._eval_gpr(i, knockouts) for i in expr.values)
            else:
                raise TypeError(f"Unsupported operation: {op.__class__.__name__}")
        elif expr is None:
            return True
        else:
            raise TypeError(f"Unsupported operation: {repr(expr)}")

    def eval(self, knockouts: Union[DictList, Set, str, Iterable] = None) -> bool:
        """Evaluate compiled ast of gene_reaction_rule with knockouts.

        This function calls _eval_gpr, but allows more flexibility in input, including
        name, and list.

        Parameters
        ----------
        knockouts
            Which gene or genes to knoc out

        Returns
        -------
        bool
            True if the gene reaction rule is true with the given knockouts
            otherwise false

        """
        if knockouts is None:
            knockouts = set()
        if knockouts is str:
            knockouts = list(knockouts)
        if self.body:
            return self._eval_gpr(self.body, knockouts=knockouts)
        else:
            return True

    def _ast2str(
        self,
        expr: Union["GPR", Expression, BoolOp, Name, list],
        level: int = 0,
        names: dict = None,
    ) -> str:
        """Convert compiled ast to gene_reaction_rule str.

        Parameters
        ----------
        expr : AST or GPR or list or Name or BoolOp
            string for a gene reaction rule, e.g "a and b"
        level : int
            internal use only
        names : dict
            Dict where each element id a gene identifier and the value is the
            gene name. Use this to get a rule str which uses names instead. This
            should be done for display purposes only. All gene_reaction_rule
            strings which are computed with should use the id.

        Returns
        ------
        string
            The gene reaction rule
        """
        if isinstance(expr, (Expression, GPR)):
            return self._ast2str(expr.body, 0, names) if expr.body else ""
        elif isinstance(expr, Name):
            return names.get(expr.id, expr.id) if names else expr.id
        elif isinstance(expr, BoolOp):
            op = expr.op
            if isinstance(op, Or):
                # noinspection PyTypeChecker
                str_exp = " or ".join(
                    self._ast2str(i, level + 1, names) for i in expr.values
                )
            elif isinstance(op, And):
                # noinspection PyTypeChecker
                str_exp = " and ".join(
                    self._ast2str(i, level + 1, names) for i in expr.values
                )
            else:
                # noinspection PyTypeChecker
                raise TypeError(f"Unsupported operation: {op.__class__.__name}")
            return f"({str_exp})" if level else str_exp
        elif expr is None or (isinstance(expr, list) and len(expr) == 0):
            return ""
        else:
            raise TypeError(f"Unsupported operation: {repr(expr)}")

    def to_string(self, names: dict = None) -> str:
        """Convert compiled ast to gene_reaction_rule str.

        Parameters
        ----------
        self : GPR
           compiled ast Module describing GPR
        names: dict
           dictionary of gene ids to gene names. If this is empty, returns gene ids

        Returns
        ------
        string
            The gene reaction rule

        Notes
        -----
        Calls _aststr()
        """
        return self._ast2str(self, names=names)

    def copy(self):
        """Copy a GPR."""
        return deepcopy(self)

    def __copy__(self) -> "GPR":
        """Ensure a correct shallow copy."""
        return self.copy()

    def __repr__(self) -> str:
        """Return the GPR with module, class, and code to recreate it."""
        return (
            f"{self.__class__.__module__}.{self.__class__.__qualname__}"
            f"({self.to_string()!r})"
        )

    def __str__(self) -> str:
        """Convert compiled ast to gene_reaction_rule str.

        Parameters
        ----------
        self : GPR
            compiled ast Module describing GPR

        Returns
        ------
        string
            The gene reaction rule
        """
        return self.to_string(names={})

    def _repr_html_(self) -> str:
        return f"""<p><strong>GPR</strong></p><p>{format_long_string(self.to_string(),
                                                                     100)}</p>"""

    def as_symbolic(
        self,
        names: dict = None,
    ) -> Union[spl.Or, spl.And, Symbol]:
        """Convert compiled ast to sympy expression.

        Parameters
        ----------
        self : GPR
           compiled ast Module describing GPR
        names: dict
           dictionary of gene ids to gene names. If this is empty,
           returns sympy expression using gene ids

        Returns
        ------
        Symbol or BooleanFunction
            SYMPY expression (Symbol or And or Or). Symbol("") if the GPR is empty

        Notes
        -----
        Calls _symbolic_gpr()
        """
        # noinspection PyTypeChecker
        if names:
            GPRGene_dict = {gid: Symbol(names[gid]) for gid in self.genes}
        else:
            GPRGene_dict = None
        return self._symbolic_gpr(self, GPRGene_dict=GPRGene_dict)

    def _symbolic_gpr(
        self,
        expr: Union["GPR", Expression, BoolOp, Name, list] = None,
        GPRGene_dict: dict = None,
    ) -> Union[spl.Or, spl.And, Symbol]:
        """Parse gpr into SYMPY using ast similar to _ast2str().

        Parameters
        ----------
        expr : AST or GPR or list or Name or BoolOp
            compiled GPR
        GPRGene_dict: dict
            dictionary from gene id to GPRGeneSymbol
        Returns
        -------
        Symbol or BooleanFunction
            SYMPY expression (Symbol or And or Or). Symbol("") if the GPR is empty
        """
        if GPRGene_dict is None:
            GPRGene_dict = {gid: Symbol(name=gid) for gid in expr.genes}
        if isinstance(expr, (Expression, GPR)):
            return (
                self._symbolic_gpr(expr.body, GPRGene_dict) if expr.body else Symbol("")
            )
        else:
            if isinstance(expr, Name):
                return GPRGene_dict.get(expr.id)
            elif isinstance(expr, BoolOp):
                op = expr.op
                if isinstance(op, Or):
                    # noinspection PyTypeChecker
                    sym_exp = spl.Or(
                        *[self._symbolic_gpr(i, GPRGene_dict) for i in expr.values]
                    )
                elif isinstance(op, And):
                    # noinspection PyTypeChecker
                    sym_exp = spl.And(
                        *[self._symbolic_gpr(i, GPRGene_dict) for i in expr.values]
                    )
                else:
                    raise TypeError("Unsupported operation " + op.__class__.__name)
                return sym_exp
            elif not expr:
                return Symbol("")
            else:
                raise TypeError("Unsupported Expression  " + repr(expr))

    @classmethod
    def from_symbolic(cls, sympy_gpr: Union[spl.BooleanFunction, Symbol]) -> "GPR":
        """Construct a GPR from a sympy expression.

        Parameters
        ----------
        sympy_gpr: sympy
            a sympy that describes the gene rules, being a Symbol for single genes
            or a BooleanFunction for AND/OR relationships

        Returns
        -------
        GPR:
            returns a new GPR while setting  self.body as
            Parsed AST tree that has the gene rules
            This function also sets self._genes with the gene ids in the AST

        """

        def _sympy_to_ast(
            sympy_expr: Union[spl.BooleanFunction, Symbol]
        ) -> Union[BoolOp, Name]:
            if sympy_expr.func is spl.Or:
                return BoolOp(
                    op=Or(), values=[_sympy_to_ast(i) for i in sympy_expr.args]
                )
            elif sympy_expr.func is spl.And:
                return BoolOp(
                    op=And(), values=[_sympy_to_ast(i) for i in sympy_expr.args]
                )
            elif not sympy_expr.args:
                return Name(id=sympy_expr.name)
            else:
                raise TypeError(f"Unsupported operation: {sympy_expr.func}")

        if not isinstance(sympy_gpr, (spl.BooleanFunction, Symbol)):
            raise TypeError(
                f"{cls.__name__}.from_symbolic "
                f"requires a sympy BooleanFunction or "
                f"Symbol argument, not {type(sympy_gpr)}."
            )
        gpr = cls()
        if sympy_gpr == Symbol(""):
            gpr.body = None
            return gpr
        try:
            tree = Expression(_sympy_to_ast(sympy_gpr))
        except SyntaxError as e:
            logger.warning(
                f"Problem with sympy expression '{sympy_gpr}' for {repr(gpr)}",
            )
            logger.warning("GPR will be empty")
            logger.warning(e.msg)
            return gpr
        gpr = cls(tree)
        gpr.update_genes()
        return gpr

    def __eq__(self, other) -> bool:
        """Check equality of GPR via symbolic equality."""
        if not self.body and not other.body:
            return True
        elif not self.body or not other.body:
            return False
        else:
            self_symb = self.as_symbolic()
            other_symb = other.as_symbolic()
            if isinstance(self_symb, Symbol) and isinstance(other_symb, Symbol):
                return self_symb == other_symb
            if isinstance(self_symb, Symbol) or isinstance(other_symb, Symbol):
                return False
            return self_symb.equals(other_symb)


def eval_gpr(expr: Union[Expression, GPR], knockouts: Union[DictList, set]) -> bool:
    """Evaluate compiled ast of gene_reaction_rule with knockouts.

    .. deprecated ::
    Use GPR().eval() in the future. Because of the GPR() class,
    this function will be removed.

    Parameters
    ----------
    expr : Expression or GPR
        The ast of the gene reaction rule
    knockouts : DictList, set
        Set of genes that are knocked out

    Returns
    -------
    bool
        True if the gene reaction rule is true with the given knockouts
        otherwise false
    """
    warn(
        "eval_gpr() will be removed soon." "Use GPR().eval(knockouts) in the future",
        DeprecationWarning,
    )
    if isinstance(expr, GPR):
        return expr.eval(knockouts=knockouts)
    else:
        return GPR(expr).eval(knockouts=knockouts)


# functions for gene reaction rules
def ast2str(expr: Union[Expression, GPR], level: int = 0, names: dict = None) -> str:
    """Convert compiled ast to gene_reaction_rule str.

    Parameters
    ----------
    expr : AST or GPR
        AST or GPR
    level : int
        internal use only. Ignored because of GPR() class, kept only for interface
        consistency with code still using ast2str.
    names : dict
        Dict where each element id a gene identifier and the value is the
        gene name. Use this to get a rule str which uses names instead. This
        should be done for display purposes only. All gene_reaction_rule
        strings which are computed with should use the id.

    Returns
    ------
    string
        The gene reaction rule
    .. deprecated ::
    Use GPR.to_string(names=) in the future. Because of the GPR() class,
    this function will be removed.
    """
    warn(
        "ast2satr() will be removed soon. Use gpr.to_string(names=names) in the future",
        DeprecationWarning,
    )
    if isinstance(expr, GPR):
        return expr.to_string(names=names)
    else:
        return GPR(expr).to_string(names=names)
