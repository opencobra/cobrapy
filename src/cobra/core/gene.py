# -*- coding: utf-8 -*-

from __future__ import absolute_import

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
from typing import FrozenSet, Iterable, Set, Tuple, Union
from warnings import warn

from cobra.core.dictlist import DictList
from cobra.core.species import Species
from cobra.util import resettable
from cobra.util.util import format_long_string


keywords = list(kwlist)
keywords.remove("and")
keywords.remove("or")
keywords.extend(("True", "False"))
keyword_re = re.compile(r"(?=\b(%s)\b)" % "|".join(keywords))
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

    def __init__(self):
        NodeVisitor.__init__(self)
        self.gene_set = set()

    def visit_Name(self, node) -> None:
        self.gene_set.add(node.id)

    def visit_BoolOp(self, node: BoolOp) -> None:
        self.generic_visit(node)
        for val in node.values:
            self.visit(val)


class GPRCleaner(NodeTransformer):
    """Parses compiled ast of a gene_reaction_rule and identifies genes.

    Parts of the tree are rewritten to allow periods in gene ID's and
    bitwise boolean operations
    """

    def __init__(self):
        NodeTransformer.__init__(self)
        self.gene_set = set()

    def visit_Name(self, node):
        if node.id.startswith("__cobra_escape__"):
            node.id = node.id[16:]
        for char, escaped in replacements:
            if escaped in node.id:
                node.id = node.id.replace(escaped, char)
        self.gene_set.add(node.id)
        return node

    def visit_BinOp(self, node):
        self.generic_visit(node)
        if isinstance(node.op, BitAnd):
            return BoolOp(And(), (node.left, node.right))
        elif isinstance(node.op, BitOr):
            return BoolOp(Or(), (node.left, node.right))
        else:
            raise TypeError("unsupported operation '%s'" % node.op.__class__.__name__)


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

    def __init__(self, id=None, name="", functional=True):
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
        Species.__init__(self, id=id, name=name)
        self._functional = functional

    @property
    def functional(self):
        """Flag indicating if the gene is functional.

        Changing the flag is reverted upon exit if executed within the model
        as context.
        """
        return self._functional

    @functional.setter
    @resettable
    def functional(self, value):
        if not isinstance(value, bool):
            raise ValueError("expected boolean")
        self._functional = value

    def knock_out(self):
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

    def remove_from_model(
        self, model=None, make_dependent_reactions_nonfunctional=True
    ):
        """Removes the association.

        Parameters
        ----------
        model : cobra model
           The model to remove the gene from
        make_dependent_reactions_nonfunctional : bool
           If True then replace the gene with 'False' in the gene
           association, else replace the gene with 'True'


        .. deprecated :: 0.4
            Use cobra.manipulation.delete_model_genes to simulate knockouts
            and cobra.manipulation.remove_genes to remove genes from
            the model.

        """
        warn("Use cobra.manipulation.remove_genes instead")
        if model is not None:
            if model != self._model:
                raise Exception(
                    "%s is a member of %s, not %s"
                    % (repr(self), repr(self._model), repr(model))
                )
        if self._model is None:
            raise Exception("%s is not in a model" % repr(self))

        if make_dependent_reactions_nonfunctional:
            gene_state = "False"
        else:
            gene_state = "True"
        the_gene_re = re.compile("(^|(?<=( |\()))%s(?=( |\)|$))" % re.escape(self.id))

        # remove reference to the gene in all groups
        associated_groups = self._model.get_associated_groups(self)
        for group in associated_groups:
            group.remove_members(self)

        self._model.genes.remove(self)
        self._model = None

        for the_reaction in list(self._reaction):
            the_reaction._gene_reaction_rule = the_gene_re.sub(
                gene_state, the_reaction.gene_reaction_rule
            )
            the_reaction._genes.remove(self)
            # Now, deactivate the reaction if its gene association evaluates
            # to False
            the_gene_reaction_relation = the_reaction.gene_reaction_rule
            for other_gene in the_reaction._genes:
                other_gene_re = re.compile(
                    "(^|(?<=( |\()))%s(?=( |\)|$))" % re.escape(other_gene.id)
                )
                the_gene_reaction_relation = other_gene_re.sub(
                    "True", the_gene_reaction_relation
                )

            if not eval(the_gene_reaction_relation):
                the_reaction.lower_bound = 0
                the_reaction.upper_bound = 0
        self._reaction.clear()

    def _repr_html_(self):
        return """
        <table>
            <tr>
                <td><strong>Gene identifier</strong></td><td>{id}</td>
            </tr><tr>
                <td><strong>Name</strong></td><td>{name}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{address}</td>
            </tr><tr>
                <td><strong>Functional</strong></td><td>{functional}</td>
            </tr><tr>
                <td><strong>In {n_reactions} reaction(s)</strong></td><td>
                    {reactions}</td>
            </tr>
        </table>""".format(
            id=self.id,
            name=self.name,
            functional=self.functional,
            address="0x0%x" % id(self),
            n_reactions=len(self.reactions),
            reactions=format_long_string(", ".join(r.id for r in self.reactions), 200),
        )


class GPR(Module):
    """A Gene Reaction rule in a cobra model, using AST as base class.

    Parameters
    ----------
    gpr_from : Expression or Module or AST
        A GPR in AST format
    """

    def __init__(self, gpr_from: Union[Expression, Module, AST] = None, **kwargs):
        super().__init__(**kwargs)
        self._genes = set()
        self.body = None
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
                # noinspection PyTypeChecker
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

        try:
            tree = ast_parse(escaped_str, "<string>", "eval")
        except (SyntaxError, TypeError) as e:
            if "AND" in string_gpr or "OR" in string_gpr:
                warn(
                    f"Uppercase AND/OR found in rule '{string_gpr}'.",
                    SyntaxWarning,
                )
                string_gpr = uppercase_AND.sub("and", string_gpr)
                string_gpr = uppercase_OR.sub("or", string_gpr)
            try:
                tree = ast_parse(string_gpr, "<string>", "eval")
            except SyntaxError as e:
                warn(
                    f"Malformed gene_reaction_rule '{string_gpr}' for {repr(gpr)}",
                    SyntaxWarning,
                )
                warn("GPR will be empty")
                warn(e.msg)
                return gpr
        return cls(tree)

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
        expr: Union[Expression, BoolOp, Name, list],
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
        Calls __aststr()
        """
        # noinspection PyTypeChecker
        return self._ast2str(self, names=names)

    def copy(self):
        """Copy a GPR."""
        return deepcopy(self)

    def __repr__(self):
        return "%s.%s(%r)" % (
            self.__class__.__module__,
            self.__class__.__qualname__,
            self.to_string(),
        )

    def __str__(self):
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

    def _repr_html__(self):
        return """<p><strong>GPR</strong></p><p>{gpr}</p>""".format(
            gpr=format_long_string(self.to_string(), 100)
        )

    # def as_symbolic(self):
    #     # ...


def eval_gpr(expr, knockouts):
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
