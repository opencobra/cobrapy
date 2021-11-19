# -*- coding: utf-8 -*-

from __future__ import absolute_import

import re
from ast import (
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
from typing import Set, Union
from warnings import warn

from _ast import AST
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


def eval_gpr(expr, knockouts):
    """evaluate compiled ast of gene_reaction_rule with knockouts

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
    .. deprecated ::
    Use GPR().eval() in the future. Because of the GPR() class,
    this function will be removed.
    """
    warn(
        "eval_gpr() will be removed soon." "Use GPR().eval(knockouts) in the future",
        DeprecationWarning,
    )
    if isinstance(expr, Expression) or isinstance(expr, GPR):
        return eval_gpr(expr.body, knockouts)
    elif isinstance(expr, Name):
        return expr.id not in knockouts
    elif isinstance(expr, BoolOp):
        op = expr.op
        if isinstance(op, Or):
            return any(eval_gpr(i, knockouts) for i in expr.values)
        elif isinstance(op, And):
            return all(eval_gpr(i, knockouts) for i in expr.values)
        else:
            raise TypeError("unsupported operation " + op.__class__.__name__)
    elif expr is None:
        return True
    else:
        raise TypeError("unsupported operation  " + repr(expr))


class GPRCleaner(NodeTransformer):
    """Parses compiled ast of a gene_reaction_rule and identifies genes

    Parts of the tree are rewritten to allow periods in gene ID's and
    bitwise boolean operations"""

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


def parse_gpr(str_expr):
    """parse gpr into AST

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
        "parse_gpr() will be removed soon." "Use GPR(string_gpr=str_expr) in the future",
        DeprecationWarning,
    )
    str_expr = str_expr.strip()
    if len(str_expr) == 0:
        return None, set()
    for char, escaped in replacements:
        if char in str_expr:
            str_expr = str_expr.replace(char, escaped)
    escaped_str = keyword_re.sub("__cobra_escape__", str_expr)
    escaped_str = number_start_re.sub("__cobra_escape__", escaped_str)
    tree = ast_parse(escaped_str, "<string>", "eval")
    cleaner = GPRCleaner()
    cleaner.visit(tree)
    eval_gpr(tree, set())  # ensure the rule can be evaluated
    return tree, cleaner.gene_set


class Gene(Species):
    """A Gene in a cobra model

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
        Species.__init__(self, id=id, name=name)
        self._functional = functional

    @property
    def functional(self):
        """A flag indicating if the gene is functional.

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
        # if value:
        #     self._gpr_gene.is_gene_functional = S.true
        # else:
        #     self._gpr_gene.is_gene_functional = S.false

    def knock_out(self):
        """Knockout gene by marking it as non-functional and setting all
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
        """Removes the association

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
    """A Gene Reaction rule in a cobra model, using AST as base class

    Parameters
    ----------
    string_gpr : string
        A GPR in string format, which will be parsed into AST
    """

    def __init__(self, string_gpr=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._genes = set()
        if string_gpr:
            self.from_string(string_gpr)

    def from_string(self, string_gpr) -> None:
        """

        Parameters
        ----------
        string_gpr: str
            a string that describes the gene rules, in a format like
            A & B

        Returns
        -------
        Does not return

            Sets self.body as
            Parsed AST tree that has the gene rules
            This function also sets self._genes with the gene ids in the AST

        """
        uppercase_AND = re.compile(r"\bAND\b")
        uppercase_OR = re.compile(r"\bOR\b")
        str_expr = string_gpr.strip()
        if len(str_expr) == 0:
            self.body = []
            return
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
                    "uppercase AND/OR found in rule '%s' for '%s'"
                    % (string_gpr, repr(self))
                )
                string_gpr = uppercase_AND.sub("and", string_gpr)
                string_gpr = uppercase_OR.sub("or", string_gpr)
            try:
                tree = ast_parse(string_gpr, "<string>", "eval")
            except SyntaxError as e:
                warn(
                    "malformed gene_reaction_rule '%s' for %s"
                    % (string_gpr, repr(self))
                )
                warn("GPR will be empty")
                warn(e.msg)
                return None
        cleaner = GPRCleaner()
        cleaner.visit(tree)
        self._genes = cleaner.gene_set
        if "" in self._genes:
            self._genes.remove("")
        self.body = tree.body
        self.eval()  # ensure the rule can be evaluated

    @property
    def geneset(self):
        self.updategenes()
        return frozenset(self._genes)

    def updategenes(self):
        if hasattr(self, "body"):
            cleaner = GPRCleaner()
            cleaner.visit(self.body)
            self._genes = cleaner.gene_set
            if "" in self._genes:
                self._genes.remove("")

    # @staticmethod
    # def from_sbml(sbml_exp):

    def __eval_gpr(self, expr, knockouts):
        """evaluate compiled ast of gene_reaction_rule with knockouts

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
        # just always call the recursions as self.__eval_gpr(a, b)
        if isinstance(expr, Expression) or isinstance(expr, GPR):
            return self.__eval_gpr(expr.body, knockouts)
        elif isinstance(expr, Name):
            return expr.id not in knockouts
        elif isinstance(expr, BoolOp):
            op = expr.op
            if isinstance(op, Or):
                return any(self.__eval_gpr(i, knockouts) for i in expr.values)
            elif isinstance(op, And):
                return all(self.__eval_gpr(i, knockouts) for i in expr.values)
            else:
                raise TypeError("unsupported operation " + op.__class__.__name__)
        elif expr is None:
            return True
        else:
            raise TypeError("unsupported operation  " + repr(expr))

    def eval(self, knockouts: Union[DictList, Set] = None):
        if knockouts is None:
            knockouts = set()
        return self.__eval_gpr(self.body, knockouts=knockouts)

    def __ast2str(self,
                  expr: Union[Expression, BoolOp, Name, list], level: int = 0,
                  names: dict = None
                  ) -> str:
        """convert compiled ast to gene_reaction_rule str

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
        if isinstance(expr, Expression) | isinstance(expr, GPR):
            return self.__ast2str(expr.body, 0, names) if hasattr(expr, "body") else ""
        elif isinstance(expr, Name):
            return names.get(expr.id, expr.id) if names else expr.id
        elif isinstance(expr, BoolOp):
            op = expr.op
            if isinstance(op, Or):
                str_exp = " or ".join(
                    self.__ast2str(i, level + 1, names) for i in expr.values)
            elif isinstance(op, And):
                str_exp = " and ".join(
                    self.__ast2str(i, level + 1, names) for i in expr.values)
            else:
                raise TypeError("unsupported operation " + op.__class__.__name)
            return "(" + str_exp + ")" if level else str_exp
        elif expr is None or (isinstance(expr, list) and len(expr) == 0):
            return ""
        else:
            raise TypeError("unsupported operation  " + repr(expr))
        pass

    def to_string(self, names: dict = None) -> str:
        """convert compiled ast to gene_reaction_rule str

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
        """
        return self.__ast2str(self, names=names)

    def __str__(self):
        """convert compiled ast to gene_reaction_rule str

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

    def __repr_html__(self):
        return """
        <table>
            <tr>
                <td><strong>GPR</strong></td><td>{gpr}</td>
            </tr>
        </table>
        """.format(
            gpr=format_long_string(self.gene_reaction_rule, 100)
        )

    def copy(self):
        """Copy a GPR"""
        return deepcopy(self)

    def __repr__(self):
        return "%s.%s(%r)" % (
            self.__class__.__module__,
            self.__class__.__qualname__,
            self.to_string(),
        )

    # def as_symbolic(self):
    #     # ...


# functions for gene reaction rules
def ast2str(
    expr: Union[Expression, GPR, BoolOp, Name, list], level: int = 0, names: dict = None
) -> str:
    """convert compiled ast to gene_reaction_rule str

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
    .. deprecated ::
    Use GPR.eval(knockouts=) in the future. Because of the GPR() class,
    this function will be removed.
    """
    warn(
        "eval_gpr() will be removed soon." "Use GPR().eval(knockouts=) in the future",
        DeprecationWarning,
    )
    if isinstance(expr, Expression) | isinstance(expr, GPR):
        return ast2str(expr.body, 0, names) if hasattr(expr, "body") else ""
    elif isinstance(expr, Name):
        return names.get(expr.id, expr.id) if names else expr.id
    elif isinstance(expr, BoolOp):
        op = expr.op
        if isinstance(op, Or):
            str_exp = " or ".join(ast2str(i, level + 1, names) for i in expr.values)
        elif isinstance(op, And):
            str_exp = " and ".join(ast2str(i, level + 1, names) for i in expr.values)
        else:
            raise TypeError("unsupported operation " + op.__class__.__name)
        return "(" + str_exp + ")" if level else str_exp
    elif expr is None or (isinstance(expr, list) and len(expr) == 0):
        return ""
    else:
        raise TypeError("unsupported operation  " + repr(expr))
