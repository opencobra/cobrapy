# -*- coding: utf-8 -*-

from __future__ import absolute_import

import re
from _ast import AST
from ast import And, BitAnd, BitOr, BoolOp, Expression, Name, NodeTransformer, Or, \
    NodeVisitor, Module
from ast import parse as ast_parse
from keyword import kwlist
from typing import Union, Set
from warnings import warn
from copy import deepcopy
# When https://github.com/symengine/symengine.py/issues/334 is resolved, change it to
# optlang.symbolics.Symbol
from sympy import Symbol as Symbol
from sympy import SympifyError, simplify_logic, symbols
from sympy.parsing.sympy_parser import parse_expr as parse_expr_sympy, \
    standard_transformations
from sympy.logic.boolalg import Or as sp_Or, And as sp_And
from sympy.core.singleton import S

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


def sympy2str(expr, names=None):
    """convert compiled sympy to gene_reaction_rule str
    Parameters
    ----------
    expr : sympy Expr
        compiled sympy expression for a gene rule
    names : dict
        Dict where each element id a gene identifier and the value is the
        gene alias. Use this to get a rule str which uses aliases instead. This
        should be done for display purposes only. All gene_reaction_rule
        strings which are computed with should use the id.
    Returns
    ------
    string
        The gene reaction rule
    """
    if names is not None:
        aliases = {id: names.get(id, id) for id in names.keys()}
        expr = expr.xreplace(aliases)
    if expr == S.true or expr == S.false:
        str_exp = ""
    else:
        str_exp = str(expr)
        str_exp = str_exp.replace('&', 'and')
        str_exp = str_exp.replace('|', 'or')
    return str_exp


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
    """
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


def eval_gpr_sympy(expr, knockouts=None):
    """evaluate compiled sympy of gene_reaction_rule with knockouts
    Parameters
    ----------
    expr : Sympy Expression
        The compiled sympy expression of the gene reaction rule
    knockouts : DictList, set
        Set of genes that are knocked out
    Returns
    -------
    bool
        True if the gene reaction rule is true with the given knockouts
        otherwise false
    """

    if knockouts is None:
        knockouts = []
    gene_knockouts = {i: S.true for i in expr.free_symbols}
    for gene in knockouts:
        gene_knockouts[gene] = S.false
    return bool(expr.subs(gene_knockouts))


def eval_gpr_sympy_xreplace(expr, knockouts=None):
    """evaluate compiled sympy of gene_reaction_rule with knockouts
    Parameters
    ----------
    expr : Sympy Expression
        The compiled sympy expression of the gene reaction rule
    knockouts : DictList, set
        Set of genes that are knocked out
    Returns
    -------
    bool
        True if the gene reaction rule is true with the given knockouts
        otherwise false
    """

    # Does eval_gpr keep track of the status of the genes
    # rxn1.gpr = 'A & B'
    # model.genes.get_by_id('A').functional=False
    # eval_gpr(rxn1) ?
    if knockouts is None:
        knockouts = []
    gene_knockouts = {gene: gene not in knockouts for gene in expr.free_symbols}
    expr.xreplace(gene_knockouts)
    return bool(expr)


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


class UpdateGenes(NodeVisitor):
    def __init__(self):
        NodeVisitor.__init__(self)
        self.gene_set = set()

    def visit_Name(self, node):
        self.gene_set.add(node.id)


# Using unescaped ":", "," or " " causes sympy to make a range or separate items
str_for_sympy_symbols_dict = {':': r'\:', ',': r'\,', ' ': r'\ '}


def fix_str_for_symbols(unescaped_str):
    escaped_str = unescaped_str
    for key, value in str_for_sympy_symbols_dict.items():
        escaped_str = escaped_str.replace(key, value)
    return escaped_str


class GPRSympifier(NodeVisitor):
    """Parses compiled ast of a gene_reaction_rule to sympy and identifies genes
    """

    def __init__(self, sym_dict: dict):
        NodeVisitor.__init__(self)
        self.gene_set = set()
        self.sym_dict = sym_dict

    def visit_Expression(self, node):
        return self.visit(node.body) if hasattr(node, "body") else None

    def visit_Name(self, node):
        if node.id.startswith("__cobra_escape__"):
            node.id = node.id[16:]
        for char, escaped in replacements:
            if escaped in node.id:
                node.id = node.id.replace(escaped, char)
        self.gene_set.add(node.id)
        return self.sym_dict[node.id]

    def visit_BinOp(self, node):
        if isinstance(node.op, BitAnd):
            return sp_And(*[node.left, node.right], evaluate=False)
        elif isinstance(node.op, BitOr):
            return sp_Or(*[node.left, node.right], evaluate=False)
        else:
            raise TypeError("unsupported operation '%s'" % node.op.__class__.__name__)

    def visit_BoolOp(self, node):
        if isinstance(node.op, And):
            return sp_And(*[self.visit(i) for i in node.values], evaluate=False)
        if isinstance(node.op, Or):
            return sp_Or(*[self.visit(i) for i in node.values], evaluate=False)


def parse_gpr_sympy(str_expr, GPRGene_dict):
    """parse gpr into SYMPY using ast and Node Visitor
    Parameters
    ----------
    str_expr : string
        string with the gene reaction rule to parse
    GPRGene_dict: dict
        dictionary from gene id to GPRGeneSymbol
    Returns
    -------
    tuple
        elements SYMPY expression and gene_ids as a set
    """
    str_expr = str_expr.strip()
    if len(str_expr) == 0:
        return None, set()
    for char, escaped in replacements:
        if char in str_expr:
            str_expr = str_expr.replace(char, escaped)
    escaped_str = keyword_re.sub("__cobra_escape__", str_expr)
    escaped_str = number_start_re.sub("__cobra_escape__", escaped_str)
    tree = ast_parse(escaped_str, "<string>", "eval")
    sympifier = GPRSympifier(GPRGene_dict)
    sympy_exp = sympifier.visit(tree)
    return sympy_exp, sympifier.gene_set


class GPRSympifierToString(NodeVisitor):
    """Parses compiled ast of a gene_reaction_rule to sympy and identifies genes
    """

    def __init__(self):
        NodeVisitor.__init__(self)
        self.gene_set = set()

    def visit_Expression(self, node):
        return self.visit(node.body) if hasattr(node, "body") else None

    def visit_Name(self, node):
        if node.id.startswith("__cobra_escape__"):
            node.id = node.id[16:]
        for char, escaped in replacements:
            if escaped in node.id:
                node.id = node.id.replace(escaped, char)
        self.gene_set.add(node.id)
        return str(node.id)

    def visit_BoolOp(self, node):
        _children = [self.visit(i) for i in node.values]
        if isinstance(node.op, And):
            return 'And(' + ", ".join(_children) + ')'
        if isinstance(node.op, Or):
            return 'Or(' + ", ".join(_children) + ')'


def parse_gpr_sympy_str(str_expr, GPRGene_dict=None):
    """parse gpr into SYMPY using ast and Node Visitor
    Parameters
    ----------
    str_expr : string
        string with the gene reaction rule to parse
    GPRGene_dict: dict
        dictionary from gene id to GPRGeneSymbol
    Returns
    -------
    tuple
        elements SYMPY expression and gene_ids as a set
    """
    str_expr = str_expr.strip()
    if len(str_expr) == 0:
        return None, set()
    if GPRGene_dict is None:
        GPRGene_dict = dict()
    for char, escaped in replacements:
        if char in str_expr:
            str_expr = str_expr.replace(char, escaped)
    escaped_str = keyword_re.sub("__cobra_escape__", str_expr)
    escaped_str = number_start_re.sub("__cobra_escape__", escaped_str)
    tree = ast_parse(escaped_str, "<string>", "eval")
    sympifier = GPRSympifierToString()
    sympy_exp = sympifier.visit(tree)
    GPRGene_dict.update({i: Symbol(fix_str_for_symbols(i)) for i in sympifier.gene_set})
    return parse_expr_sympy(sympy_exp, local_dict=GPRGene_dict)


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
    """
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
                warn("malformed gene_reaction_rule '%s' for %s" %
                     (string_gpr, repr(self)))
                warn('GPR will be empty')
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

    @property
    def gene_reaction_rule(self):
        return str(self)

    # @staticmethod
    # def from_sbml(sbml_exp):

    def eval(self, knockouts: Union[DictList, Set] = None):
        if knockouts is None:
            knockouts = set()
        return eval_gpr(self.body, knockouts=knockouts)

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
        return ast2str(self)

    def str_names(self, names: dict = None) -> str:
        """convert compiled ast to gene_reaction_rule str

         Parameters
         ----------
         self : GPR
            compiled ast Module describing GPR
         names: dict
            dictionary of gene ids to gene names

         Returns
         ------
         string
             The gene reaction rule
         """
        return ast2str(self, names=names)

    #
    def __repr_html__(self):
        return """
        <table>
            <tr>
                <td><strong>GPR</strong></td><td>{gpr}</td>
            </tr>
        </table>
        """.format(gpr=format_long_string(self.gene_reaction_rule, 100))

    def copy(self):
        """Copy a GPR
        """
        return deepcopy(self)

    def __repr__(self):
        return '%s.%s(%r)' % (self.__class__.__module__,
                              self.__class__.__qualname__,
                              ast2str(self))

    # def as_symbolic(self):
    #     # ...


# functions for gene reaction rules
def ast2str(expr: Union[Expression, GPR, BoolOp, Name, list], level: int = 0,
            names: dict = None) -> str:
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
