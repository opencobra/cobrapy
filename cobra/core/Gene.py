import re
from warnings import warn
from ast import parse as ast_parse, Name, And, Or, BitOr, BitAnd, BinOp, \
    BoolOp, Attribute, Expression, NodeTransformer
from keyword import kwlist

from .Species import Species


keywords = list(kwlist)
keywords.remove("and")
keywords.remove("or")
keywords.extend(("True", "False"))
keyword_re = re.compile(r"(?=\b(%s)\b)" % "|".join(keywords))


# functions for gene reaction rules
def ast2str(expr, level=0):
    """convert compiled ast to gene_reaction_rule str"""
    if isinstance(expr, Expression):
        return ast2str(expr.body, 0) if hasattr(expr, "body") else ""
    elif isinstance(expr, Name):
        return expr.id
    elif isinstance(expr, BoolOp):
        op = expr.op
        if isinstance(op, Or):
            str_exp = " or ".join(ast2str(i, level + 1) for i in expr.values)
        elif isinstance(op, And):
            str_exp = " and ".join(ast2str(i, level + 1) for i in expr.values)
        else:
            raise TypeError("unsupported operation " + op.__class__.__name)
        return "(" + str_exp + ")" if level else str_exp
    elif expr is None:
        return ""
    else:
        raise TypeError("unsupported operation  " + repr(expr))


def eval_gpr(expr, knockouts):
    """evaluate compiled ast of gene_reaction_rule with knockouts"""
    if isinstance(expr, Expression):
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
        node.id = node.id.replace("__COBRA_COLON__", ":").\
            replace("__COBRA_SQUOTE__", "'").replace("__COBRA_DQUOTE__", '"')
        self.gene_set.add(node.id)
        return node

    def visit_BinOp(self, node):
        self.generic_visit(node)
        if isinstance(node.op, BitAnd):
            return BoolOp(And(), (node.left, node.right))
        elif isinstance(node.op, BitOr):
            return BoolOp(Or(), (node.left, node.right))
        else:
            raise TypeError("unsupported operation '%s'" %
                            node.op.__class__.__name__)

    def visit_Attribute(self, node):
        # handle periods in gene id
        gene_id = node.attr
        while isinstance(node.value, Attribute):
            node = node.value
            gene_id = node.attr + "." + gene_id
        gene_id = node.value.id + "." + gene_id
        self.gene_set.add(gene_id)
        return Name(id=gene_id)


class GeneRemover(NodeTransformer):
    def __init__(self, target_genes):
        NodeTransformer.__init__(self)
        self.target_genes = {str(i) for i in target_genes}

    def visit_Name(self, node):
        return None if node.id in self.target_genes else node

    def visit_BoolOp(self, node):
        original_n = len(node.values)
        self.generic_visit(node)
        if len(node.values) == 0:
            return None
        # AND with any entities removed
        if len(node.values) < original_n and isinstance(node.op, And):
            return None
        # if one entity in an OR was removed, just that entity passed up
        if len(node.values) == 1:
            return node.values[0]
        return node


def parse_gpr(str_expr):
    """parse gpr into AST

    returns: (ast_tree, {gene_ids})"""
    if len(str_expr.strip()) == 0:
        return None, set()
    str_expr = str_expr.replace(":", "__COBRA_COLON__").\
        replace("'", "__COBRA_SQUOTE__").replace('"', "__COBRA_DQUOTE__")
    escaped_str = keyword_re.sub("__cobra_escape__", str_expr)
    tree = ast_parse(escaped_str, "<string>", "eval")
    cleaner = GPRCleaner()
    cleaner.visit(tree)
    eval_gpr(tree, set())  # ensure the rule can be evaluated
    return tree, cleaner.gene_set


class Gene(Species):

    def __init__(self, id, name=None, functional=True):
        """
        id: A string.

        name: String.  A human readable name.

        functional: Boolean.  Indicate whether the gene is functional.  If it
        is not functional then it cannot be used in an enzyme complex nor
        can its products be used.

        """
        Species.__init__(self, id, name=name)
        self.functional = functional

    def remove_from_model(self, model=None,
                          make_dependent_reactions_nonfunctional=True):
        """Removes the association

        make_dependent_reactions_nonfunctional: Boolean.  If True then replace
        the gene with 'False' in the gene association, else replace the gene
        with 'True'

        .. deprecated :: 0.4
            Use cobra.manipulation.delete_model_genes to simulate knockouts
            and cobra.manipulation.remove_genes to remove genes from
            the model.

        """
        warn("Use cobra.manipulation.remove_genes instead")
        if model is not None:
            if model != self._model:
                raise Exception("%s is a member of %s, not %s" %
                                (repr(self), repr(self._model), repr(model)))
        if self._model is None:
            raise Exception('%s is not in a model' % repr(self))

        if make_dependent_reactions_nonfunctional:
            gene_state = 'False'
        else:
            gene_state = 'True'
        the_gene_re = re.compile('(^|(?<=( |\()))%s(?=( |\)|$))' %
                                 re.escape(self.id))

        self._model.genes.remove(self)
        self._model = None

        for the_reaction in list(self._reaction):
            the_reaction._gene_reaction_rule = the_gene_re.sub(
                gene_state, the_reaction.gene_reaction_rule)
            the_reaction._genes.remove(self)
            # Now, deactivate the reaction if its gene association evaluates
            # to False
            the_gene_reaction_relation = the_reaction.gene_reaction_rule
            for other_gene in the_reaction._genes:
                other_gene_re = re.compile('(^|(?<=( |\()))%s(?=( |\)|$))' %
                                           re.escape(other_gene.id))
                the_gene_reaction_relation = other_gene_re.sub(
                    'True',
                    the_gene_reaction_relation)

            if not eval(the_gene_reaction_relation):
                the_reaction.lower_bound = 0
                the_reaction.upper_bound = 0
        self._reaction.clear()
