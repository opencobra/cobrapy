"""Provide functions to modify model components."""

from ast import NodeTransformer
from itertools import chain
from typing import TYPE_CHECKING, Dict

from ..core.gene import ast2str
from .delete import get_compiled_gene_reaction_rules


if TYPE_CHECKING:
    from cobra import Gene, Model

# Set of tuples of operators and their corresponding textual form
_renames = (
    (".", "_DOT_"),
    ("(", "_LPAREN_"),
    (")", "_RPAREN_"),
    ("-", "__"),
    ("[", "_LSQBKT"),
    ("]", "_RSQBKT"),
    (",", "_COMMA_"),
    (":", "_COLON_"),
    (">", "_GT_"),
    ("<", "_LT"),
    ("/", "_FLASH"),
    ("\\", "_BSLASH"),
    ("+", "_PLUS_"),
    ("=", "_EQ_"),
    (" ", "_SPACE_"),
    ("'", "_SQUOT_"),
    ('"', "_DQUOT_"),
)


def _escape_str_id(id_str: str) -> str:
    """Make a single string ID SBML compliant.

    Parameters
    ----------
    id_str: str
        The ID string to operate on.

    Returns
    -------
    str
        The SBML compliant ID string.

    """
    for c in ("'", '"'):
        if id_str.startswith(c) and id_str.endswith(c) and id_str.count(c) == 2:
            id_str = id_str.strip(c)
    for char, escaped_char in _renames:
        id_str = id_str.replace(char, escaped_char)
    return id_str


class _GeneEscaper(NodeTransformer):
    """Class to represent a gene ID escaper."""

    def visit_Name(self, node: "Gene") -> "Gene":
        """Escape string ID.

        Parameters
        ----------
        node: cobra.Gene
            The gene object to work on.

        Returns
        -------
        cobra.Gene
            The gene object whose ID has been escaped.

        """
        node.id = _escape_str_id(node.id)
        return node


def escape_ID(model: "Model") -> None:
    """Make all model component object IDs SBML compliant.

    Parameters
    ----------
    model: cobra.Model
        The model to operate on.

    """
    for x in chain([model], model.metabolites, model.reactions, model.genes):
        x.id = _escape_str_id(x.id)
    model.repair()
    gene_renamer = _GeneEscaper()
    for rxn, rule in get_compiled_gene_reaction_rules(model).items():
        if rule is not None:
            rxn._gene_reaction_rule = ast2str(gene_renamer.visit(rule))


class _Renamer(NodeTransformer):
    """
    Class to represent a gene renamer.

    Parameters
    ----------
    rename_dict: dict of {str: str}
        The dictionary having keys as old gene names
        and value as new gene names.

    """

    def __init__(self, rename_dict: Dict[str, str], **kwargs) -> None:
        """Initialize a new object.

        Other Parameters
        ----------------
        **kwargs:
            Further keyword arguments are passed on to the parent class.

        """
        super().__init__(**kwargs)
        self.rename_dict = rename_dict

    def visit_Name(self, node: "Gene") -> "Gene":
        """Rename a gene.

        Parameters
        ----------
        node: cobra.Gene
            The gene to rename.

        Returns
        -------
        cobra.Gene
            The renamed gene object.

        """
        node.id = self.rename_dict.get(node.id, node.id)
        return node


def rename_genes(model: "Model", rename_dict: Dict[str, str]) -> None:
    """Rename genes in a model from the `rename_dict`.

    Parameters
    ----------
    model: cobra.Model
        The model to operate on.
    rename_dict: dict of {str: str}
        The dictionary having keys as old gene names
        and value as new gene names.

    """
    recompute_reactions = set()  # need to recompute related genes
    remove_genes = []
    for old_name, new_name in rename_dict.items():
        # undefined if there a value matches a different key
        try:
            gene_index = model.genes.index(old_name)
        except ValueError:
            gene_index = None
        old_gene_present = gene_index is not None
        new_gene_present = new_name in model.genes
        if old_gene_present and new_gene_present:
            old_gene = model.genes.get_by_id(old_name)
            # Added in case not renaming some genes:
            if old_gene is not model.genes.get_by_id(new_name):
                remove_genes.append(old_gene)
                recompute_reactions.update(old_gene._reaction)
        elif old_gene_present and not new_gene_present:
            # rename old gene to new gene
            gene = model.genes[gene_index]
            # trick DictList into updating index
            model.genes._dict.pop(gene.id)  # ugh
            gene.id = new_name
            model.genes[gene_index] = gene
        elif not old_gene_present and new_gene_present:
            pass
        else:  # if not old gene_present and not new_gene_present
            # the new gene's _model will be set by repair
            # This would add genes from rename_dict
            # that are not associated with a rxn
            # cobra_model.genes.append(Gene(new_name))
            pass
    model.repair()

    gene_renamer = _Renamer(rename_dict)
    for rxn, rule in get_compiled_gene_reaction_rules(model).items():
        if rule is not None:
            rxn._gene_reaction_rule = ast2str(gene_renamer.visit(rule))

    for rxn in recompute_reactions:
        rxn.gene_reaction_rule = rxn._gene_reaction_rule
    for i in remove_genes:
        model.genes.remove(i)
