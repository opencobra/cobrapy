"""Provide functions to modify model components."""

from ast import NodeTransformer
from functools import partial
from itertools import chain
from typing import TYPE_CHECKING, Dict

from cobra.util import get_context


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
    gene_renamer = _GeneEscaper()
    for rxn in model.reactions:
        if rxn.gpr is not None:
            gene_renamer.visit(rxn.gpr)
    model.repair()


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

    # That's not right
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
    remove_genes = set()
    context = get_context(model)

    # Needs to be added first since the history is executed from the tail and
    # this has to be run last to repair the Gene <-> Reaction mapping
    if context:
        for rxn in model.reactions:
            context(partial(rxn.update_genes_from_gpr))

    for old_name, new_name in rename_dict.items():
        # undefined if there a value matches a different key
        try:
            gene_index = model.genes.index(old_name)
        except ValueError:
            continue
        new_gene_present = new_name in model.genes
        if new_gene_present:
            old_gene = model.genes.get_by_id(old_name)
            # Added in case not renaming some genes:
            if old_gene is not model.genes.get_by_id(new_name):
                remove_genes.add(old_gene)
                recompute_reactions.update(old_gene.reactions)
        else:
            # rename old gene to new gene
            gene = model.genes[gene_index]
            gene.id = new_name
            model.genes._generate_index()
            recompute_reactions.update(gene.reactions)
            if context:
                context(model.genes._generate_index)
                context(partial(setattr, gene, "id", old_name))

    gene_renamer = _Renamer(rename_dict)
    for rxn in recompute_reactions:
        if rxn.gpr is not None:
            old_gpr = rxn.gpr.copy()
            gene_renamer.visit(rxn.gpr)
            if context:
                context(partial(setattr, rxn, "_gpr", old_gpr))

    model.repair()

    for i in remove_genes:
        model.genes.remove(i)
        i._model = None
        if context:
            context(partial(model.genes.add, i))
            context(partial(setattr, i, "_model", model))
