import re
from warnings import warn

from .Species import Species


class Gene(Species):

    """A Gene is a special class of metabolite.


    TODO: Make design decisions about TUs and such
    """

    def __init__(self, id, formula=None,
                 name=None, compartment=None, strand='+',
                 locus_start=0, locus_end=0, functional=True):
        """
        id: A string.

        formula: cobra.Formula or a chemical formula str.  Defaults to None
        to save time in pickling and such.

        name: String.  A human readable name.

        compartment: None or a dictionary indicating the cellular location
        of the metabolite.  Used when in a cobra.Reaction or Model
        object

        strand: '+' or '-' to indicate forward or reverse strand for DNA.

        locus_start: Int.  The index of the starting base for the gene.

        locus_end: Int. The index of the last base for the gene.

        functional: Boolean.  Indicate whether the gene is functional.  If it
        is not functional then it cannot be used in an enzyme complex nor
        can its products be used.

        """
        Species.__init__(self, id, formula=formula,
                         name=name, compartment=compartment)
        self.locus_start = locus_start
        self.locus_end = locus_end
        self.strand = strand
        self.functional = functional

    def remove_from_model(self, model=None,
                          make_dependent_reactions_nonfunctional=True):
        """Removes the association

        make_dependent_reactions_nonfunctional: Boolean.  If True then replace
        the gene with 'False' in the gene association, else replace the gene
        with 'True'

        .. note :: Simulating gene knockouts is much better handled by
                   cobra.manipulation.delete_model_genes

        """
        if model is not None:
            warn("passing the model in is unnecessary and deprecated")
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
