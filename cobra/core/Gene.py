import re
from warnings import warn

from .Species import Species


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
