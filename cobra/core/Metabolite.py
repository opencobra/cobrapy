from warnings import warn

from .Species import Species


class Metabolite(Species):
    """Metabolite is a class for holding information regarding
    a metabolite in a cobra.Reaction object.

    """

    def __init__(self, id=None, formula=None,
                 name=None, compartment=None):
        """
        id: str

        formula: cobra.Formula or String  of a chemical formula.

        name: str
            A human readable name.

        compartment: None or a dictionary indicating the cellular location
        of the metabolite.  Used when in a cobra.Reaction or Model
        object

        """
        Species.__init__(self, id, formula, name, compartment)
        self._constraint_sense = 'E'
        self._bound = 0.

    @property
    def y(self):
        """The shadow price for the metabolite in the most recent solution

        Shadow prices are computed from the dual values of the bounds in
        the solution.

        """
        try:
            return self._model.solution.y_dict[self.id]
        except Exception as e:
            if self._model is None:
                raise Exception("not part of a model")
            if not hasattr(self._model, "solution") or \
                    self._model.solution is None or \
                    self._model.solution.status == "NA":
                raise Exception("model has not been solved")
            if self._model.solution.status != "optimal":
                raise Exception("model solution was not optimal")
            raise e  # Not sure what the exact problem was

    def remove_from_model(self, method='subtractive', **kwargs):
        """Removes the association from self.model

        method: 'subtractive' or 'destructive'.
            If 'subtractive' then the metabolite is removed from all
            associated reactions.  If 'destructive' then all associated
            reactions are removed from the Model.

        """
        # why is model being taken in as a parameter? This plays
        # back to the question of allowing a Metabolite to be associated
        # with multiple Models
        if "model" in kwargs:
            warn("model argument deprecated")

        self._model.metabolites.remove(self)
        self._model = None
        if method.lower() == 'subtractive':
            for the_reaction in list(self._reaction):
                the_coefficient = the_reaction._metabolites[self]
                the_reaction.subtract_metabolites({self: the_coefficient})
        elif method.lower() == 'destructive':
            for x in self._reaction():
                x.remove_from_model()
        else:
            raise Exception(method + " is not 'subtractive' or 'destructive'")
