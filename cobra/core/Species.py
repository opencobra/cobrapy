from warnings import warn
from copy import deepcopy
from .Formula import Formula
from .Object import Object


class Species(Object):
    """Species is a class for holding information regarding
    a chemical Species


    """

    def __init__(self, id=None, formula=None,
                 name=None, compartment=None, mnx_id=None):
        """
        id: A string.

        formula: cobra.Formula or str of a chemical formula.  Defaults to None
        to save time in pickling and such.

        name: String.  A human readable name.

        compartment: None or a dictionary indicating the cellular location
        of the metabolite.  Used when in a cobra.Reaction or Model
        object

        mnx_id: None or a String of the MetaNetX.org ID for the object.

        """
        Object.__init__(self, id, mnx_id=mnx_id)
        self.name = name
        if not name:
            self.name = self.id
        if isinstance(formula, str):
            formula = Formula(formula)

        self.formula = formula
        self.parse_composition()
        # because in a Model a metabolite may participate in multiple Reactions
        self.compartment = compartment
        # self.model is None or refers to the cobra.Model that
        # contains self
        self._model = self.charge = None
        # references to reactions that operate on this species
        self._reaction = set()

    @property
    def reactions(self):
        return frozenset(self._reaction)

    def parse_composition(self):
        """Breaks the chemical formula down by element.
        Useful for making sure Reactions are balanced.'

        """
        if isinstance(self.formula, Formula):
            self.formula.parse_composition()
        elif isinstance(self.formula, str):
            self.formula = Formula(self.formula)

    def __getstate__(self):
        """Remove the references to container reactions when serializing to
        avoid problems associated with recursion.

        """
        state = Object.__getstate__(self)
        state['_reaction'] = set()
        return state

    def copy(self):
        """When copying a reaction, it is necessary to deepcopy the
        components so the list references aren't carried over.

        Additionally, a copy of a reaction is no longer in a cobra.Model.

        This should be fixed with self.__deecopy__ if possible
        """
        return deepcopy(self)

    def guided_copy(self, the_model):
        """.. deprecated :: 0.3 Use copy directly"""
        warn("deprecated")
        the_copy = Object.guided_copy(self)
        # Copy the more complex objects in a faster fashion
        the_copy.formula = deepcopy(self.formula)
        the_copy._model = the_model
        the_copy._reaction = set()
        return(the_copy)

    def get_reaction(self):
        """Returns a list of Reactions that contain this Species"""
        warn("deprecated, used species.reactions instead")
        return list(self._reaction)

    def get_model(self):
        """Returns the Model object that contain this Object"""
        print("get_model is deprecated. used model property instead")
        return self._model

    @property
    def model(self):
        return(self._model)
