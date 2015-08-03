from warnings import warn
from copy import deepcopy
from .Object import Object


class Species(Object):
    """Species is a class for holding information regarding
    a chemical Species


    """

    def __init__(self, id=None, name=None):
        """
        id: A string.

        name: String.  A human readable name.

        """
        Object.__init__(self, id, name)
        self._model = None
        # references to reactions that operate on this species
        self._reaction = set()

    @property
    def reactions(self):
        return frozenset(self._reaction)

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

    @property
    def model(self):
        return(self._model)
