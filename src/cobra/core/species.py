"""Define the Species class, used as a base for Gene and Metabolite."""


from copy import deepcopy
from typing import TYPE_CHECKING, FrozenSet, Optional

from ..core.object import Object


if TYPE_CHECKING:
    from .. import Model


class Species(Object):
    """Species is a base class in Cobrapy.

    Species is a class for holding information regarding
    a chemical Species

    Parameters
    ----------
    id : string
       An identifier for the chemical species
    name : string
       A human readable name.
    """

    # noinspection PyShadowingBuiltins
    def __init__(
        self, id: Optional[str] = None, name: Optional[str] = None, **kwargs
    ) -> None:
        """Initialize a Species.

        Parameters
        ----------
        id : string, optional, default None
            An identifier for the chemical species
        name : string, optional, default None
            A human readable name.

        A species also contains a _model, reference to a cobra.model (initialized as
        None) and a self._reaction, a set of cobra.reactions (initialized as empty set).
        """
        super().__init__(id=id, name=name, **kwargs)
        self._model = None
        # references to reactions that operate on this species
        self._reaction = set()

    @property
    def reactions(self) -> FrozenSet:
        """Return a frozenset of reactions.

        Returns
        -------
        FrozenSet
            A frozenset that includes the reactions of the species.
        """
        return frozenset(self._reaction)

    def __getstate__(self) -> dict:
        """Return the state of the species.

        Remove the references to container reactions when serializing to
        avoid problems associated with recursion.

        Returns
        -------
        dict
            A dictionary describing the state, without the self._reaction to avoid
            recursion.
        """
        state = Object.__getstate__(self)
        state["_reaction"] = set()
        return state

    def copy(self) -> "Species":
        """Copy a species.

        When copying a reaction, it is necessary to deepcopy the
        components so the list references aren't carried over.

        Additionally, a copy of a reaction is no longer in a cobra.Model.

        This should be fixed with self.__deepcopy__ if possible

        Returns
        -------
        Species
            A copy of the species.
        """
        return deepcopy(self)

    @property
    def model(self) -> Optional["Model"]:
        """Return the model.

        Returns
        -------
        model
            Returns the cobra model that the species is associated with. None if there
            is no model associated with this species.
        """
        return self._model
