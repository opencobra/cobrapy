"""Define base Object class in Cobra."""

from typing import Optional


class Object:
    """Defines common behavior of object in cobra.core."""

    def __init__(self, id: Optional[str] = None, name: str = "") -> None:
        """Initialize a simple object with an identifier.

        Parameters
        ----------
        id: string, optional
            the identifier to associate with the object
        name: string, optional
            The name to associate with the object. Default "".

        Objects will have notes and _annotation as dicitionaries, initialized as empty
        dictionaries.
        """
        self._id = id
        self.name = name

        self.notes = {}
        self._annotation = {}

    @property
    def id(self) -> str:
        """Get the Object id.

        Returns
        -------
        id: str
        """
        return getattr(self, "_id", None)

    @id.setter
    def id(self, value) -> None:
        """Set the id to value.

        Parameters
        ----------
        value: str
            The string to set the id to.

        Raises
        ------
        TypeError if value is not a string.
        """
        if value == self.id:
            pass
        elif not isinstance(value, str):
            raise TypeError("ID must be a string")
        elif getattr(self, "_model", None) is not None:
            self._set_id_with_model(value)
        else:
            self._id = value

    def _set_id_with_model(self, value) -> None:
        """Set id with model.

        This appears to be a stub so it can be modified in dependant classes.

        Parameters
        ----------
        value: str
            The string to set the id to.
        """
        self._id = value

    @property
    def annotation(self) -> dict:
        """Get annotation dictionary.

        Returns
        -------
        _annotation: dict
            Returns _annotation as a dictionary.
        """
        return self._annotation

    @annotation.setter
    def annotation(self, annotation):
        """Set annotation.

        Parameters
        ----------
        annotation: dict
            Annotation dictionary to set _annotation to. Will raise error if not dict.

        Raises
        ------
        TypeError if annotation not a dict.
        """
        if not isinstance(annotation, dict):
            raise TypeError("Annotation must be a dict")
        else:
            self._annotation = annotation

    def __getstate__(self) -> dict:
        """Get state of annotation.

        To prevent excessive replication during deepcopy, ignores _model in state.

        Returns
        -------
        state: dict
            Dictionary of state, excluding _model.
        """
        state = self.__dict__.copy()
        if "_model" in state:
            state["_model"] = None
        return state

    def __repr__(self) -> str:
        """Return string representation of Object, with class.

        Returns
        -------
        str
            Composed of class.name, id and hexadecimal of id.
        """
        return f"<{self.__class__.__name__} {self.id} at {id(self):#x}>"

    def __str__(self) -> str:
        """Return string representation of object.

        Returns
        -------
        str
            Object.id as string.
        """
        return str(self.id)
