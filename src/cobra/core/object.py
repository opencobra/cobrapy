"""Define base Object class in Cobra."""

from typing import TYPE_CHECKING, Iterable, Optional, Union, Tuple


if TYPE_CHECKING:
    from cobra.core.metadata import (
        CustomAnnotation,
        MetaData,
        StandardizedAnnotation,
        Identifier,
    )


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
        from cobra.core.metadata import MetaData

        self._id = id
        self.name = name

        self.notes = {}
        self._annotations = MetaData()

    @property
    def id(self) -> Optional[str]:
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
    def annotations(self) -> "MetaData":
        """Get annotation dictionary.

        Returns
        -------
        _annotation: dict
            Returns _annotation as a dictionary.
        """
        # TODO: Fix doc
        return self._annotations

    @annotations.setter
    def annotations(self, annotations: Optional["MetaData"]):
        """Set annotations.

        Parameters
        ----------
        annotation: dict
            Annotation dictionary to set _annotation to. Will raise error if not dict.

        Raises
        ------
        TypeError if annotation not a dict.
        """
        # TODO: Fix doc
        from cobra.core.metadata import MetaData

        if annotations is None:
            self._annotations = MetaData()
        elif isinstance(annotations, MetaData):
            self._annotations = annotations
        else:
            raise TypeError(
                f"The data passed for annotation must be inside "
                f"a dictionary or MetaData: {annotations}"
            )

    def add_annotations(
        self,
        annotations: Union[
            str,
            "Identifier",
            Tuple[str, str],
            "StandardizedAnnotation",
            "CustomAnnotation",
            Iterable[
                Union[
                    str,
                    "Identifier",
                    Tuple[str, str],
                    "StandardizedAnnotation",
                    "CustomAnnotation",
                ]
            ],
        ],
    ):
        from cobra.core.metadata import (
            CustomAnnotation,
            MetaData,
            StandardizedAnnotation,
            Identifier,
        )

        if isinstance(
            annotations,
            (str, tuple, Identifier, StandardizedAnnotation, CustomAnnotation),
        ):
            annotations = [annotations]

        if self._annotations is None:
            self._annotations = MetaData()

        for annotation in annotations:
            if isinstance(annotation, (str, tuple, Identifier)):
                self._annotations.simplified.add(annotation)
            if isinstance(annotation, StandardizedAnnotation):
                self._annotations.standardized.add([annotation])
            elif isinstance(annotation, CustomAnnotation):
                self._annotations.custom.add(annotation)

    def remove_annotations(
        self,
        annotations: Union[
            "StandardizedAnnotation",
            "CustomAnnotation",
            Iterable[Union["StandardizedAnnotation", "CustomAnnotation"]],
        ],
    ):
        from cobra.core.metadata import CustomAnnotation, StandardizedAnnotation

        if isinstance(annotations, StandardizedAnnotation) or isinstance(
            annotations, CustomAnnotation
        ):
            annotations = [annotations]

        if self._annotations is None:
            raise ValueError(
                "Cannot remove annotations, because there are no annotations"
                "associated with this object."
            )

        for annotation in annotations:
            if isinstance(annotation, StandardizedAnnotation):
                self._annotations.standardized.remove(annotation)
            elif isinstance(annotation, CustomAnnotation):
                self._annotations.custom.remove(annotation)

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
