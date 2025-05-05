"""Define base Object class in Cobra."""

from typing import TYPE_CHECKING, Dict, Iterable, List, Optional, Tuple, Union


if TYPE_CHECKING:
    from cobra.core.metadata import (
        CustomAnnotation,
        Metadata,
        Resource,
        StandardizedAnnotation,
    )
    from cobra.core.metadata.standardized import SimplifiedAnnotationInterface


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

        Objects will have notes and metadata as empty attributes.
        """
        self._id = id
        self.name = name

        self.notes = {}
        self.metadata = None

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
    def metadata(self) -> "Metadata":
        """Get objects metadata.

        Returns
        -------
        Metadata
            Metadata object containing annotations and object creator and history.
        """
        from cobra.core.metadata import Metadata

        if self._metadata is None:
            self.metadata = Metadata()
        return self._metadata

    @metadata.setter
    def metadata(self, metadata: Optional["Metadata"]):
        """Set the metadata of the object.

        Parameters
        ----------
        metadata: Metadata
            Metadata object containing annotations and object history.
        """
        from cobra.core.metadata import Metadata
        from cobra.core.metadata.standardized import SimplifiedAnnotationInterface

        if metadata is None:
            self._metadata = None
            self._annotation = None
        elif isinstance(metadata, Metadata):
            self._metadata = metadata
            self._annotation = SimplifiedAnnotationInterface(self._metadata)
        else:
            raise TypeError(
                f"The data passed for annotation must be inside "
                f"a dictionary or Metadata: {metadata}"
            )

    @property
    def annotation(self) -> "SimplifiedAnnotationInterface":
        """Access standardized annotations through a dict-like interface.

        Warnings
        --------
        This attribute is in place to retain compatibility with older cobrapy versions.
        For new code, it is recommended to directly use the methods of
        StandardizedAnnotationStore, which can be accessed through
        `object.metadata.standardized`.

        Returns
        -------
        SimplifiedAnnotationInterface

        See Also
        --------
        StandardizedAnnotationStore
        SimplifiedAnnotationInterface
        """
        from cobra.core.metadata.standardized import SimplifiedAnnotationInterface

        if self._annotation is None:
            self._annotation = SimplifiedAnnotationInterface(self.metadata)
        return self._annotation

    @annotation.setter
    def annotation(
        self,
        value: Optional[
            Union[
                Dict,
                "SimplifiedAnnotationInterface",
                List[Union[str, "Resource", Tuple[str, Union[str, List[str]]]]],
            ]
        ],
    ):
        """Set the standardized annotations using a dict-like object.

        This method removes all standardized annotations and adds the ones provided as
        argument.

        Warnings
        --------
        This attribute is in place to retain compatibility with older cobrapy versions.
        For new code, it is recommended to directly use the methods of
        StandardizedAnnotationStore, which can be accessed through
        `object.metadata.standardized`.

        Parameters
        ----------
        value: dict, SimplifiedAnnotationInterface, list of str, Resource or tuples.
            Sets resources as annotations, using default qualifiers. Tuples are
            interpreted as namespace-identifiers pairs, strings should be valid URIs and
            if a dictionary is provided, its keys should represent namespaces and its
            values identifiers.

        See Also
        --------
        SimplifiedAnnotationInterface.clear
        SimplifiedAnnotationInterface.add
        """
        self.annotation.clear()
        if value is not None:
            self.annotation.add(value)

    def add_annotations(
        self,
        annotations: Union[
            str,
            "Resource",
            Tuple[str, str],
            "StandardizedAnnotation",
            "CustomAnnotation",
            Iterable[
                Union[
                    str,
                    "Resource",
                    Tuple[str, str],
                    "StandardizedAnnotation",
                    "CustomAnnotation",
                ]
            ],
        ],
    ):
        """Add annotations to the metadata object by inferring the annotation type.

        If a StandardizedAnnotation or CustomAnnotation object is present, it is added
        to the `standardized` or `custom` attribute, respectively. Strings, tuples, and
        Resource objects are added using the `SimplifiedAnnotationInterface`, which
        means that they are added to the `standardized` attribute using a default
        qualifier (typically Qualifier.Biological_is).

        Parameters
        ----------
        annotations: (list of) StandardizedAnnotation, CustomAnnotation, str, tuple
            Annotations to add to the metadata.

        See Also
        --------
        StandardizedAnnotationStore.add
        CustomAnnotationStore.add
        SimplifiedAnnotationInterface.add
        """
        from cobra.core.metadata import (
            CustomAnnotation,
            Resource,
            StandardizedAnnotation,
        )

        if isinstance(
            annotations,
            (str, tuple, Resource, StandardizedAnnotation, CustomAnnotation),
        ):
            annotations = [annotations]

        for annotation in annotations:
            if isinstance(annotation, (str, tuple, Resource)):
                self.annotation.add(annotation)
            elif isinstance(annotation, StandardizedAnnotation):
                self.metadata.standardized.add([annotation])
            elif isinstance(annotation, CustomAnnotation):
                self.metadata.custom.add(annotation)
            else:
                raise TypeError(
                    "Could not convert object to annotation: "
                    f"{annotation} ({type(annotation)})"
                )

    def remove_annotations(
        self,
        annotations: Union[
            "StandardizedAnnotation",
            "CustomAnnotation",
            Iterable[Union["StandardizedAnnotation", "CustomAnnotation"]],
        ],
    ):
        """Remove an annotation from the Metadata object.

        This method only accepts StandardizedAnnotation and CustomAnnotation objects,
        or a list thereof, and removes them from the corresponding annotation stores at
        Metadata.standardized and Metadata.custom.

        Parameters
        ----------
        annotations: (list of) StandardizedAnnotation or CustomAnnotation objects
            Annotations to remove from metadata.

        See Also
        --------
        StandardizedAnnotationStore.remove
        CustomAnnotationStore.remove
        """
        from cobra.core.metadata import CustomAnnotation, StandardizedAnnotation

        if isinstance(annotations, StandardizedAnnotation) or isinstance(
            annotations, CustomAnnotation
        ):
            annotations = [annotations]

        for annotation in annotations:
            if isinstance(annotation, StandardizedAnnotation):
                self.metadata.standardized.remove(annotation)
            elif isinstance(annotation, CustomAnnotation):
                self.metadata.custom.remove(annotation)
            else:
                raise TypeError(
                    "All annotations should be of type StandardizedAnnotation "
                    f"or CustomAnnotation, not: {type(annotation)}."
                )

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
