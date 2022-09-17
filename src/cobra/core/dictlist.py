"""Define the DictList class."""

import re
from itertools import islice
from typing import (
    Any,
    Callable,
    Iterable,
    Iterator,
    List,
    Optional,
    Pattern,
    Tuple,
    Type,
    Union,
)

import numpy as np

from .object import Object


class DictList(list):
    """
    Define a combined dict and list.

    This object behaves like a list, but has the O(1) speed
    benefits of a dict when looking up elements by their id.

    """

    def __init__(self, *args):
        """Instantiate a combined dict and list.

        Parameters
        ----------
        args : iterable
            iterable as single argument to create new DictList from

        """
        if len(args) > 2:
            raise TypeError(f"takes at most 1 argument ({len(args):d} given)")
        super(DictList, self).__init__(self)
        self._dict = {}
        if len(args) == 1:
            other = args[0]
            if isinstance(other, DictList):
                list.extend(self, other)
                self._dict = other._dict.copy()
            else:
                self.extend(other)

    # noinspection PyShadowingBuiltins
    def has_id(self, id: Union[Object, str]) -> bool:
        """Check if id is in DictList."""
        return id in self._dict

    # noinspection PyShadowingBuiltins
    def _check(self, id: Union[Object, str]) -> None:
        """Make sure duplicate id's are not added.

        This function is called before adding in elements.

        """
        if id in self._dict:
            raise ValueError(f"id {str(id)} is already present in list")

    def _generate_index(self) -> None:
        """Rebuild the _dict index."""
        self._dict = {v.id: k for k, v in enumerate(self)}

    # noinspection PyShadowingBuiltins
    def get_by_id(self, id: Union[Object, str]) -> Object:
        """Return the element with a matching id."""
        return list.__getitem__(self, self._dict[id])

    def list_attr(self, attribute: str) -> list:
        """Return a list of the given attribute for every object."""
        return [getattr(i, attribute) for i in self]

    def get_by_any(self, iterable: List[Union[str, Object, int]]) -> list:
        """Get a list of members using several different ways of indexing.

        Parameters
        ----------
        iterable : list (if not, turned into single element list)
            list where each element is either int (referring to an index in
            in this DictList), string (a id of a member in this DictList) or
            member of this DictList for pass-through

        Returns
        -------
        list
            a list of members
        """

        def get_item(item: Any) -> Any:
            if isinstance(item, int):
                return self[item]
            elif isinstance(item, str):
                return self.get_by_id(item)
            elif item in self:
                return item
            else:
                raise TypeError(f"item in iterable cannot be '{type(item)}'")

        if not isinstance(iterable, list):
            iterable = [iterable]
        return [get_item(item) for item in iterable]

    def query(
        self,
        search_function: Union[str, Pattern, Callable],
        attribute: Union[str, None] = None,
    ) -> "DictList":
        """Query the list.

        Parameters
        ----------
        search_function : a string, regular expression or function
            Used to find the matching elements in the list.
            - a regular expression (possibly compiled), in which case the
            given attribute of the object should match the regular expression.
            - a function which takes one argument and returns True for
            desired values

        attribute : string or None
            the name attribute of the object to passed as argument to the
            `search_function`. If this is None, the object itself is used.

        Returns
        -------
        DictList
            a new list of objects which match the query

        Examples
        --------
        >>> from cobra.io import load_model
        >>> model = load_model('textbook')
        >>> model.reactions.query(lambda x: x.boundary)
        >>> import re
        >>> regex = re.compile('^g', flags=re.IGNORECASE)
        >>> model.metabolites.query(regex, attribute='name')
        """

        def select_attribute(x: Optional[Any]) -> Any:
            if attribute is None:
                return x
            else:
                return getattr(x, attribute)

        try:
            # if the search_function is a regular expression
            regex_searcher = re.compile(search_function)

            if attribute is not None:
                matches = (
                    i for i in self if regex_searcher.findall(select_attribute(i)) != []
                )

            else:
                # Don't regex on objects
                matches = (i for i in self if regex_searcher.findall(i.id) != [])

        except TypeError:
            matches = (i for i in self if search_function(select_attribute(i)))

        results = self.__class__()
        results._extend_nocheck(matches)
        return results

    def _replace_on_id(self, new_object: Object) -> None:
        """Replace an object by another with the same id."""
        the_id = new_object.id
        the_index = self._dict[the_id]
        list.__setitem__(self, the_index, new_object)

    # overriding default list functions with new ones
    def append(self, entity: Object) -> None:
        """Append object to end."""
        the_id = entity.id
        self._check(the_id)
        self._dict[the_id] = len(self)
        list.append(self, entity)

    def union(self, iterable: Iterable[Object]) -> None:
        """Add elements with id's not already in the model."""
        _dict = self._dict
        append = self.append
        for i in iterable:
            if i.id not in _dict:
                append(i)

    def extend(self, iterable: Iterable[Object]) -> None:
        """Extend list by appending elements from the iterable.

        Sometimes during initialization from an older pickle, _dict
        will not have initialized yet, because the initialization class was
        left unspecified. This is an issue because unpickling calls
        DictList.extend, which requires the presence of _dict. Therefore,
        the issue is caught and addressed here.

        Parameters
        ----------
        iterable : Iterable
        """
        if getattr(self, "_dict", None) is None:
            self._dict = {}
        _dict = self._dict
        current_length = len(self)
        list.extend(self, iterable)
        for i, obj in enumerate(islice(self, current_length, None), current_length):
            the_id = obj.id
            if the_id not in _dict:
                _dict[the_id] = i
            else:
                # undo the extend and raise an error
                self = self[:current_length]
                self._check(the_id)
                # if the above succeeded, then the id must be present
                # twice in the list being added
                raise ValueError(
                    f"id '{str(the_id)}' at index {i :d} is non-unique. "
                    f"Is it present twice?"
                )

    def _extend_nocheck(self, iterable: Iterable[Object]) -> None:
        """Extend without checking for uniqueness.

        This function should only be used internally by DictList when it
        can guarantee elements are already unique (as in when coming from
        self or other DictList). It will be faster because it skips these
        checks.

        Parameters
        ----------
        iterable : Iterable

        """
        current_length = len(self)
        list.extend(self, iterable)
        _dict = self._dict
        if not current_length:
            self._generate_index()
            return
        for i, obj in enumerate(islice(self, current_length, None), current_length):
            _dict[obj.id] = i

    def __sub__(self, other: Iterable[Object]) -> "DictList":
        """Remove a value or values, and returns the new DictList.

        x.__sub__(y) <==> x - y

        Parameters
        ----------
        other : iterable
            other must contain only unique id's present in the list
        Returns
        -------
        total: DictList
            new DictList with item(s) removed
        """
        total = DictList()
        total.extend(self)
        for item in other:
            total.remove(item)
        return total

    def __isub__(self, other: Iterable[Object]) -> "DictList":
        """Remove a value or values in place.

        x.__sub__(y) <==> x -= y

        Parameters
        ----------
        other : iterable
            other must contain only unique id's present in the list
        """
        for item in other:
            self.remove(item)
        return self

    def __add__(self, other: Iterable[Object]) -> "DictList":
        """Add item while returning a new DictList.

        x.__add__(y) <==> x + y

        Parameters
        ----------
        other : iterable
            other must contain only unique id's which do not intersect
            with self
        """
        total = DictList()
        total.extend(self)
        total.extend(other)
        return total

    def __iadd__(self, other: Iterable[Object]) -> "DictList":
        """Add item while returning the same DictList.

        x.__iadd__(y) <==> x += y

        Parameters
        ----------
        other : iterable
            other must contain only unique id's whcih do not intersect
            with self

        """
        self.extend(other)
        return self

    def __reduce__(self) -> Tuple[Type["DictList"], Tuple, dict, Iterator]:
        """Return a reduced version of DictList.

        This reduced version details the class, an empty Tuple, a dictionary of the
        state and an iterator to go over the DictList.
        """
        return self.__class__, (), self.__getstate__(), self.__iter__()

    def __getstate__(self) -> dict:
        """Get internal state.

        This is only provided for backwards compatibility so older
        versions of cobrapy can load pickles generated with cobrapy. In
        reality, the "_dict" state is ignored when loading a pickle
        """
        return {"_dict": self._dict}

    def __setstate__(self, state: dict) -> None:
        """Pretend to set internal state. Actually recalculates.

        Ignore the passed in state and recalculate it. This is only for
        compatibility with older pickles which did not correctly specify
        the initialization class
        """
        self._generate_index()

    # noinspection PyShadowingBuiltins
    def index(self, id: Union[str, Object], *args) -> int:
        """Determine the position in the list.

        Parameters
        ----------
        id: A string or a :class:`~cobra.core.Object.Object`

        """
        # because values are unique, start and stop are not relevant
        if isinstance(id, str):
            try:
                return self._dict[id]
            except KeyError:
                raise ValueError(f"{id} not found")
        try:
            i = self._dict[id.id]
            if self[i] is not id:
                raise ValueError(
                    f"Another object with the identical id ({id.id}) found"
                )
            return i
        except KeyError:
            raise ValueError(f"{str(id)} not found")

    def __contains__(self, entity: Union[str, Object]) -> bool:
        """Ask if the DictList contain an entity.

        DictList.__contains__(entity) <==> entity in DictList

        Parameters
        ----------
        entity: str or :class:`~cobra.core.Object.Object`

        """
        if hasattr(entity, "id"):
            the_id = entity.id
        # allow to check with the object itself in addition to the id
        else:
            the_id = entity
        return the_id in self._dict

    def __copy__(self) -> "DictList":
        """Copy the DictList into a new one."""
        the_copy = DictList()
        list.extend(the_copy, self)
        the_copy._dict = self._dict.copy()
        return the_copy

    def insert(self, index: int, entity: Object) -> None:
        """Insert entity before index."""
        self._check(entity.id)
        list.insert(self, index, entity)
        # all subsequent entries now have been shifted up by 1
        _dict = self._dict
        for i, j in _dict.items():
            if j >= index:
                _dict[i] = j + 1
        _dict[entity.id] = index

    def pop(self, *args) -> Object:
        """Remove and return item at index (default last)."""
        value = list.pop(self, *args)
        index = self._dict.pop(value.id)
        # If the pop occurred from a location other than the end of the list,
        # we will need to subtract 1 from every entry afterwards
        if len(args) == 0 or args == [-1]:  # removing from the end of the list
            return value
        _dict = self._dict
        for i, j in _dict.items():
            if j > index:
                _dict[i] = j - 1
        return value

    def add(self, x: Object) -> None:
        """Opposite of `remove`. Mirrors set.add."""
        self.extend([x])

    def remove(self, x: Union[str, Object]) -> None:
        """.. warning :: Internal use only.

        Each item is unique in the list which allows this
        It is much faster to do a dict lookup than n string comparisons
        """
        self.pop(self.index(x))

    # these functions are slower because they rebuild the _dict every time
    def reverse(self) -> None:
        """Reverse *IN PLACE*."""
        list.reverse(self)
        self._generate_index()

    def sort(
        self, cmp: Callable = None, key: Callable = None, reverse: bool = False
    ) -> None:
        """Stable sort *IN PLACE*.

        cmp(x, y) -> -1, 0, 1

        """
        if key is None:

            def key(i):
                return i.id

        list.sort(self, key=key, reverse=reverse)

        self._generate_index()

    def __getitem__(
        self, i: Union[int, slice, Iterable, Object, "DictList"]
    ) -> Union["DictList", Object]:
        """Get item from DictList."""
        if isinstance(i, int):
            return list.__getitem__(self, i)
        elif isinstance(i, slice):
            selection = self.__class__()
            selection._extend_nocheck(list.__getitem__(self, i))
            return selection
        elif hasattr(i, "__len__"):
            if len(i) == len(self) and isinstance(i[0], (bool, np.bool)):
                selection = self.__class__()
                result = (o for j, o in enumerate(self) if i[j])
                selection._extend_nocheck(result)
                return selection
            else:
                return self.__class__(list.__getitem__(self, i))
        else:
            return list.__getitem__(self, i)

    def __setitem__(self, i: Union[slice, int], y: Union[list, Object]) -> None:
        """Set an item via index or slice.

        Parameters
        ----------
        i : slice, int
            i can be slice or int. If i is a slice, y needs to be a list
        y: list, Object
            Object to set as
        """
        if isinstance(i, slice):
            # In this case, y needs to be a list. We will ensure all
            # the id's are unique
            for obj in y:  # need to be setting to a list
                self._check(obj.id)
                # Insert a temporary placeholder so we catch the presence
                # of a duplicate in the items being added
                self._dict[obj.id] = None
            list.__setitem__(self, i, y)
            self._generate_index()
            return
        # in case a rename has occurred
        if self._dict.get(self[i].id) == i:
            self._dict.pop(self[i].id)
        the_id = y.id
        self._check(the_id)
        list.__setitem__(self, i, y)
        self._dict[the_id] = i

    def __delitem__(self, index: int) -> None:
        """Remove item from DictList."""
        removed = self[index]
        list.__delitem__(self, index)
        if isinstance(removed, list):
            self._generate_index()
            return
        _dict = self._dict
        _dict.pop(removed.id)
        for i, j in _dict.items():
            if j > index:
                _dict[i] = j - 1

    def __getslice__(self, i: int, j: int) -> "DictList":
        """Get a slice from it to j of DictList."""
        return self.__getitem__(slice(i, j))

    def __setslice__(self, i: int, j: int, y: Union[list, Object]) -> None:
        """Set slice, where y is an iterable."""
        self.__setitem__(slice(i, j), y)

    def __delslice__(self, i: int, j: int) -> None:
        """Remove slice."""
        self.__delitem__(slice(i, j))

    def __getattr__(self, attr: Any) -> Any:
        """Get an attribute by id."""
        try:
            return DictList.get_by_id(self, attr)
        except KeyError:
            raise AttributeError(f"DictList has no attribute or entry {attr}")

    def __dir__(self) -> list:
        """Directory of the DictList.

        Override this to allow tab complete of items by their id.

        Returns
        -------
        attributes: list
            A list of attributes/entities.
        """
        attributes = dir(self.__class__)
        attributes.append("_dict")
        attributes.extend(self._dict.keys())
        return attributes
