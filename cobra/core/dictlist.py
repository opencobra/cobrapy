# -*- coding: utf-8 -*-

from __future__ import absolute_import

import re
from itertools import islice

from numpy import bool_
from six import PY3, iteritems, string_types


class DictList(list):
    """A combined dict and list

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
            raise TypeError("takes at most 1 argument (%d given)" % len(args))
        super(DictList, self).__init__(self)
        self._dict = {}
        if len(args) == 1:
            other = args[0]
            if isinstance(other, DictList):
                list.extend(self, other)
                self._dict = other._dict.copy()
            else:
                self.extend(other)

    def has_id(self, id):
        return id in self._dict

    def _check(self, id):
        """make sure duplicate id's are not added.
        This function is called before adding in elements.

        """
        if id in self._dict:
            raise ValueError("id %s is already present in list" % str(id))

    def _generate_index(self):
        """rebuild the _dict index"""
        self._dict = {v.id: k for k, v in enumerate(self)}

    def get_by_id(self, id):
        """return the element with a matching id"""
        return list.__getitem__(self, self._dict[id])

    def list_attr(self, attribute):
        """return a list of the given attribute for every object"""
        return [getattr(i, attribute) for i in self]

    def get_by_any(self, iterable):
        """
        Get a list of members using several different ways of indexing

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
        def get_item(item):
            if isinstance(item, int):
                return self[item]
            elif isinstance(item, string_types):
                return self.get_by_id(item)
            elif item in self:
                return item
            else:
                raise TypeError("item in iterable cannot be '%s'" % type(item))

        if not isinstance(iterable, list):
            iterable = [iterable]
        return [get_item(item) for item in iterable]

    def query(self, search_function, attribute=None):
        """Query the list

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
        >>> import cobra.test
        >>> model = cobra.test.create_test_model('textbook')
        >>> model.reactions.query(lambda x: x.boundary)
        >>> import re
        >>> regex = re.compile('^g', flags=re.IGNORECASE)
        >>> model.metabolites.query(regex, attribute='name')
        """
        def select_attribute(x):
            if attribute is None:
                return x
            else:
                return getattr(x, attribute)

        try:
            # if the search_function is a regular expression
            regex_searcher = re.compile(search_function)

            if attribute is not None:
                matches = (
                    i for i in self if
                    regex_searcher.findall(select_attribute(i)) != [])

            else:
                # Don't regex on objects
                matches = (
                    i for i in self if
                    regex_searcher.findall(getattr(i, 'id')) != [])

        except TypeError:
            matches = (
                i for i in self if search_function(select_attribute(i)))

        results = self.__class__()
        results._extend_nocheck(matches)
        return results

    def _replace_on_id(self, new_object):
        """Replace an object by another with the same id."""
        the_id = new_object.id
        the_index = self._dict[the_id]
        list.__setitem__(self, the_index, new_object)

    # overriding default list functions with new ones
    def append(self, object):
        """append object to end"""
        the_id = object.id
        self._check(the_id)
        self._dict[the_id] = len(self)
        list.append(self, object)

    def union(self, iterable):
        """adds elements with id's not already in the model"""
        _dict = self._dict
        append = self.append
        for i in iterable:
            if i.id not in _dict:
                append(i)

    def extend(self, iterable):
        """extend list by appending elements from the iterable"""
        # Sometimes during initialization from an older pickle, _dict
        # will not have initialized yet, because the initialization class was
        # left unspecified. This is an issue because unpickling calls
        # DictList.extend, which requires the presence of _dict. Therefore,
        # the issue is caught and addressed here.
        if not hasattr(self, "_dict") or self._dict is None:
            self._dict = {}
        _dict = self._dict
        current_length = len(self)
        list.extend(self, iterable)
        for i, obj in enumerate(islice(self, current_length, None),
                                current_length):
            the_id = obj.id
            if the_id not in _dict:
                _dict[the_id] = i
            else:
                # undo the extend and raise an error
                self = self[:current_length]
                self._check(the_id)
                # if the above succeeded, then the id must be present
                # twice in the list being added
                raise ValueError("id '%s' at index %d is non-unique. "
                                 "Is it present twice?" % (str(the_id), i))

    def _extend_nocheck(self, iterable):
        """extends without checking for uniqueness

        This function should only be used internally by DictList when it
        can guarantee elements are already unique (as in when coming from
        self or other DictList). It will be faster because it skips these
        checks.

        """
        current_length = len(self)
        list.extend(self, iterable)
        _dict = self._dict
        if current_length is 0:
            self._generate_index()
            return
        for i, obj in enumerate(islice(self, current_length, None),
                                current_length):
            _dict[obj.id] = i

    def __sub__(self, other):
        """x.__sub__(y) <==> x - y

        Parameters
        ----------
        other : iterable
            other must contain only unique id's present in the list
        """
        total = DictList()
        total.extend(self)
        for item in other:
            total.remove(item)
        return total

    def __isub__(self, other):
        """x.__sub__(y) <==> x -= y

        Parameters
        ----------
        other : iterable
            other must contain only unique id's present in the list
        """

        for item in other:
            self.remove(item)
        return self

    def __add__(self, other):
        """x.__add__(y) <==> x + y

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

    def __iadd__(self, other):
        """x.__iadd__(y) <==> x += y

        Parameters
        ----------
        other : iterable
            other must contain only unique id's whcih do not intersect
            with self

        """
        self.extend(other)
        return self

    def __reduce__(self):
        return (self.__class__, (), self.__getstate__(), self.__iter__())

    def __getstate__(self):
        """gets internal state

        This is only provided for backwards compatibility so older
        versions of cobrapy can load pickles generated with cobrapy. In
        reality, the "_dict" state is ignored when loading a pickle"""
        return {"_dict": self._dict}

    def __setstate__(self, state):
        """sets internal state

        Ignore the passed in state and recalculate it. This is only for
        compatibility with older pickles which did not correctly specify
        the initialization class"""
        self._generate_index()

    def index(self, id, *args):
        """Determine the position in the list

        id: A string or a :class:`~cobra.core.Object.Object`

        """
        # because values are unique, start and stop are not relevant
        if isinstance(id, string_types):
            try:
                return self._dict[id]
            except KeyError:
                raise ValueError("%s not found" % id)
        try:
            i = self._dict[id.id]
            if self[i] is not id:
                raise ValueError(
                    "Another object with the identical id (%s) found" % id.id)
            return i
        except KeyError:
            raise ValueError("%s not found" % str(id))

    def __contains__(self, object):
        """DictList.__contains__(object) <==> object in DictList

        object: str or :class:`~cobra.core.Object.Object`

        """
        if hasattr(object, "id"):
            the_id = object.id
        # allow to check with the object itself in addition to the id
        else:
            the_id = object
        return the_id in self._dict

    def __copy__(self):
        the_copy = DictList()
        list.extend(the_copy, self)
        the_copy._dict = self._dict.copy()
        return the_copy

    def insert(self, index, object):
        """insert object before index"""
        self._check(object.id)
        list.insert(self, index, object)
        # all subsequent entries now have been shifted up by 1
        _dict = self._dict
        for i, j in iteritems(_dict):
            if j >= index:
                _dict[i] = j + 1
        _dict[object.id] = index

    def pop(self, *args):
        """remove and return item at index (default last)."""
        value = list.pop(self, *args)
        index = self._dict.pop(value.id)
        # If the pop occured from a location other than the end of the list,
        # we will need to subtract 1 from every entry afterwards
        if len(args) == 0 or args == [-1]:  # removing from the end of the list
            return value
        _dict = self._dict
        for i, j in iteritems(_dict):
            if j > index:
                _dict[i] = j - 1
        return value

    def add(self, x):
        """Opposite of `remove`. Mirrors set.add"""
        self.extend([x])

    def remove(self, x):
        """.. warning :: Internal use only"""
        # Each item is unique in the list which allows this
        # It is much faster to do a dict lookup than n string comparisons
        self.pop(self.index(x))

    # these functions are slower because they rebuild the _dict every time
    def reverse(self):
        """reverse *IN PLACE*"""
        list.reverse(self)
        self._generate_index()

    def sort(self, cmp=None, key=None, reverse=False):
        """stable sort *IN PLACE*

        cmp(x, y) -> -1, 0, 1

        """
        if key is None:
            def key(i):
                return i.id
        if PY3:
            list.sort(self, key=key, reverse=reverse)
        else:
            list.sort(self, cmp=cmp, key=key, reverse=reverse)
        self._generate_index()

    def __getitem__(self, i):
        if isinstance(i, int):
            return list.__getitem__(self, i)
        elif isinstance(i, slice):
            selection = self.__class__()
            selection._extend_nocheck(list.__getitem__(self, i))
            return selection
        elif hasattr(i, "__len__"):
            if len(i) == len(self) and isinstance(i[0], (bool, bool_)):
                selection = self.__class__()
                result = (o for j, o in enumerate(self) if i[j])
                selection._extend_nocheck(result)
                return selection
            else:
                return self.__class__(list.__getitem__(self, i))
        else:
            return list.__getitem__(self, i)

    def __setitem__(self, i, y):
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
        # in case a rename has occured
        if self._dict.get(self[i].id) == i:
            self._dict.pop(self[i].id)
        the_id = y.id
        self._check(the_id)
        list.__setitem__(self, i, y)
        self._dict[the_id] = i

    def __delitem__(self, index):
        removed = self[index]
        list.__delitem__(self, index)
        if isinstance(removed, list):
            self._generate_index()
            return
        _dict = self._dict
        _dict.pop(removed.id)
        for i, j in iteritems(_dict):
            if j > index:
                _dict[i] = j - 1

    def __getslice__(self, i, j):
        return self.__getitem__(slice(i, j))

    def __setslice__(self, i, j, y):
        self.__setitem__(slice(i, j), y)

    def __delslice__(self, i, j):
        self.__delitem__(slice(i, j))

    def __getattr__(self, attr):
        try:
            return DictList.get_by_id(self, attr)
        except KeyError:
            raise AttributeError("DictList has no attribute or entry %s" %
                                 attr)

    def __dir__(self):
        # override this to allow tab complete of items by their id
        attributes = dir(self.__class__)
        attributes.append("_dict")
        attributes.extend(self._dict.keys())
        return attributes
