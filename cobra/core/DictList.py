from copy import copy, deepcopy
import re
from ..external.six import string_types, iteritems

class DictList(list):
    """A combined dict and list that feels like a list, but has
    the speed benefits of a dict.  This may be eventually
    replaced by collections.OrderedDict.

    This was written to address the performance issues associated
    with searching, accessing, or iterating over a list in python
    that resulted in notable performance decays with COBRA for
    python.

    """
    _dict = {}
    def __init__(self, *args, **kwargs):
        list.__init__(self, *args, **kwargs)
        self._generate_index()

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
        return self[self._dict[id]]

    def list_attr(self, attribute):
        """return a list of the given attribute for every object

        """
        return [getattr(i, attribute) for i in self]

    def query(self, search_function, attribute="id"):
        """query the list

        search_function: this will be used to select which objects to return
        This can be:
            - a string, in which case any object.attribute containing
              the string will be returned
            - a compiled regular expression
            - a boolean function which takes one argument and returns True
              for desired values

        attribute: the attribute to be searched for (default is 'id').
                   If this is None, the object itself is used.

        returns: a list of objects which match the query
        """
        if attribute == None:
            select_attribute = lambda x : x
        else:
            select_attribute = lambda x: getattr(x, attribute)

        # if the search_function is a regular expression
        if isinstance(search_function, str):
            search_function = re.compile(search_function)
        if hasattr(search_function, "findall"):
            matches = [i for i in self
                if search_function.findall(select_attribute(i)) != []]
        else:
            matches = [i for i in self
                if search_function(select_attribute(i))]
        return DictList(matches)


    # overriding default list functions with new ones
    def __setitem__(self, i, y):
        if isinstance(i, slice):
            for obj in y:  # need to be setting to a list
                self._check(obj.id)
            list.__setitem__(self, i, y)
            self._generate_index()
            return
        self._dict.pop(self[i].id)
        the_id = y.id
        self._check(the_id)
        list.__setitem__(self, i, y)
        self._dict[the_id] = i

    def _replace_on_id(self, new_object):
        """Allows one to replace an object by one with
        the same id.
        
        """
        the_id = new_object.id
        the_index = self._dict[the_id]
        list.__setitem__(self, the_index, new_object)

    def append(self, object):
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
        append = self.append
        for i in iterable:
            append(i)

    def __add__(self, other, should_deepcopy=False):
        """
        other: an DictList

        should_deepcopy: Boolean. 

        """
        if should_deepcopy:
            from copy import deepcopy
            total = deepcopy(self)
        else:
            total = DictList()
            total.extend(self)
        total.extend(other)
        return total

    def __iadd__(self, other):
        self.extend(other)
        return self

    def index(self, id):
        """
        id: A string or a :class:`~cobra.core.Object`
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
            raise ValueError("%s not found" % str(i))

    def __contains__(self, object):
        """DictList.__contains__(object) <==> object in DictList
        object can either be the object to search for itself, or 
        simply the id

        """
        if hasattr(object, "id"):
            the_id = object.id
        # allow to check with the object itself in addition to the id
        else:
            the_id = object
        return the_id in self._dict

    def __copy__(self):
        self._dict.clear()
        the_copy = copy(super(DictList, self))
        self._generate_index()
        the_copy._generate_index()
        return the_copy

    def __deepcopy__(self, *args, **kwargs):
        _dict = self._dict
        self._dict = {}
        new = deepcopy(super(DictList, self), *args, **kwargs)
        self._dict = _dict
        return new

    def insert(self, index, object):
        self._check(object.id)
        list.insert(self, index, object)
        # all subsequent entries now have been shifted up by 1
        _dict = self._dict
        for i in _dict:
            j = _dict[i]
            if j >= index:
                _dict[i] = j + 1
        _dict[object.id] = index

    def pop(self, *args):
        value = list.pop(self, *args)
        index = self._dict.pop(value.id)
        # If the pop occured from a location other than the end of the list,
        # we will need to subtract 1 from every entry afterwards
        if len(args) == 0 or args == [-1]: # removing from the end of the list
            return
        _dict = self._dict
        for i, j in iteritems(_dict):
            if j > index:
                _dict[i] = j - 1
        return value

    def remove(self, x):
        # Each item is unique in the list which allows this
        # It is much faster to do a dict lookup than n string comparisons
        self.pop(self.index(x))

    # these functions are slower because they rebuild the _dict every time
    def reverse(self, *args, **kwargs):
        list.reverse(self, *args, **kwargs)
        self._generate_index()

    def sort(self, key=None, **kwargs):
        if key is None:
            key = lambda i: i.id
        list.sort(self, key=key, **kwargs)
        self._generate_index()

    def __setslice__(self, *args, **kwargs):
        list.__setslice__(self, *args, **kwargs)
        self._generate_index()

    def __delslice__(self, *args, **kwargs):
        list.__delslice__(self, *args, **kwargs)
        self._generate_index()

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

    def __getattr__(self, attr):
        try:
            return self[self.__dict__["_dict"][attr]]
        except KeyError:
            raise AttributeError("DictList has no attribute or entry %s" % \
                (attr))

    def __dir__(self):
        # override this to allow tab complete of items by their id
        attributes = dir(self.__class__)
        attributes.extend(self.__dict__.keys())
        attributes.extend(self._dict.keys())
        return attributes
