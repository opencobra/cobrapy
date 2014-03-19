from copy import copy, deepcopy
import re
from ..external.six import string_types

class DictList(list):
    """A combined dict and list that feels like a list, but has
    the speed benefits of a dict.  This may be eventually
    replaced by collections.OrderedDict.

    This was written to address the performance issues associated
    with searching, accessing, or iterating over a list in python
    that resulted in notable performance decays with COBRA for
    python.

    """
    def __init__(self, *args, **kwargs):
        list.__init__(self, *args, **kwargs)
        self._dict = {}
        self._object_dict = {}
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
        self._object_dict = {v.id: v for v in self}


    def get_by_id(self, id):
        """return the element with a matching id

        """
        return self._object_dict[id]

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
        the_id = y.id
        self._check(the_id)
        super(DictList, self).__setitem__(i, y)
        self._dict[the_id] = i
        self._object_dict[the_id] = y

    def _replace_on_id(self, new_object):
        """Allows one to replace an object by one with
        the same id.
        
        """
        the_id = new_object.id
        the_index = self._dict[the_id]
        super(DictList, self).__setitem__(the_index, new_object)
        self._object_dict[the_id] = new_object

    def append(self, object):
        the_id = object.id
        self._check(the_id)
        self._dict[the_id] = len(self)
        super(DictList, self).append(object)
        self._object_dict[the_id] = object

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
        self._object_dict.clear()
        the_copy = copy(super(DictList, self))
        self._generate_index()
        the_copy._generate_index()
        return the_copy


    def __deepcopy__(self, *args, **kwargs):
        _dict = self._dict
        _object_dict = self._object_dict
        self._dict = {}
        self._object_dict = {}
        new = deepcopy(super(DictList, self), *args, **kwargs)
        self._dict = _dict
        self._object_dict = _object_dict
        return new

    def insert(self, index, object):
        self._check(object.id)
        super(DictList, self).insert(index, object)
        self._object_dict[object.id] = object
        # all subsequent entries now have been shifted up by 1
        _dict = self._dict
        for i in _dict:
            j = _dict[i]
            if j >= index:
                _dict[i] = j + 1
        _dict[object.id] = index

    def pop(self, *args):
        value = super(DictList, self).pop(*args)
        index = self._dict.pop(value.id)
        self._object_dict.pop(value.id)
        # If the pop occured from a location other than the end of the list,
        # we will need to subtract 1 from every entry afterwards
        if len(args) == 0 or args == [-1]: # removing from the end of the list
            return
        _dict = self._dict
        for i in _dict:
            j = _dict[i]
            if j > index:
                _dict[i] = j - 1
        return value

    def remove(self, x):
        # Each item is unique in the list which allows this
        # It is much faster to do a dict lookup than n string comparisons
        self.pop(self.index(x))

    # these functions are slower because they rebuild the _dict every time
    def reverse(self, *args, **kwargs):
        super(DictList, self).reverse(*args, **kwargs)
        self._generate_index()

    def sort(self, *args, **kwargs):
        super(DictList, self).sort(*args, **kwargs)
        self._generate_index()

    def __setslice__(self, *args, **kwargs):
        super(DictList, self).__setslice__(*args, **kwargs)
        self._generate_index()

    def __delslice__(self, *args, **kwargs):
        super(DictList, self).__delslice__(*args, **kwargs)
        self._generate_index()

    def __delitem__(self, *args, **kwargs):
        super(DictList, self).__delitem__(*args, **kwargs)
        self._generate_index()

    #def __getattr__(self, attr):
    #    # makes items attributes as well
    #    # this runs after __getattribute__ has already run
    #    try:
    #        return self.get_by_id(attr)
    #    except KeyError:
    #        raise AttributeError("DictList has no attribute or key " + (attr))

    #def __dir__(self):
    #    # override this to allow tab complete of items by their id
    #    attributes = dir(self.__class__)
    #    attributes.extend(self.__dict__.keys())
    #    attributes.extend(self._dict.keys())
    #    return attributes
