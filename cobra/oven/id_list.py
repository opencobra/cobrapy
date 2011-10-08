import copy
import re


def id(object):
    """return an id for the object"""
    return object.id


class id_list(list):
    """an extension to the list class which allows faster searching
    This list is good for storing a list of objects with a unique
    id attribute
    """
    def __init__(self, *args, **kwargs):
        if "type" in kwargs:
            _object_type = kwargs[type]
            kwargs.pop(type)
        else:
            object_type = None
        list.__init__(self, *args, **kwargs)
        self._dict = {}
        self._generate_index()

    def _check(self, id):
        """make sure duplicate id's are not added.
        This function is called before adding in elements
        """
        if id in self._dict:
            raise ValueError("id %s is already present in list" % id)

    def _generate_index(self):
        """rebuild the _dict index"""
        self._dict.clear()
        for i, obj in enumerate(self):
            self._dict[id(obj)] = i

    def get_by_id(self, id):
        """return the element with a matching id"""
        return self[self._dict[id]]

    def list_attr(self, attribute):
        """return a list of the given attribute for every object"""
        return [getattr(object, attribute) for object in self]

    def query(self, search_function, attribute="id"):
        """query the list

        search_function: this will be used to select which objects to return
        This can be:
            - a string, in which case any object.attribute containing
              the string will be returned
            - a compiled regular expression
            - a boolean function which takes one argument and returns True
              for desired values

        attribute: the attribute to be searched for (default is "id").
                   If this is None, the object itself is used.

        returns: a list of objects which match the query
        """
        if attribute is None:
            attr_select = lambda x: x
        else:
            attr_select = lambda x: getattr(obj, attribute)

        # if the search_function is a regular expression
        match_list = id_list()
        if type(search_function) is str:
            search_function = re.compile(search_function)
        if hasattr(search_function, "findall"):
            for obj in self:
                if search_function.findall(attr_select(obj)) != []:
                    match_list.append(obj)
        else:
            for obj in self:
                if search_function(attr_select(obj)):
                    match_list.append(obj)
        return match_list

    # overriding default list functions with new ones
    def __setitem__(self, i, y):
        self._check(id(y))
        super(id_list, self).__setitem__(i, y)
        self._dict[id(y)] = i

    def append(self, object):
        self._check(id(object))
        super(id_list, self).append(object)
        self._dict[id(object)] = len(self)

    def union(self, iterable):
        """adds elements with id's not already in the model"""
        for i in iterable:
            if not i.id in self._dict:
                self.append(i)

    def extend(self, iterable):
        for i in iterable:
            self.append(i)

    def __add__(self, other):
        sum = copy.deepcopy(other)  # should this be deepcopy or shallow?
        sum.extend(other)
        sum._generate_index()
        return sum

    def __iadd__(self, other):
        self.extend(other)
        self._generate_index()
        return self

    def index(self, id):
        # because values are unique, start and stop are not relevant
        return self._dict[id]

    def __contains__(self, object):
        """id_list.__contains__(object) <==> object in id_list
        object can either be the object to search for itself, or
         simply the id"""
        if hasattr(object, "id"):
            id = id(object)
        # allow to check with the object itself in addition to the id
        else:
            id = object
        return id in self._dict

    def __copy__(self):
        self._dict.clear()
        c = copy.copy(super(id_list, self))
        self._generate_index()
        c._generate_index()
        return c

    def __deepcopy__(self, *args, **kwargs):
        self._dict.clear()
        c = copy.deepcopy(super(id_list, self), *args, **kwargs)
        self._generate_index()
        c._generate_index()
        return c

    # these functions are slower because they rebuild the _dict every time
    # TODO: speed up
    def insert(index, object):
        self._check(id(object))
        super(id_list, self).insert(index, object)
        self._generate_index()

    def pop(self, *args, **kwargs):
        value = super(id_list, self).pop(*args, **kwargs)
        self._generate_index()
        return value

    def remove(self, *args, **kwargs):
        super(id_list, self).remove(*args, **kwargs)
        self._generate_index()

    def reverse(self, *args, **kwargs):
        super(id_list, self).reverse(*args, **kwargs)
        self._generate_index()

    def sort(self, *args, **kwargs):
        super(id_list, self).sort(*args, **kwargs)
        self._generate_index()

    def __setslice__(self, *args, **kwargs):
        super(id_list, self).__setslice__(*args, **kwargs)
        self._generate_index()

    def __delslice__(self, *args, **kwargs):
        super(id_list, self).__delslice__(*args, **kwargs)
        self._generate_index()

    def __delitem__(self, *args, **kwargs):
        super(id_list, self).__delitem__(*args, **kwargs)
        self._generate_index()

    # def __getattribute__(self, *args, **kwargs):
        """This will allow the id_list class to list attributes of the
        containing objects. For example, if the_list is an instance
        of id_list, calling the_list.id will return a list of all the id's"""
        # name = args[0]
        # if hasattr(super(id_list, self), name):
            # return getattr(super(id_list, self), name)
        # if name[-1] == "_" or not hasattr(, name):
            # raise AttributeError(name)
        # if self.object_type is not None:
            # pass  # TODO - finish
