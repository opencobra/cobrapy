from copy import copy, deepcopy
import re

def get_id(object):
    """return an id for the object

    This allows the function to be generalize to non-cobra.core objects,
    however, this added function call slows things down.

    """
    return object.id

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
        
    
    def _check(self, the_id):
        """make sure duplicate id's are not added.
        This function is called before adding in elements.

        """
        if the_id in self._dict:
            raise ValueError, "the_id is already present in list"
    
    def _generate_index(self):
        """rebuild the _dict index

        """
        self._dict = {}
        self._object_dict = {}
        [(self._dict.update({v.id: k}),
          self._object_dict.update({v.id: v}))
         for k, v in enumerate(self)]
            
    def get_by_id(self, the_id):
        """return the element with a matching id

        """
        return self._object_dict[the_id]
        
    def list_attr(self, attribute):
        """return a list of the given attribute for every object

        """
        return [getattr(i, attribute)
                for i in self]
        
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
            select_attribute = lambda x: getattr(the_object, attribute)
        
        # if the search_function is a regular expression
        match_list = DictList()
        if isinstance(search_function, str):
            search_function = re.compile(search_function)
        if hasattr(search_function, "findall"):
            for the_object in self:
                if search_function.findall(select_attribute(the_object)) != []:
                    match_list.append(the_object)
        else:
            for the_object in self:
                if search_function(select_attribute(the_object)):
                    match_list.append(the_object)
        return match_list
    
    # overriding default list functions with new ones
    def __setitem__(self, i, y):
        the_id = get_id(y)
        self._check(the_id)
        super(DictList, self).__setitem__(i, y)
        self._dict[the_id] = i
        self._object_dict[the_id] = y
        
    def append(self, the_object):
        the_id = get_id(the_object)
        self._check(the_id)
        super(DictList, self).append(the_object)
        self._dict[the_id] = len(self)
        self._object_dict[the_id] = the_object
    
    def union(self, iterable):
        """adds elements with id's not already in the model"""
        [self.append(i)
         for i in iterable
         if get_id(i) not in self._dict]

    
    def extend(self, iterable):
        [self.append(i) for i in iterable]

            
    def __add__(self, other, deepcopy=True):
        """
        other: an DictList

        deepcopy: Boolean.  Allow for shallow copying, however,
        this can cause problems if one doesn't know that the
        items are referenceable from different id

        
        """
        if deepcopy:
            sum = deepcopy(other) # should this be deepcopy or shallow?
        else:
            sum = other
        sum.extend(other)
        sum._generate_index()
        return sum
    
    def __iadd__(self, other):
        self.extend(other)
        self._generate_index()
        return self
    
    
    def index(self, the_id):
        # because values are unique, start and stop are not relevant
        return self._dict[the_id]
        
    def __contains__(self, the_object):
        """DictList.__contains__(object) <==> object in DictList
        object can either be the object to search for itself, or 
        simply the id

        """
        if hasattr(object, "id"):
            the_id = get_id(object)
        # allow to check with the object itself in addition to the id
        else:
            the_id = the_object
        return self._dict.has_key(the_id)
        
    def __copy__(self):
        self._dict.clear()
        self._object_dict.clear()
        the_copy = copy(super(DictList, self))
        self._generate_index()
        the_copy._generate_index()
        return the_copy
    
    def __deepcopy__(self, *args, **kwargs):
        self._dict.clear()
        self._object_dict.clear()
        the_copy = deepcopy(super(DictList, self), *args, **kwargs)
        self._generate_index()
        the_copy._generate_index()
        return the_copy
    
    # these functions are slower because they rebuild the _dict every time
    # TODO: speed up
    def insert(index, object):
        self._check(get_id(object))
        super(DictList, self).insert(index, object)
        self._generate_index()
    
    def pop(self, *args, **kwargs):
        value = super(DictList, self).pop(*args, **kwargs)
        self._generate_index()
        return value
    
    def remove(self, *args, **kwargs):
        super(DictList, self).remove(*args, **kwargs)
        self._generate_index()
   
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


if __name__ == '__main__':
    from cPickle import load
    from time import time
    from numpy import round
    from copy import deepcopy
    from cobra import Model
    from collections import defaultdict
    import cobra.core.Model_2
    reload(cobra.core.Model_2)
    from cobra.core.Model_2 import Model as Model_2
    model_dict = {'list': Model,
                  'DictList': Model_2}
    
    from cobra.manipulation import initialize_growth_medium
    test_directory = '../../test/data/'
    with open(test_directory + 'salmonella.pickle') as in_file:
        cobra_model = load(in_file)
    #
    original_model = deepcopy(cobra_model)
    print 'testing DictList'
    start_time = time()
    the_list = DictList(cobra_model.metabolites)
    print 'made a DictList in %f seconds from cobra_model.metabolites'%(time()-start_time)
    time_results = defaultdict(dict)
    the_models =  {}
    for the_kind, Model in model_dict.items():
        print '\nTesting %s'%the_kind
        start_time = time()
        the_copy = Model('asdfd')
        the_models[the_kind] = the_copy
        the_test = 'add_reactions'
        the_copy.add_reactions(cobra_model.reactions)
        the_time = time_results[the_test][the_kind] = time()-start_time
        print '%s %s time: %f seconds'%(the_kind, the_test,
                                        the_time)
        start_time = time()
        the_test = 'copy'
        asdf = the_copy.copy(print_time=True)
        the_time = time_results[the_test][the_kind] = time()-start_time
        print '%s %s time: %f seconds'%(the_kind, the_test,
                                                  the_time)


    print '\n%% improvement with DictList (- is bad)'
    for the_test, the_results in time_results.items():
        print '\t%s: %1.1f'%(the_test,
                             100*(1-the_results['DictList']/
                                  the_results['list']))

