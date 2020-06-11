# -*- coding: utf-8 -*-

"""
Define the Controlled Vocabulary term class for refering to external
resources
"""

from __future__ import absolute_import


class CVTerm(dict):
    """
    Class representation of Controlled Vocabulary term inside Annotation.
    It will look similar to a dictionary, but will have restrictions on the
    type of keys and values

    Parameters
    ----------
    cvterm : dict
        A dictionary that maps the provider to its corresponding CVList

    Attributes
    ----------
        The providers are set as keys for accessing their CVLists

    """

    def __init__(self, cvterm={}):
        if not isinstance(cvterm, dict):
            raise TypeError("The annotation data must be in a dict form")
        else:
            for key, value in cvterm.items():
                if not isinstance(key, str):
                    raise TypeError("the provider must be of type string")
                if isinstance(value, list):
                    dict.__setitem__(self, key, self.CVList(value))
                else:
                    dict.__setitem__(self, key, value)

    def __getitem__(self, key):
        return dict.__getitem__(self, key)

    def __setitem__(self, key, value):
        """Make sure that key passed is of type string and value
           passed confirms to CVList type (CVList or list)
        """
        if not isinstance(key, str):
            raise TypeError("The key passed must be a string")
        if isinstance(value, list):
            dict.__setitem__(self, key, self.CVList(value))
        elif isinstance(value, self.CVList):
            dict.__setitem__(self, key, value)
        else:
            raise TypeError("The value passed does not confirm to CVList type")

    def __delitem__(self, key):
        dict.__delitem__(self, key)

    def __iter__(self):
        return dict.__iter__(self)

    def __len__(self):
        return dict.__len__(self)

    def __contains__(self, x):
        return dict.__contains__(self, x)

    class CVList(list):
        """
        Class representation of all sets of resources and their nested
        annotation corresponding to a given qualifier. It have similar
        structure like that of a list but has only restricted type of
        entries (of type ExternalResources) within it

        Parameters
        ----------
        cvlist : list
            a list containing entries confirming to ExternalResources structure

        """

        def __init__(self, cvlist=[]):
            if not isinstance(cvlist, list):
                raise TypeError("The resources passed must be inside a list")
            for item in cvlist:
                if isinstance(item, CVTerm.ExternalResources):
                    list.append(self, item)
                elif isinstance(item, dict):
                    list.append(self, CVTerm.ExternalResources(item))
                else:
                    raise TypeError("All items must confirm to "
                                    "ExternalResources structure")

        def __len__(self):
            return list.__len__(self)

        def __delitem__(self, index):
            list.__delitem__(self, index)

        def insert(self, index, value):
            list.insert(self, index, CVTerm.ExternalResources(value))

        def append(self, value):
            list.append(self, CVTerm.ExternalResources(value))

        def __setitem__(self, index, value):
            list.__setitem__(self, index, CVTerm.ExternalResources(value))

        def __getitem__(self, index):
            return list.__getitem__(self, index)

    class ExternalResources(dict):
        """
        Class representation of a single set of resources and its nested
        annotation. Its a special type of dict with restricted keys and
        values

        Parameters
        ----------
        data : dict
            A dictionary containing the resources and nested annotation
            {
                "resources" : [],
                "nested_data" : CVTerm
             }

        Allowed Keys
        ----------
        "resources" : string
            for accessing the mapped resources
        "nested_data" : string
            for accessing the nested annotation data

        """

        ANNOTATION_KEYS = ['resources', 'nested_data']

        def __init__(self, data={}):
            self._resources = []
            self._nested_data = None
            if not isinstance(data, dict):
                raise TypeError("The value passed must be of type dict.")
            for key, value in data.items():
                if key not in self.ANNOTATION_KEYS:
                    raise ValueError("Key '%s' is not allowed. Only "
                                     "allowed keys are 'resource', "
                                     "'nested_data'." % key)
                if key == 'resources':
                    if not isinstance(value, list):
                        raise TypeError("Resources must be put in a list")
                    dict.__setitem__(self, key, value)
                if key == 'nested_data':
                    if isinstance(value, CVTerm):
                        dict.__setitem__(self, key, value)
                    elif isinstance(value, dict):
                        dict.__setitem__(self, key, CVTerm(value))
                    else:
                        raise TypeError("The nested data structure does "
                                        "not have valid CVTerm format")

        def __getitem__(self, key):
            if key not in self.ANNOTATION_KEYS:
                raise ValueError("Key %s is not allowed. Only allowed keys are"
                                 " 'resources', 'nested_data'." % key)
            return dict.__getitem__(self, key)

        def __setitem__(self, key, value):
            """Restricting the keys and values that can be set.
               Only allowed keys are 'resources' and 'nested_data'
            """
            if key not in self.ANNOTATION_KEYS:
                raise ValueError("Key %s is not allowed. Only allowed "
                                 "keys are 'resources', 'nested_data'."
                                 % key)
            if key == 'resources':
                if not isinstance(value, list):
                    raise TypeError("Resources must be put in a list")
                dict.__setitem__(self, key, value)
            elif key == 'nested_data':
                if isinstance(value, CVTerm):
                    dict.__setitem__(self, key, value)
                elif isinstance(value, dict):
                    dict.__setitem__(self, key, CVTerm(value))
                else:
                    raise TypeError("The value passed has invalid format.")
                dict.__setitem__(self, key, value)

        def __delitem__(self, key):
            dict.__delitem__(self, key)

        def __iter__(self):
            return dict.__iter__(self)

        def __len__(self):
            return dict.__len__(self)

        def __contains__(self, x):
            return dict.__contains__(self, x)
