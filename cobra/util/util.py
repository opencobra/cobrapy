# -*- coding: utf-8 -*-

from __future__ import absolute_import

import sympy
import logging
import operator
from collections import Mapping, OrderedDict

from six.moves import reduce


class Frozendict(dict):
    def __init__(self, iterable, **kwargs):
        super(Frozendict, self).__init__(iterable, **kwargs)

    def popitem(self):
        raise AttributeError("'Frozendict' object has no attribute 'popitem")

    def pop(self, k, d=None):
        raise AttributeError("'Frozendict' object has no attribute 'pop")

    def __setitem__(self, key, value):
        raise AttributeError(
            "'Frozendict' object has no attribute '__setitem__")

    def setdefault(self, k, d=None):
        raise AttributeError(
            "'Frozendict' object has no attribute 'setdefault")

    def __delitem__(self, key):
        raise AttributeError(
            "'Frozendict' object has no attribute '__delitem__")

    def __hash__(self):
        return hash(tuple(sorted(self.items())))

    def update(self, E=None, **F):
        raise AttributeError("'Frozendict' object has no attribute 'update")


class FrozenOrderedDict(Mapping):
    """
    Frozen OrderedDict.
    """

    def __init__(self, *args, **kwargs):
        self.__dict = OrderedDict(*args, **kwargs)
        self.__hash = None

    def __getitem__(self, item):
        return self.__dict[item]

    def __iter__(self):
        return iter(self.__dict)

    def __len__(self):
        return len(self.__dict)

    def __hash__(self):
        if self.__hash is None:
            self.__hash = reduce(operator.xor, map(hash, self.iteritems()), 0)

        return self.__hash

    def __repr__(self):
        return '{}({!r})'.format(self.__class__.__name__, self.items())

    def copy(self, *args, **kwargs):
        new_dict = self.__dict.copy()

        if args or kwargs:
            new_dict.update(OrderedDict(*args, **kwargs))

        return self.__class__(new_dict)


class AutoVivification(dict):
    """Implementation of perl's autovivification feature. Checkout
    http://stackoverflow.com/a/652284/280182 """

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
