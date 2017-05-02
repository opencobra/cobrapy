# -*- coding: utf-8 -*-

from __future__ import absolute_import

import sympy
import logging
import operator
from collections import Mapping

from six import iteritems


class FrozenDict(Mapping):
    def __init__(self, *args, **kwargs):
        super(FrozenDict, self).__init__()
        self.__dict = dict(*args, **kwargs)
        self.__hash = None
        self.__items = None

    def __getitem__(self, item):
        return self.__dict[item]

    def __iter__(self):
        return iter(self.__dict)

    def __len__(self):
        return len(self.__dict)

    def __hash__(self):
        if self.__items is None:
            self.__items = sorted(iteritems(self.__dict))
        if self.__hash is None:
            self.__hash = hash(self.__items)
        return self.__hash

    def __repr__(self):
        if self.__items is None:
            self.__items = sorted(iteritems(self.__dict))
        return '{}({!r})'.format(self.__class__.__name__, self.__items)

    def copy(self, *args, **kwargs):
        new_dict = self.__dict.copy()

        if args or kwargs:
            new_dict.update(dict(*args, **kwargs))

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
