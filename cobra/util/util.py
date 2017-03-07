# -*- coding: utf-8 -*-

from __future__ import absolute_import, print_function

import logging

logger = logging.getLogger(__name__)


def _is_positive(n):
    """Robustly test if n is positive, yielding True on Exceptions"""
    try:
        if n >= 0:
            return True
        else:
            return False
    except Exception:
        return True


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


class AutoVivification(dict):
    """Implementation of perl's autovivification feature. Checkout
    http://stackoverflow.com/a/652284/280182 """

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
