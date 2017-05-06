# -*- coding: utf-8 -*-

from __future__ import absolute_import


class AutoVivification(dict):
    """Implementation of perl's autovivification feature. Checkout
    http://stackoverflow.com/a/652284/280182 """

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
