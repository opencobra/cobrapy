# -*- coding: utf-8 -*-

from __future__ import absolute_import


def format_long_string(string, max_length=50):
    if len(string) > max_length:
        string = string[:max_length - 3]
        string += '...'
    return string


class AutoVivification(dict):
    """Implementation of perl's autovivification feature. Checkout
    http://stackoverflow.com/a/652284/280182 """

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
