# -*- coding: utf-8 -*-

from __future__ import absolute_import

from depinfo import print_dependencies
from six import string_types


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


def show_versions():
    """Print dependency information."""
    print_dependencies("cobra")


def is_not_sane(id):
    """Check if a id is sane to be used for cobra components."""
    return not isinstance(id, string_types) or \
        len(id) < 1 or any(c.isspace() for c in id)
