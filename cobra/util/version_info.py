# -*- coding: utf-8 -*-

# Adapated from:
# https://github.com/pandas-dev/pandas/blob/master/pandas/util/_print_versions.py
# which is published under a BSD license.

from __future__ import absolute_import, print_function

from builtins import dict

import platform

__all__ = ("show_versions",)

SYS_ORDER = [
    "OS",
    "OS-release",
    "Python"
]


def get_sys_info():
    """Returns system information as a dict."""
    blob = dict()
    blob["OS"] = platform.system()
    blob["OS-release"] = platform.release()
    blob["Python"] = platform.python_version()
    return blob


def show_versions():
    """Print the formatted information to standard out."""
    info = get_sys_info()
    format_str = "{:<%d} {:>%d}" % (max(map(len, info)),
                                    max(map(len, info.values())))
    print("\nSystem Information")
    print("==================")
    for name in SYS_ORDER:
        print(format_str.format(name, info[name]))
