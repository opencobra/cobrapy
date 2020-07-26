# -*- coding: utf-8 -*-

"""Define the singleton meta class."""

from __future__ import absolute_import


class Singleton(type):
    """Implementation of the singleton pattern as a meta class."""

    _instances = {}

    def __call__(cls, *args, **kwargs):
        """Override an inheriting class' call."""
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args,
                                                                 **kwargs)
        return cls._instances[cls]
