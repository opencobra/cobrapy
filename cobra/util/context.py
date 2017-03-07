# -*- coding: utf-8 -*-

from __future__ import absolute_import

from functools import partial


class HistoryManager(object):
    """Record a list of actions to be taken at a later time. Used to
    implement context managers that allow temporary changes to a
    :class:`~cobra.core.Model`.

    """

    def __init__(self):

        # self._history just acts as a stack
        self._history = []

    def __call__(self, operation):
        """Add the corresponding method to the history stack.

        Parameters
        ----------
        operation : `function`
            A function to be called at a later time

        """

        self._history.append(operation)

    def reset(self):
        """Trigger executions for all items in the stack in reverse order"""
        while self._history:
            entry = self._history.pop()
            entry()


def get_context(obj):
    """Search for a context manager"""
    try:
        return obj._contexts[-1]
    except (AttributeError, IndexError):
        pass

    try:
        return obj._model._contexts[-1]
    except (AttributeError, IndexError):
        pass

    return None


def resettable(f):
    """A decorator to simplify the context management of simple object
    attributes. Gets the value of the attribute prior to setting it, and stores
    a function to set the value to the old value in the HistoryManager.
    """

    def wrapper(self, new_value):
        context = get_context(self)
        if context:
            old_value = getattr(self, f.__name__)
            # Don't clutter the context with unchanged variables
            if old_value == new_value:
                return
            context(partial(f, self, old_value))

        f(self, new_value)

    return wrapper
