from collections import OrderedDict
from uuid import uuid1
from functools import partial


class HistoryManager(object):
    """Record a list of actions to be taken at a later time. Used to
    implement context managers that allow temporary changes to a
    :class:`~cobra.core.Model`.

    """

    def __init__(self):
        self._history = OrderedDict()

    def __call__(self, operation, bookmark=None):
        """ Add the corresponding method to the history stack.

        Parameters
        ----------
        operation: `function`
            A function to be called at a later time
        bookmark: string or None
            The index of the operation in `_history`. If the given index is
            already present, the previous entry will be removed.

        Returns
        -------
        str: the newly created uuid

        """

        if bookmark is None:
            entry_id = uuid1()
        else:
            entry_id = bookmark

        # make sure that entry is added to the end of history
        self._history.pop(entry_id, None)
        self._history[entry_id] = operation

        return entry_id

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.reset()

    def execute(self, bookmark=None):
        """ Execute up to the given entry in the history.
        If no entry is specified, executes the last entry.

        Parameters
        ----------
        bookmark: `str`
            The index of the operation to execute.

        ..NOTE: I'm not sure we want this to work this way. I'm wondering if
        bookmark should instead check to see if that particular attribute has
        been previously set, and if so, don't add this undo function.
        """
        if bookmark is None:
            try:
                (uuid, entry) = self._history.popitem()
                entry()
            except KeyError:  # history is empty
                pass

        elif bookmark in list(self._history.keys()):
            uuid = False
            while uuid is not bookmark:
                (uuid, entry) = self._history.popitem()
                entry()
        else:
            raise Exception('Provided bookmark %s cannot be found.')

    def reset(self):
        """Trigger executions for all items in the stack in reverse order"""
        if self._history:  # history is not empty
            self.execute(bookmark=list(self._history.keys())[0])

    @property
    def last_uuid(self):
        """Return the key for the most recent item in history"""
        return next(reversed(self._history))


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
