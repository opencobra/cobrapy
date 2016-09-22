# Copyright 2013 Novo Nordisk Foundation Center for Biosustainability, DTU.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import absolute_import, print_function

import colorsys
import itertools
import logging
import platform
import re
from collections import OrderedDict
from datetime import datetime
from functools import partial
import inspect
from itertools import islice
from time import time
from uuid import uuid1

import pip
import six
from six.moves import range

logger = logging.getLogger(__name__)


class frozendict(dict):
    def __init__(self, iterable, **kwargs):
        super(frozendict, self).__init__(iterable, **kwargs)

    def popitem(self):
        raise AttributeError("'frozendict' object has no attribute 'popitem")

    def pop(self, k, d=None):
        raise AttributeError("'frozendict' object has no attribute 'pop")

    def __setitem__(self, key, value):
        raise AttributeError(
            "'frozendict' object has no attribute '__setitem__")

    def setdefault(self, k, d=None):
        raise AttributeError(
            "'frozendict' object has no attribute 'setdefault")

    def __delitem__(self, key):
        raise AttributeError(
            "'frozendict' object has no attribute '__delitem__")

    def __hash__(self):
        return hash(tuple(sorted(self.items())))

    def update(self, E=None, **F):
        raise AttributeError("'frozendict' object has no attribute 'update")


class ProblemCache(object):
    """
    Variable and constraint cache for models.

    To be used in complex methods that require many extra variables when one
    must run multiple simulations.

    It allows rollback to the previous state in case one iteration fails to
    build the problem or generates an invalid state.

    """

    def __init__(self, model):
        self.time_machine = None
        self._model = model
        self.variables = {}
        self.constraints = {}
        self.objective = None
        self.original_objective = model.objective
        self.time_machine = TimeMachine()
        self.transaction_id = None

    def begin_transaction(self):
        """
        Creates a time point. If rollback is called, the variables and
        constrains will be reverted to this point.
        """
        self.transaction_id = uuid1()
        self.time_machine(do=int, undo=int, bookmark=self.transaction_id)

    @property
    def model(self):
        return self._model

    def _append_constraint(self, constraint_id, create, *args, **kwargs):
        self.constraints[constraint_id] = create(self.model, constraint_id,
                                                 *args, **kwargs)
        self._model.solver._add_constraint(self.constraints[constraint_id])

    def _remove_constraint(self, constraint_id):
        constraint = self.constraints.pop(constraint_id)
        self.model.solver._remove_constraint(constraint)

    def _append_variable(self, variable_id, create, *args, **kwargs):
        self.variables[variable_id] = create(self.model, variable_id, *args,
                                             **kwargs)
        self.model.solver._add_variable(self.variables[variable_id])

    def _remove_variable(self, variable_id):
        variable = self.variables.pop(variable_id)
        self.model.solver._remove_variable(variable)

    def _rebuild_variable(self, variable):
        (type, lb, ub,
         name) = variable.type, variable.lb, variable.ub, variable.name

        def rebuild():
            self.model.solver._remove_variable(variable)
            new_variable = self.model.solver.interface.Variable(name, lb=lb,
                                                                ub=ub,
                                                                type=type)
            self.variables[name] = variable
            self.model.solver._add_variable(new_variable, sloppy=True)

        return rebuild

    def add_constraint(self, constraint_id, create, update, *args, **kwargs):
        """
        Adds a new cached constraint.

        The create and update functions must have the following signatures:
        >>> create(model, constraint_id, *args)
        >>> update(model, constraint, *args)

        "args" in the first example must match args on the second example.

        Arguments
        ---------
        constraint_id: str
            The identifier of the constraint
        create: function
            A function that creates an optlang.interface.Constraint
        update: function
            a function that updates an optlang.interface.Constraint
        """
        if constraint_id in self.constraints:
            if update is not None:
                update(self.model, self.constraints[constraint_id], *args,
                       **kwargs)
        else:
            self.time_machine(
                do=partial(self._append_constraint, constraint_id, create,
                           *args, **kwargs),
                undo=partial(self._remove_constraint, constraint_id)
            )

    def add_variable(self, variable_id, create, update, *args, **kwargs):
        """
        Adds a new cached variable.

        The create and update functions must have the following signatures:
        >>> create(model, variable_id, *args)
        >>> update(model, variable, *args)

        "args" in the first example must match args on the second example.

        Arguments
        ---------
        constraint_id: str
            The identifier of the constraint
        create: function
            A function that creates an optlang.interface.Variable
        update: function
            a function that updates an optlang.interface.Variable
        """
        if variable_id in self.variables:
            if update is not None:
                self.time_machine(
                    do=partial(update, self.model, self.variables[variable_id],
                               *args, **kwargs),
                    undo=self._rebuild_variable(self.variables[variable_id]))
        else:
            self.time_machine(
                do=partial(self._append_variable, variable_id, create, *args,
                           **kwargs),
                undo=partial(self._remove_variable, variable_id)
            )

    def add_objective(self, create, update, *args):
        if self.objective is None:
            self.objective = create(self.model, *args)
            self.time_machine(
                do=partial(setattr, self.model, 'objective', self.objective),
                undo=partial(setattr, self.model, 'objective',
                             self.model.objective)
            )
        else:
            if update:
                self.objective = update(self.model, *args)
                self.time_machine(
                    do=partial(setattr, self.model, 'objective',
                               self.objective),
                    undo=partial(setattr, self.model, 'objective',
                                 self.model.objective)
                )

    def reset(self):
        """
        Removes all constraints and variables from the cache.
        """
        self.model.solver._remove_constraints(self.constraints.values())
        self.model.solver._remove_variables(self.variables.values())
        self.model.objective = self.original_objective
        self.variables = {}
        self.objective = None
        self.constraints = {}
        self.transaction_id = None
        self.time_machine.history.clear()

    def rollback(self):
        """
        Returns to the previous transaction start point.
        """
        if self.transaction_id is None:
            raise RuntimeError(
                "Start transaction must be called before rollback")
        self.time_machine.undo(self.transaction_id)
        self.transaction_id = None

    def __enter__(self):
        """
        Allows problem cache to be used with a _with_ statement.

        Examples
        --------
        You want to run room/lmoma for every single knockout.
        >>> with ProblemCache(model) as cache:
        >>>     for reaction in reactions:
        >>>         result = lmoma(model, reference=reference, cache=cache)

        Returns
        -------
        ProblemCache
            returns itself
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.reset()


class Singleton(object):
    """
    Singleton class to be extended
    """
    _instance = None

    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(Singleton, cls).__new__(cls)
        return cls._instance


class AutoVivification(dict):
    """Implementation of perl's autovivification feature. Checkout
    http://stackoverflow.com/a/652284/280182 """

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


class TimeMachine(object):
    """Travel back and forth in time."""

    def __init__(self):
        super(TimeMachine, self).__init__()
        self.history = OrderedDict()

    def __call__(self, do=None, undo=None, bookmark=None):
        output = do()
        current_time = time()
        if bookmark is None:
            entry_id = uuid1()
        else:
            entry_id = bookmark
        # make sure that entry is added to the end of history
        self.history.pop(entry_id, None)
        self.history[entry_id] = {'unix_epoch': current_time, 'undo': undo,
                                  'redo': do}
        return entry_id, output

    def __str__(self):
        info = '\n'
        for item in six.iteritems(self.history):
            info += self._history_item_to_str(item)
        return info

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.reset()

    @staticmethod
    def _history_item_to_str(item):
        info = ''
        uuid, entry = item
        info += datetime.fromtimestamp(entry['unix_epoch']).strftime(
            '%Y-%m-%d %H:%M:%S') + '\n'
        undo_entry = entry['undo']
        try:
            # partial  (if .keywords is None print {} instead)
            elements = undo_entry.func, undo_entry.args, \
                       undo_entry.keywords or {}
            info += 'undo: ' + ' '.join(
                [str(elem) for elem in elements]) + '\n'
        except AttributeError:  # normal python function
            info += 'undo: ' + undo_entry.__name__ + '\n'

        redo_entry = entry['redo']
        try:
            # partial
            elements = redo_entry.func, redo_entry.args, \
                       redo_entry.keywords or {}
            info += 'redo: ' + ' '.join(
                [str(elem) for elem in elements]) + '\n'
        except AttributeError:
            info += 'redo: ' + redo_entry.__name__ + '\n'
        return info

    def undo(self, bookmark=None):
        if bookmark is None:
            try:
                (uuid, entry) = self.history.popitem()
                entry['undo']()
            except KeyError:  # history is empty
                pass
        elif bookmark in list(self.history.keys()):
            uuid = False
            while uuid is not bookmark:
                (uuid, entry) = self.history.popitem()
                entry['undo']()
        else:
            raise Exception(
                'Provided bookmark %s cannot be found in the time machine.')

    def redo(self):
        raise NotImplementedError

    def reset(self):
        if self.history:  # history is not empty
            self.undo(bookmark=list(self.history.keys())[0])


class Timer(object):
    """Taken from http://stackoverflow.com/a/5849861/280182"""

    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print('[%s]' % self.name, end=' ')
        print('Elapsed: %s' % (time() - self.tstart))


class IntelliContainer(object):
    def __init__(self, **kwargs):
        self._dict = dict(**kwargs)

    def __getattr__(self, value):
        return self._dict.get(value)

    def __setitem__(self, key, value):
        self._dict[key] = value

    def __iter__(self):
        return six.itervalues(self._dict)

    def __dir__(self):
        return list(self._dict.keys())


def inheritdocstring(name, bases, attrs):
    """Use as metaclass to inherit class and method docstrings from parent.

    Adapted from http://stackoverflow.com/questions/13937500/inherit-a
    -parent-class-docstring-as-doc-attribute """
    temp = type('temporaryclass', bases, {})
    if '__doc__' not in attrs or not attrs["__doc__"]:
        # create a temporary 'parent' to (greatly) simplify the MRO search
        for cls in inspect.getmro(temp):
            if cls.__doc__ is not None:
                attrs['__doc__'] = cls.__doc__
                break

    for attr_name, attr in attrs.items():
        if not attr.__doc__:
            for cls in inspect.getmro(temp):
                try:
                    if getattr(cls, attr_name).__doc__ is not None:
                        attr.__doc__ = getattr(cls, attr_name).__doc__
                        break
                except (AttributeError, TypeError):
                    continue

    return type(name, bases, attrs)


def partition_(lst, n):
    """Partition a list into n bite size chunks."""
    division = len(lst) / float(n)
    return [lst[int(round(division * i)): int(round(division * (i + 1)))] for i
            in range(n)]


def partition(ite, n):
    """Partition an iterable into n bite size chunks."""
    try:
        length = len(ite)
    except TypeError:
        ite = list(ite)
        length = len(ite)
    division = length / float(n)
    iterator = iter(ite)
    return [list(
        islice(iterator, 0, round(division * (i + 1)) - round(division * i)))
            for i in range(n)]


def flatten(l):
    return [item for sublist in l for item in sublist]


def generate_colors(n):
    hsv_tuples = [(v * 1.0 / n, 0.5, 0.5) for v in range(n)]
    color_map = {}
    for i in range(n):
        rgb = colorsys.hsv_to_rgb(*hsv_tuples[i])
        color = tuple(int(channel * 256) for channel in rgb)
        color_map[i] = '#%02x%02x%02x' % color
    return color_map


def memoize(function, memo={}):
    def wrapper(*args):
        if args in memo:
            return memo[args]
        else:
            rv = function(*args)
            memo[args] = rv
            return rv

    return wrapper


def get_system_info():
    # pip freeze (adapted from http://stackoverflow.com/a/24322465/280182)
    package_info = list()
    for dist in pip.get_installed_distributions():
        req = str(dist.as_requirement())
        package_info.append(req)
    return dict(package_info=package_info,
                platform=platform.platform(),
                machine=platform.machine(),
                system=platform.system())


def str_to_valid_variable_name(s):
    """Adapted from http://stackoverflow.com/a/3303361/280182"""

    # Remove invalid characters
    s = re.sub('[^0-9a-zA-Z_]', '_', s)

    # Remove leading characters until we find a letter or underscore
    s = re.sub('^[^a-zA-Z_]+', '', s)

    return s


def zip_repeat(long_iter, short_iter):
    """
    Zips two iterable objects but repeats the second one if it is shorter
    than the first one.

    Parameters
    ----------
    long_iter: iterable

    short_iter: iterable

    Returns
    -------
    generator
    """
    for i, j in zip(long_iter, itertools.cycle(short_iter)):
        yield i, j
