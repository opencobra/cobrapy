# -*- coding: utf-8 -*-

"""Provide a class for compartments."""

from __future__ import absolute_import

from copy import deepcopy
from cobra.util import format_long_string, is_not_sane
from six import string_types

from cobra.core.object import Object


class Compartment(Object):
    """
    Compartment is a class for holding information regarding
    a compartment in a cobra.Model object

    Parameters
    ----------
    id : string
       An identifier for the compartment
    name : string
       A human readable name.

    """
    def __init__(self, id=None, name=""):
        super(Compartment, self).__init__(id=id, name=name)
        self._id = None
        self.id = id

    def __contains__(self, metabolite):
        return metabolite.compartment is self

    def __eq__(self, other):
        if self is other:
            return True
        if isinstance(other, string_types):
            return self._id == other
        if isinstance(other, Compartment):
            return self._id == other.id
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self._id)

    def copy(self):
        return deepcopy(self)

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, value):
        if is_not_sane(value):
            raise TypeError("The compartment ID must be a non-empty string")
        self._id = value

    def _repr_html_(self):
        return """
        <table>
            <tr>
                <td><strong>Compartment identifier</strong></td><td>{id}</td>
            </tr><tr>
                <td><strong>Name</strong></td><td>{name}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{address}</td>
            </tr>
        </table>""".format(id=self.id, name=format_long_string(self.name),
                           address='0x0%x' % id(self))
