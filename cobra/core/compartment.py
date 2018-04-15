# -*- coding: utf-8 -*-

from __future__ import absolute_import

from copy import deepcopy
from cobra.util.util import format_long_string

from cobra.core.object import Object


class Compartment(Object):
    """Compartment is a class for holding information regarding
        a compartment in a cobra.Model object

        Parameters
        ----------
        id : string
           An identifier for the compartment
        name : string
           A human readable name.
        """
    def __init__(self, id=None, name=None):
        Object.__init__(self, id, name)
        self._model = None

    def copy(self):
        return deepcopy(self)

    @property
    def model(self):
        return (self._model)

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
