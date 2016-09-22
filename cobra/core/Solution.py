from __future__ import absolute_import, print_function

from collections import OrderedDict

import time
import datetime

import cobra

try:
    import pandas
except ImportError:
    pandas = None

from cobra.exceptions import UndefinedSolution
from warnings import warn
import logging

logger = logging.getLogger(__name__)


class SolutionBase(object):
    def __new__(cls, *args, **kwargs):
        # this is a cobrapy compatibility hack
        if len(args) == 1 and not isinstance(args[0], cobra.core.Model):
            cobrapy_solution = super(SolutionBase, cls).__new__(LegacySolution)
            cobrapy_solution.__init__(*args, **kwargs)
            return cobrapy_solution
        else:
            return super(SolutionBase, cls).__new__(cls)

    def __init__(self, model, *args, **kwargs):
        self.f = None
        self.model = model
        self._x = None
        self._y = None
        self._x_dict = None
        self._y_dict = None

    def get_primal_by_id(self, reaction_id):
        """Return a flux/primal value for a reaction.

        Parameters
        ----------
        reaction_id : str
            A reaction ID.
        """
        return self.x_dict[reaction_id]

    @property
    def x_dict(self):
        warn("use solution.fluxes instead", DeprecationWarning)
        if self._x_dict is None:
            return self.fluxes
        else:
            return self._x_dict

    @x_dict.setter
    def x_dict(self, value):
        warn("not used", DeprecationWarning)
        self._x_dict = value

    @property
    def x(self):
        warn("use solution.fluxes.values instead", DeprecationWarning)
        if self._x is None:
            return self.fluxes.values()
        else:
            return self._x

    @x.setter
    def x(self, value):
        warn("not used", DeprecationWarning)
        self._x = value

    @property
    def y_dict(self):
        warn("use solution.reduced_costs instead", DeprecationWarning)
        if self._y_dict is None:
            return self.reduced_costs
        else:
            return self._y_dict

    @y_dict.setter
    def y_dict(self, value):
        self._y_dict = value

    @property
    def y(self):
        warn("use solution.reduced_costs.values instead", DeprecationWarning)
        if self._y is None:
            return self.reduced_costs.values()
        else:
            return self._y

    @y.setter
    def y(self, value):
        warn("not used", DeprecationWarning)
        self._y = value

    @property
    def objective_value(self):
        return self.f

    def __repr__(self):
        if self.f is None:
            return "<Solution '%s' at 0x%x>" % (self.status, id(self))
        return "<Solution %.2f at 0x%x>" % (self.f, id(self))


class Solution(SolutionBase):
    """Stores the solution from optimizing a cobra.Model. This is
    used to provide a single interface to results from different
    solvers that store their values in different ways.

    Attributes
    ----------
    fluxes : OrderedDict
        A dictionary of flux values.
    reduced_costs : OrderedDict
        A dictionary of reduced costs.

    Notes
    -----
    See also documentation for cobra.core.Solution.Solution for an extensive
    list of inherited attributes.
    """

    def __init__(self, model, *args, **kwargs):
        """
        Parameters
        ----------
        model : SolverBasedModel
        """
        super(Solution, self).__init__(model, *args, **kwargs)
        self.f = model.solver.objective.value
        self.fluxes = OrderedDict()
        self.shadow_prices = OrderedDict()
        self.reduced_costs = OrderedDict()
        for reaction in model.reactions:
            self.fluxes[reaction.id] = reaction.flux
            self.reduced_costs[reaction.id] = reaction.reduced_cost
        for metabolite in model.metabolites:
            self.shadow_prices[metabolite.id] = self.model.solver.constraints[
                metabolite.id].dual
        self.status = model.solver.status
        self._reaction_ids = [r.id for r in self.model.reactions]
        self._metabolite_ids = [m.id for m in self.model.metabolites]

    def __dir__(self):
        # Hide deprecated attributes and methods from user.
        fields = sorted(dir(type(self)) + list(self.__dict__.keys()))
        fields.remove('x')
        fields.remove('y')
        fields.remove('x_dict')
        fields.remove('y_dict')
        return fields

    def _repr_html_(self):
        return "%s: %f" % (self.model.solver.objective.expression, self.f)


class LazySolution(SolutionBase):
    """This class implements a lazy evaluating version of the cobrapy
    Solution class.

    Attributes
    ----------
    model : SolverBasedModel
    fluxes : OrderedDict
        A dictionary of flux values.
    reduced_costs : OrderedDict
        A dictionary of reduced costs.

    Notes
    -----
    See also documentation for cobra.core.Solution.Solution for an extensive
    list of inherited attributes.

    """

    def __init__(self, model, *args, **kwargs):
        """
        Parameters
        ----------
        model : SolverBasedModel
        """
        super(LazySolution, self).__init__(model, *args, **kwargs)
        if self.model._timestamp_last_optimization is not None:
            self._time_stamp = self.model._timestamp_last_optimization
        else:
            self._time_stamp = time.time()
        self._f = None

    @property
    def data_frame(self):
        if pandas:
            return pandas.DataFrame(
                {'fluxes': pandas.Series(self.fluxes),
                 'reduced_costs': pandas.Series(self.reduced_costs)})
        else:
            warn("pandas not available")

    def _repr_html_(self):
        if pandas:
            return self.data_frame._repr_html_()
        else:
            warn("pandas not available")

    def _check_freshness(self):
        """Raises an exceptions if the solution might have become invalid
        due to re-optimization of the attached model.

        Raises
        ------
        UndefinedSolution
            If solution has become invalid.
        """
        if self._time_stamp != self.model._timestamp_last_optimization:
            def timestamp_formatter(timestamp):
                datetime.datetime.fromtimestamp(timestamp).strftime(
                    "%Y-%m-%d %H:%M:%S:%f")

            raise UndefinedSolution(
                'The solution (captured around %s) has become invalid as the '
                'model has been re-optimized recently (%s).' % (
                    timestamp_formatter(self._time_stamp),
                    timestamp_formatter(
                        self.model._timestamp_last_optimization))
            )

    @property
    def status(self):
        self._check_freshness()
        return self.model.solver.status

    @property
    def f(self):
        self._check_freshness()
        if self._f is None:
            return self.model.solver.objective.value
        else:
            return self._f

    @f.setter
    def f(self, value):
        self._f = value

    @property
    def fluxes(self):
        self._check_freshness()
        primals = OrderedDict()
        for reaction in self.model.reactions:
            primals[reaction.id] = reaction.flux
        return primals

    @property
    def reduced_costs(self):
        self._check_freshness()
        duals = OrderedDict()
        for reaction in self.model.reactions:
            duals[reaction.id] = reaction.reduced_cost
        return duals

    @property
    def shadow_prices(self):
        self._check_freshness()
        duals = OrderedDict()
        for metabolite in self.model.metabolites:
            duals[metabolite.id] = self.model.solver.constraints[
                metabolite.id].dual
        return duals

    def get_primal_by_id(self, reaction_id):
        """Return a flux/primal value for a reaction.

        Parameters
        ----------
        reaction_id : str
            A reaction ID.
        """
        self._check_freshness()
        return self.model.reactions.get_by_id(reaction_id).flux


class LegacySolution(object):
    """Stores the solution from optimizing a cobra.Model. This is
    used to provide a single interface to results from different
    solvers that store their values in different ways.

    f: The objective value

    solver: A string indicating which solver package was used.

    x: List or Array of the values from the primal.

    x_dict: A dictionary of reaction ids that maps to the primal values.

    y: List or Array of the values from the dual.

    y_dict: A dictionary of reaction ids that maps to the dual values.

    """

    def __init__(self, f, x=None,
                 x_dict=None, y=None, y_dict=None,
                 solver=None, the_time=0, status='NA'):
        self.solver = solver
        self.f = f
        self.x = x
        self.x_dict = x_dict
        self.status = status
        self.y = y
        self.y_dict = y_dict

    def dress_results(self, model):
        """.. warning :: deprecated"""
        warn("unnecessary to call this deprecated function",
             DeprecationWarning)
