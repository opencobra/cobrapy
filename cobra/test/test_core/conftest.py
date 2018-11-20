# -*- coding: utf-8 -*-

"""Module level fixtures"""

from __future__ import absolute_import

import pytest

from cobra.util.solver import solvers


solver_trials = ['glpk',
                 pytest.param('cplex',
                              marks=pytest.mark.skipif(
                                  'cplex' not in solvers,
                                  reason='No CPLEX found on PYTHONPATH')),
                 pytest.param('gurobi',
                              marks=pytest.mark.skipif(
                                  'gurobi' not in solvers,
                                  reason='No Gurobi found on PYTHONPATH'))]


@pytest.fixture(scope="function", params=solver_trials)
def solved_model(request, model):
    model.solver = request.param
    solution = model.optimize()
    return solution, model
