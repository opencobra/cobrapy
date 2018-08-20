# -*- coding: utf-8 -*-

"""Define module level fixtures."""


from __future__ import absolute_import

import pytest
from pandas import Series

import cobra.util.solver as sutil
from cobra.core import Metabolite, Model, Reaction, Solution

# The scipy interface is currently unstable and may yield errors or infeasible
# solutions.


@pytest.fixture(scope="session",
                params=[s for s in ["glpk", "cplex", "gurobi"]
                        if s in sutil.solvers])
def all_solvers(request):
    """Return the avaialble solvers."""
    return request.param


@pytest.fixture(scope="session",
                params=[s for s in ["cplex", "gurobi"]
                        if s in sutil.solvers])
def qp_solvers(request):
    """Return the available QP solvers."""
    return request.param


@pytest.fixture(scope="module")
def room_model():
    """
    Generate ROOM model as described in [1]_

    References
    ----------
    .. [1] Tomer Shlomi, Omer Berkman and Eytan Ruppin, "Regulatory on/off
     minimization of metabolic flux changes after genetic perturbations",
     PNAS 2005 102 (21) 7695-7700; doi:10.1073/pnas.0406346102

    """
    test_model = Model("papin_2003")

    v_1 = Reaction("v1")
    v_2 = Reaction("v2")
    v_3 = Reaction("v3")
    v_4 = Reaction("v4")
    v_5 = Reaction("v5")
    v_6 = Reaction("v6", upper_bound=0.0)
    b_1 = Reaction("b1", upper_bound=10.0, lower_bound=0.0)
    b_2 = Reaction("b2")
    b_3 = Reaction("b3")

    test_model.add_reactions([v_1, v_2, v_3, v_4, v_5, v_6, b_1, b_2, b_3])

    v_1.reaction = "A -> B"
    v_2.reaction = "2 B -> C + byp"
    v_3.reaction = "2 B + cof -> D"
    v_4.reaction = "D -> E + cof"
    v_5.reaction = "C + cof -> D"
    v_6.reaction = "C -> E"
    b_1.reaction = "-> A"
    b_2.reaction = "E ->"
    b_3.reaction = "byp ->"

    test_model.objective = 'b2'

    return test_model


@pytest.fixture(scope="module")
def room_solution():
    """
    Generate ROOM solution as described in [1]_

    References
    ----------
    .. [1] Tomer Shlomi, Omer Berkman and Eytan Ruppin, "Regulatory on/off
     minimization of metabolic flux changes after genetic perturbations",
     PNAS 2005 102 (21) 7695-7700; doi:10.1073/pnas.0406346102

    """
    fluxes = Series({'b1': 10.0, 'b2': 5.0, 'b3': 5.0, 'v1': 10.0, 'v2': 5.0,
                     'v3': 0.0, 'v4': 0.0, 'v5': 0.0, 'v6': 5.0})
    reduced_costs = Series({'b1': 0.0, 'b2': 0.0, 'b3': 0.0, 'v1': 0.0,
                            'v2': 0.0, 'v3': 0.0, 'v4': 0.0, 'v5': 0.0,
                            'v6': 0.0})
    shadow_prices = Series({'b1': 0.0, 'b2': 0.0, 'b3': 0.0, 'v1': 0.0,
                            'v2': 0.0, 'v3': 0.0, 'v4': 0.0, 'v5': 0.0,
                            'v6': 0.0})
    sol = Solution(objective_value=5.000, status='optimal',
                   fluxes=fluxes,
                   reduced_costs=reduced_costs,
                   shadow_prices=shadow_prices)
    return sol
