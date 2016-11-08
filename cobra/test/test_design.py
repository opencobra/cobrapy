import pytest

from cobra.design import *
from cobra.design.design_algorithms import _add_decision_variable
from cobra.solvers import get_solver_name
from .conftest import model

try:
    solver = get_solver_name(mip=True)
except ImportError:
    no_mip_solver = True
else:
    no_mip_solver = False


class TestDesignAlgorithms:
    """Test functions in cobra.design"""
    def test_dual(self, model):
        assert abs(model.optimize("maximize").f - 0.874) < 0.001
        dual = dual_problem(model)
        assert abs(dual.optimize("minimize").f - 0.874) < 0.001

    def test_dual_integer_vars_as_lp(self, model):
        var = _add_decision_variable(model, "AKGDH")
        assert abs(model.optimize("maximize").f - 0.874) < 0.001
        # as lp: make integer continuous, set to 1
        dual = dual_problem(model, "maximize", [var.id], copy=True)
        r = dual.reactions.get_by_id(var.id)
        r.variable_kind = "continuous"
        r.lower_bound = r.upper_bound = 1
        assert abs(dual.optimize("minimize").f - 0.874) < 0.001
        r.lower_bound = r.upper_bound = 0
        assert abs(dual.optimize("minimize").f - 0.858) < 0.001

    @pytest.mark.skipif(no_mip_solver, reason="no MILP solver found")
    def test_dual_integer_vars_as_mip(self, model):
        # mip
        var = _add_decision_variable(model, "AKGDH")
        dual = dual_problem(model, "maximize", [var.id], copy=True)
        var_in_dual = dual.reactions.get_by_id(var.id)
        # minimization, so the optimal value state is to turn off AKGDH
        assert abs(dual.optimize("minimize").f - 0.858) < 0.001
        # turn off AKGDH in dual
        var_in_dual.lower_bound = var_in_dual.upper_bound = 1
        assert abs(dual.optimize("minimize").f - 0.874) < 0.001
        # turn on AKGDH in dual
        var_in_dual.lower_bound = var_in_dual.upper_bound = 0
        assert abs(dual.optimize("minimize").f - 0.858) < 0.001

    @pytest.mark.skipif(no_mip_solver, reason="no MILP solver found")
    def test_optknock(self, model):
        model.reactions.get_by_id("EX_o2_e").lower_bound = 0
        knockable_reactions = ["ACKr", "AKGDH", "ACALD", "LDH_D"]
        optknock_problem = set_up_optknock(model, "EX_lac__D_e",
                                           knockable_reactions, n_knockouts=2,
                                           copy=False)
        solution = run_optknock(optknock_problem, tolerance_integer=1e-9)
        assert "ACKr" in solution.knockouts
        assert "ACALD" in solution.knockouts
        assert abs(solution.f - 17.891) < 0.001
