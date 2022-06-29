"""
Some documentation:

Reactions: v1; v2; v3;
Stochiometric Contraints (metabolic nework, SBML Reactions):
v1 - v2 = 0

UserDefinedConstraints:
v1 + v2 < 2

Non-constant parameters (x1: not reactions; )

expression =: v1 + v2 + 1.2 x1 < 2
v1 + v2 + 1.2 x1*x1 < 2

-> examples:
-> optimize: (v1, v2, x1); how to access the optimal solution?

Test Cases:
------------
1. add user constraint to simple model
A -v1-> B -v2-> C
optimize: v2; bounds v1[0; 10]; v2[0, 5]
optimal value: 5;

UserConstraint: v1 + v2 <= 4
optimal value: 2;

2. only user constraint;
empty cobra model

x1 <= 2
maximize: x1
flux distribution optimal: x1=2
"""

from os.path import join

import pytest
from optlang import Variable, Constraint

from cobra.core.model import Model
from cobra.io import load_json_model, read_sbml_model, save_json_model


def ex_model(data_directory):
    model = read_sbml_model(join(data_directory, "userConstraint.xml"))
    return model


def test_user_defined_constraints(data_directory):
    cons_model = ex_model(data_directory)
    solution1 = cons_model.optimize()
    assert solution1.objective_value == pytest.approx(5.0)

    v1 = Variable('v1')
    v2 = Variable('v2')
    c1 = Constraint(v1 + v2, lb=0, ub=4, name='c1')
    cons_model.add_user_defined_constraints(c1)

    solution2 = cons_model.optimize()
    assert solution2.objective_value == pytest.approx(2.00)


def test_user_defined_constraints_on_single_flux(data_directory):
    cons_model = ex_model(data_directory)
    solution1 = cons_model.optimize()
    assert solution1.objective_value == pytest.approx(5.0)

    v2 = Variable('v2')
    const_1 = Constraint(v2, lb=0, ub=3, name='const_1')
    cons_model.add_user_defined_constraints([const_1])
    solution2 = cons_model.optimize()
    assert solution2.objective_value == pytest.approx(3.00)


def test_user_defined_constraints_on_single_variable():
    # an empty model
    model = Model("model_abc")
    cc1 = Variable("new_var")
    c1 = Constraint(cc1, lb=0, ub=2, name='c1')
    model.add_user_defined_constraints([c1])
    model.objective = model.variables.new_var
    solution = model.optimize()
    assert solution.objective_value == pytest.approx(2.00)
    # save json
    # read json, make sure it is till objective


def test_json_reading_writing(model, tmp_path):
    cc1 = Variable("FBA")
    cc2 = Variable('NH4t')
    cc3 = Variable('difference')
    c1 = Constraint(cc1 - cc2 - cc3, lb=0, ub=0, name='c1')

    c2 = Constraint(cc1 + cc2, lb=0, ub=10, name='c2')
    model.add_user_defined_constraints([c1, c2])
    solution1 = model.optimize()

    path_to_json = join(tmp_path, "userConstraint.json")
    save_json_model(model, path_to_json)

    model = load_json_model(path_to_json)
    assert model.user_defined_const is not None
    assert len(model.user_defined_const) == 2
    const_1 = model.user_defined_const[0]
    assert const_1.name == "c1"
    assert len(const_1.variables) == 5
    variable_names = {var.name for var in const_1.variables}
    assert 'FBA' in variable_names
    solution2 = model.optimize()
    assert solution1 == pytest.approx(solution2)


def test_user_defined_constraints_read_write_json(data_directory, tmp_path):
    cons_model = ex_model(data_directory)
    solution1 = cons_model.optimize()
    assert solution1.objective_value == pytest.approx(5.0)

    v1 = Variable('v1')
    v2 = Variable('v2')
    c1 = Constraint(v1 + v2, lb=0, ub=4, name='c1')
    cons_model.add_user_defined_constraints(c1)

    solution2 = cons_model.optimize()
    assert solution2.objective_value == pytest.approx(2.00)
    path_to_json = join(tmp_path, "userConstraint.json")
    save_json_model(cons_model, path_to_json)

    cons_model2 = load_json_model(path_to_json)
    solution3 = cons_model2.optimize()
    assert solution3.objective_value == pytest.approx(2.00)
    assert solution3.objective_value == pytest.approx(solution2.objective_value)


def test_user_defined_constraints_on_single_variable_read_write_json(tmp_path):
    # an empty model
    model = Model("model_abc")
    cc1 = Variable("new_var")
    c1 = Constraint(cc1, lb=0, ub=2, name='c1')
    model.add_user_defined_constraints([c1])
    model.objective = model.variables.new_var
    solution = model.optimize()
    assert solution.objective_value == pytest.approx(2.00)
    path_to_json = join(tmp_path, "userConstraint.json")
    save_json_model(model, path_to_json)

    model2 = load_json_model(path_to_json)
    solution2 = model2.optimize()
    assert solution2.objective_value == pytest.approx(2.00)
    assert solution2.objective_value == pytest.approx(solution.objective_value)


def test_user_defined_constraints_documented(model):
    solution1 = model.optimize()
    assert solution1.objective_value == pytest.approx(0.87392, 0.0001)

    FBA = Variable("FBA")
    NH4t = Variable("NH4t")
    c1 = Constraint(FBA - NH4t, lb=0, ub=0, name='c1')
    model.add_user_defined_constraints([c1])
    solution2 = model.optimize()
    assert solution2.fluxes["FBA"] == pytest.approx(4.66274, 0.0001)
    assert solution2.fluxes["NH4t"] == pytest.approx(4.66274, 0.0001)
    assert solution2.objective_value == pytest.approx(0.85511, 0.0001)


def test_user_defined_constraints_with_variable_documented(model):
    solution1 = model.optimize()
    assert solution1.objective_value == pytest.approx(0.87392, 0.0001)

    EX_glc__D_e = Variable("EX_glc__D_e")
    EX_nh4_e = Variable("EX_nh4_e")
    cc3 = Variable("difference")
    c1 = Constraint(EX_glc__D_e - EX_nh4_e - cc3, lb=0, ub=0)
    model.add_user_defined_constraints([c1])
    solution2 = model.optimize()
    assert solution2.objective_value == pytest.approx(0.87392, 0.0001)

    reaction1 = model.reactions[0]
    reaction1.knock_out()
    model.optimize()
    assert model.solver.variables.difference.primal == pytest.approx(-5.23468, 0.0001)

    reaction2 = model.reactions[1]
    reaction2.knock_out()
    model.optimize()
    assert model.solver.variables.difference.primal == pytest.approx(-5.23468, 0.0001)

    reaction3 = model.reactions[2]
    reaction3.knock_out()
    model.optimize()
    assert model.solver.variables.difference.primal == pytest.approx(-5.23468, 0.0001)

    reaction4 = model.reactions[3]
    reaction4.knock_out()
    model.optimize()
    assert model.solver.variables.difference.primal == pytest.approx(-1.86444, 0.0001)

    reaction5 = model.reactions[4]
    reaction5.knock_out()
    model.optimize()
    assert model.solver.variables.difference.primal == pytest.approx(-1.86444, 0.0001)
