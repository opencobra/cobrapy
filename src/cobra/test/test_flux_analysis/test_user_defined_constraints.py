# -*- coding: utf-8 -*-

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

from cobra.core.model import Model
from cobra.core.user_defined_constraints import (
    ConstraintComponent, UserDefinedConstraint)
from cobra.io import load_json_model, read_sbml_model, save_json_model


def ex_model(data_directory):
    model = read_sbml_model(join(data_directory, "userConstraint.xml"))
    return model


def textbook_model(data_directory):
    model = read_sbml_model(join(data_directory, "textbook.xml.gz"))
    return model


def test_user_defined_constraints(data_directory):
    cons_model = ex_model(data_directory)
    solution1 = cons_model.optimize()
    assert solution1.objective_value == pytest.approx(5.0)

    cons_comp_1 = ConstraintComponent(id='c1', variable="v1")
    cons_comp_2 = ConstraintComponent(id='c2', variable="v2")
    const_1 = UserDefinedConstraint(id="c1", lower_bound=0, upper_bound=4,
                                    const_comps=[cons_comp_1, cons_comp_2])
    cons_model.add_user_defined_constraints([const_1])
    solution2 = cons_model.optimize()
    assert solution2.objective_value == pytest.approx(2.00)


def test_user_defined_constraints_on_single_flux(data_directory):
    cons_model = ex_model(data_directory)
    solution1 = cons_model.optimize()
    assert solution1.objective_value == pytest.approx(5.0)

    cons_comp_1 = ConstraintComponent(id="cc2", variable="v2")
    const_1 = UserDefinedConstraint(id="c1", lower_bound=0, upper_bound=3,
                                    const_comps=[cons_comp_1])
    cons_model.add_user_defined_constraints([const_1])
    solution2 = cons_model.optimize()
    assert solution2.objective_value == pytest.approx(3.00)


def test_user_defined_constraints_on_single_variable():
    # an empty model
    model = Model("model_abc")
    cc1 = ConstraintComponent(id="cc1", variable="new_var")
    c1 = UserDefinedConstraint(id="c1", lower_bound=0, upper_bound=2,
                               const_comps=[cc1])
    model.add_user_defined_constraints([c1])

    model.objective = model.variables.new_var
    solution = model.optimize()
    assert solution.objective_value == pytest.approx(2.00)


def test_json_reading_writing(model, tmp_path):
    cc1 = ConstraintComponent(id="cc1", variable="FBA")
    cc2 = ConstraintComponent(id="cc2", variable="NH4t", coefficient=-1)
    cc3 = ConstraintComponent(id="cc3", variable="difference", coefficient=-1)
    c1 = UserDefinedConstraint("c1", lower_bound=0, upper_bound=0, const_comps=[cc1, cc2, cc3])

    cc4 = ConstraintComponent(id="cc4", variable="FBA")
    cc5 = ConstraintComponent(id="cc5", variable="NH4t")
    c2 = UserDefinedConstraint(id="c2", lower_bound=0, upper_bound=10, const_comps=[cc4, cc5])
    model.add_user_defined_constraints([c1, c2])

    path_to_json = join(tmp_path, "userConstraint.json")
    save_json_model(model, path_to_json)

    model = load_json_model(path_to_json)
    assert model.user_defined_const is not None
    assert len(model.user_defined_const) == 2
    const_1 = model.user_defined_const[0]
    assert const_1.id == "c1"
    assert len(const_1.constraint_comps) == 3
    assert const_1.constraint_comps[0].variable == "FBA"


def test_user_defined_constraints_documented(model):
    solution1 = model.optimize()
    assert solution1.objective_value == pytest.approx(0.87392, 0.0001)

    cc1 = ConstraintComponent(id="cc1", variable="FBA")
    cc2 = ConstraintComponent(id="cc2", variable="NH4t", coefficient=-1)
    c1 = UserDefinedConstraint(id="c1", lower_bound=0, upper_bound=0, const_comps=[cc1, cc2])
    model.add_user_defined_constraints([c1])
    solution2 = model.optimize()
    assert solution2.fluxes['FBA'] == pytest.approx(4.66274, 0.0001)
    assert solution2.fluxes['NH4t'] == pytest.approx(4.66274, 0.0001)
    assert solution2.objective_value == pytest.approx(0.85511, 0.0001)


def test_user_defined_constraints_with_variable_documented(data_directory):
    model = textbook_model(data_directory)
    solution1 = model.optimize()
    assert solution1.objective_value == pytest.approx(0.87392, 0.0001)

    cc1 = ConstraintComponent(id="cc1", variable="EX_glc__D_e")
    cc2 = ConstraintComponent(id="cc2", variable="EX_nh4_e", coefficient=-1)
    cc3 = ConstraintComponent(id="cc3", variable="difference", coefficient=-1)
    c1 = UserDefinedConstraint(id="c1", lower_bound=0, upper_bound=0, const_comps=[cc1, cc2, cc3])
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
