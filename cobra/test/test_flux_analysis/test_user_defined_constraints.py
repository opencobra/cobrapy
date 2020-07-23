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
    UserDefinedConstraintComponents, UserDefinedConstraints)
from cobra.io import load_json_model, read_sbml_model, save_json_model


@pytest.fixture(scope="session")
def ex_model(data_directory):
    model = read_sbml_model(join(data_directory, "userConstraint.xml"))
    return model


@pytest.fixture(scope="function")
def cons_model(ex_model):
    return ex_model.copy()


def test_user_defined_constraints(cons_model):
    solution1 = cons_model.optimize()
    assert solution1.objective_value == 5.0

    cons_comp_1 = UserDefinedConstraintComponents("cc1", "comp1",
                                                  "v1", 1, "linear")
    cons_comp_2 = UserDefinedConstraintComponents("cc2", "comp2",
                                                  "v2", 1, "linear")
    const_1 = UserDefinedConstraints("c1", "constraint1",
                                     0, 4, [cons_comp_1, cons_comp_2])
    cons_model.add_user_defined_constraints([const_1])
    solution2 = cons_model.optimize()
    assert solution2.objective_value == 2.00


def test_user_defined_constraints_on_single_flux(cons_model):
    solution1 = cons_model.optimize()
    assert solution1.objective_value == 5.0

    cons_comp_1 = UserDefinedConstraintComponents("cc2", "comp2",
                                                  "v2", 1, "linear")
    const_1 = UserDefinedConstraints("c1", "constraint1", 0, 3, [cons_comp_1])
    cons_model.add_user_defined_constraints([const_1])
    solution2 = cons_model.optimize()
    assert solution2.objective_value == 3.00


def test_user_defined_constraints_on_single_variable():
    # an empty model
    model = Model("model_abc")
    cc1 = UserDefinedConstraintComponents("cc1", "comp1",
                                          "new_var", 1, "linear")
    c1 = UserDefinedConstraints("c1", "constraint1", 0, 2, [cc1])
    model.add_user_defined_constraints([c1])

    model.objective = model.variables.new_var
    solution = model.optimize()
    assert solution.objective_value == 2.00


def test_json_reading_writing(small_model, tmp_path):
    cc1 = UserDefinedConstraintComponents("cc1", "comp1", "FBA", 1, "linear")
    cc2 = UserDefinedConstraintComponents("cc2", "comp2", "NH4t", -1, "linear")
    cc3 = UserDefinedConstraintComponents("cc3", "comp3",
                                          "difference", -1, "linear")
    c1 = UserDefinedConstraints("c1", "constraint1", 0, 0, [cc1, cc2, cc3])

    cc4 = UserDefinedConstraintComponents("cc4", "comp4", "FBA", 1, "linear")
    cc5 = UserDefinedConstraintComponents("cc5", "comp5", "NH4t", 1, "linear")
    c2 = UserDefinedConstraints("c2", "constraint2", 0, 10, [cc4, cc5])
    small_model.add_user_defined_constraints([c1, c2])

    path_to_json = join(tmp_path, "userConstraint.json")
    save_json_model(small_model, path_to_json)

    model = load_json_model(path_to_json)
    assert model.user_defined_const is not None
    assert len(model.user_defined_const) == 2
    const_1 = model.user_defined_const[0]
    assert const_1.id == "c1"
    assert len(const_1.constraint_comps) == 3
    assert const_1.constraint_comps[0].ref_var == "FBA"


def test_user_defined_constraints_documented(small_model):
    solution1 = small_model.optimize()
    assert solution1.objective_value == 0.8739215069684307

    cc1 = UserDefinedConstraintComponents("cc1", "comp1", "FBA", 1, "linear")
    cc2 = UserDefinedConstraintComponents("cc2", "comp2", "NH4t", -1, "linear")
    c1 = UserDefinedConstraints("c1", "constraint1", 0, 0, [cc1, cc2])
    small_model.add_user_defined_constraints([c1])
    solution2 = small_model.optimize()
    assert solution2.fluxes['FBA'] == 4.662749047738146
    assert solution2.fluxes['NH4t'] == 4.662749047738146
    assert solution2.objective_value == 0.8551109609261563


def test_user_defined_constraints_with_variable_documented(small_model):
    solution1 = small_model.optimize()
    assert solution1.objective_value == 0.8739215069684307

    cc1 = UserDefinedConstraintComponents("cc1", "comp1", "FBA", 1, "linear")
    cc2 = UserDefinedConstraintComponents("cc2", "comp2", "NH4t", -1, "linear")
    cc3 = UserDefinedConstraintComponents("cc3", "comp3",
                                          "difference", -1, "linear")
    c1 = UserDefinedConstraints("c1", "constraint1", 0, 0, [cc1, cc2, cc3])
    small_model.add_user_defined_constraints([c1])
    # solution2 = small_model.optimize()
    # assert solution2.objective_value == 0.8739215069684302

    # reaction1 = small_model.reactions[0]
    # reaction1.knock_out()
    # small_model.optimize()
    # assert small_model.solver.variables.difference.primal == -5.234680806802543
    #
    # reaction2 = small_model.reactions[1]
    # reaction2.knock_out()
    # small_model.optimize()
    # assert small_model.solver.variables.difference.primal == -5.2346808068025386
    #
    # reaction3 = small_model.reactions[2]
    # reaction3.knock_out()
    # small_model.optimize()
    # assert small_model.solver.variables.difference.primal == -5.234680806802525
    #
    # reaction4 = small_model.reactions[3]
    # reaction4.knock_out()
    # small_model.optimize()
    # assert small_model.solver.variables.difference.primal == -1.8644444444444337
    #
    # reaction5 = small_model.reactions[4]
    # reaction5.knock_out()
    # small_model.optimize()
    # assert small_model.solver.variables.difference.primal == -1.8644444444444466
