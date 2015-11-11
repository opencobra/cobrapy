from ..core import Model, Reaction, Metabolite
from ..manipulation.modify import convert_to_irreversible

from six import iteritems
from collections import defaultdict


def dual_problem(model, already_irreversible=False, copy=True,
                 variable_max=1000):
    """Return a new model representing the dual of the model.

    Make the problem irreversible, then take the dual. Convert the problem:

        Maximize (c^T)x subject to Ax <= b, x >= 0

    which is something like this in COBRApy:

        Maximize [objective_coefficients] * [reactions]
            subject to
                (assuming met._sense is 'E')
                [coefficients by reaction] * [reactions] <= met._bound
                - [coefficients by reaction] * [reactions] <= - met._bound
                [reactions] <= [UB]
                - [reactions] <= - [LB]
                [reactions] >= 0

    to the problem:

        Minimize (b^T)y subject to (A^T)y >= c, y >= 0

    which is something like this in COBRApy:

        Minimize [met._bound's, - met._bound's, UB's, - LB's] * [dual_variables]
            subject to
                [coefficients by metabolite, - coefficients by metabolite, diag, - diag] * [dual_variables] >= [objective_coefficients]
                [dual_variables] >= 0


    model : :class:`~cobra.core.Model` object.

    already_irreversible: Boolean. If True, then do not convert the model to
    irreversible.

    copy: Boolean. If True, then make a copy of the model before modifying
    it. This is not necessary if already_irreversible is True.

    variable_max: Number. The upper bound for dual variables.

    """

    # make irreversible and copy
    irrev_model = model
    if not already_irreversible:
        if copy:
            irrev_model = irrev_model.copy()
        convert_to_irreversible(irrev_model)

    # new model
    dual = Model("%s_dual" % irrev_model.id)

    # keep track of variables
    dual_vars_for_met = defaultdict(list)

    # add variables and objective coefficients
    for metabolite in irrev_model.metabolites:
        # S v = 0 constraints
        if metabolite._constraint_sense != "E":
            raise NotImplementedError("Only metabolite._constraint_sense='E' is supported")
        b_i = metabolite._bound

        # stoichiometric matrix dual variables
        for sense, factor in ("LE", 1), ("GE", -1):
            var = Reaction("%s_dual_%s" % (metabolite.id, sense))
            var.objective_coefficient = factor * b_i
            # [dual_variables] >= 0
            var.lower_bound = 0
            var.upper_bound = variable_max
            dual.add_reaction(var)
            # remember
            dual_vars_for_met[metabolite.id].append(var)

    # add constraints and upper & lower bound variables
    for reaction in irrev_model.reactions:
        # objective function -> constraints
        constr = Metabolite("%s_dual_constraint" % reaction.id)
        constr._constraint_sense = "G"
        constr._bound = reaction.objective_coefficient
        for met, coeff in iteritems(reaction.metabolites):
            var_le, var_ge = dual_vars_for_met[met.id]
            var_le.add_metabolites({ constr: coeff })
            var_ge.add_metabolites({ constr: - coeff })

        for bound_name, factor in ("upper_bound", 1), ("lower_bound", -1):
            # upper/lower bound constraints -> variables
            var_bound = Reaction("%s_dual_%s" % (reaction.id, bound_name))
            var_bound.objective_coefficient = factor * getattr(reaction, bound_name)
            # [dual_variables] >= 0
            var_bound.lower_bound = 0
            var_bound.upper_bound = variable_max
            # add bound dual variables to dual constraints
            var_bound.add_metabolites({ constr: factor })
            dual.add_reaction(var_bound)

    return dual
