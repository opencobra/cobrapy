from ..core import Model, Reaction, Metabolite
from ..manipulation.modify import convert_to_irreversible

from six import iteritems
from collections import defaultdict


def _add_decision_variable(model, reaction_id):
    """Add an integer decision variable for the given reaction."""
    reaction = model.reactions.get_by_id(reaction_id)
    # add integer variable
    var = Reaction("%s_decision_var" % reaction_id)
    var.lower_bound = 0
    var.upper_bound = 1
    var.variable_kind = "integer"
    model.add_reaction(var)
    # add constraints
    # v <= ub * y  -->  v - ub * y <= 0
    ub_constr = Metabolite("%s_upper_bound" % var.id)
    ub_constr._constraint_sense = "L"
    # v >= lb * y  -->  v - lb * y >= 0
    lb_constr = Metabolite("%s_lower_bound" % var.id)
    lb_constr._constraint_sense = "G"
    reaction.add_metabolites({ lb_constr: 1, ub_constr: 1 })
    var.add_metabolites({ lb_constr: - reaction.lower_bound,
                          ub_constr: - reaction.upper_bound })
    return var


def optknock(model, chemical_objective, knockable_reactions,
             biomass_objective=None, n_knockouts=5, n_knockouts_required=True,
             dual_maximum=1000, copy=True):
    """Run the OptKnock algorithm described by Burgard et al., 2003:

        Burgard AP, Pharkya P, Maranas CD. Optknock: a bilevel programming
        framework for identifying gene knockout strategies for microbial strain
        optimization.  Biotechnol Bioeng. 2003;84(6):647–57. doi:10.1002/bit.10803.


    model : :class:`~cobra.core.Model` object.

    chemical_objective: str. The ID of the reaction to Maximize in the outer
    problem.

    knockable_reactions: [str]. A list of reaction IDs that can be knocked out.

    biomass_objective: str. The ID of the reaction to Maximize in the inner
    problem. By default, this is the existing objective function in the passed
    model.

    n_knockouts: int. The number of knockouts allowable.

    n_knockouts_required: bool. Require exactly the number of knockouts
    specified by n_knockouts.

    dual_maximum: float or int. The upper bound for dual variables.

    copy: bool. Copy the model before making any modifications.


    Zachary King 2015

    """

    if copy:
        model = model.copy()

    # add the integer decision variables
    decision_variable_ids = [_add_decision_variable(model, r_id).id
                             for r_id in knockable_reactions]

    # inner problem
    inner_problem = model.copy()
    if biomass_objective:
        found = False
        for reaction in inner_problem.reactions:
            obj = reaction.id == biomass_objective
            reaction.objective_coefficient = 1 if obj else 0
            if obj:
                found = True
        if not found:
            raise Exception("Could not find biomass_objective %s in model" % biomass_objective)

    # dual of inner problem
    inner_dual = dual_problem(inner_problem, integer_vars_to_maintain=decision_variable_ids,
                              already_irreversible=False, copy=False,
                              dual_maximum=dual_maximum)

    # add constraints and variables from inner problem to outer problem
    model.add_reactions(inner_dual.reactions)

    # constraint to set outer and inner objectives to be equal
    equal_objectives_constr = Metabolite("equal_objectives_constr")
    equal_objectives_constr._constraint_sense = "E"

    raise NotImplementedError


def dual_problem(model, integer_vars_to_maintain=[], already_irreversible=False,
                 copy=True, dual_maximum=1000):
    """Return a new model representing the dual of the model.

    Make the problem irreversible, then take the dual. Convert the problem:

        Maximize (c^T)x subject to Ax <= b, x >= 0

    which is something like this in COBRApy:

        Maximize [objective_coefficients] * [reactions]
            subject to
                (assuming met._sense is "E")
                [coefficients by reaction] * [reactions] <= met._bound
                - [coefficients by reaction] * [reactions] <= - met._bound
                [reactions] <= [UB]
                - [reactions] <= - [LB]
                [reactions] >= 0

    to the problem:

        Minimize (b^T)y subject to (A^T)y >= c, y >= 0

    which is something like this in COBRApy:

        Minimize [met._bound's, - met._bound's, UB's, - LB's] * [dual_vars]
            subject to
                [coefficients by metabolite, - coefficients by metabolite, diag, - diag] * [dual_vars] >= [objective_coefficients]
                [dual_vars] >= 0


    Arguments
    ---------

    model : :class:`~cobra.core.Model` object.

    iteger_vars_to_maintain: [str]. A list of IDs for Boolean integer variables
    to be maintained in the dual problem. See 'Maintaining integer variables'
    below for more details

    already_irreversible: Boolean. If True, then do not convert the model to
    irreversible.

    copy: bool. If True, then make a copy of the model before modifying
    it. This is not necessary if already_irreversible is True.

    dual_maximum: float or int. The upper bound for dual variables.


    Maintaining integer variables
    -----------------------------

    The argument integer_vars_to_maintain can be used to specify certin Boolean
    integer variables that will be maintained in the dual problem. This makes it
    possible to join outer and inner problems in a bi-level MILP. The method for
    maintaining integer variables is described by Tepper and Shlomi, 2010:

        Tepper N, Shlomi T. Predicting metabolic engineering knockout strategies
        for chemical production: accounting for competing pathways. Bioinformatics.
        2010;26(4):536–43. doi:10.1093/bioinformatics/btp704.

    In COBRApy, this roughly translates to transforming:

        (1) Maximize [objective_coefficients] * [reactions]
              subject to
        (2)       [coefficients by reaction] * [reactions] <=   met._bound ("L" or "E")
        (3)     - [coefficients by reaction] * [reactions] <= - met._bound ("G" or "E")
        (4)       [reactions] <=   [UB]
        (5)     - [reactions] <= - [LB]   (these are still necessary because they are attributes of Reaction objects)
        (6)       [reactions] <=   [decision_vars] .* [UB]   (optional)
        (7)     - [reactions] <= - [decision_vars] .* [LB]
        (8)       [reactions decision_vars] >= 0

    to the problem:

        (9) Minimize [met._bound's, - met._bound's, UB's, - LB's] * [dual_vars] + [aux_to_decision] * [auxiliary_vars]
              subject to
       (10)       [coefficients by metabolite, - coefficients by metabolite, diag, - diag] * [dual_vars] >= [objective_coefficients]
       (11)       [auxiliary_vars] - dual_maximum * corresponding([decision_vars]) <= 0
       (12)     - [auxiliary_vars] + corresponding([dual_vars]) >= 0
       (13)       [auxiliary_vars] - (corresponding([dual_vars]) - dual_maximum * (1 - correspond([decision_vars]))) >= 0
       (14)       [dual_vars, auxiliary_vars] >= 0

        aux_to_decision gives the coefficients of the decision_variables that
        correspond to each auxiliary variable. Auxiliary variables exist for
        each pair of decision variable and dual for the decision variable
        constraint.



    Zachary King 2015

    """

    # make irreversible and copy
    irrev_model = model
    if not already_irreversible:
        if copy:
            irrev_model = irrev_model.copy()
        convert_to_irreversible(irrev_model)

    # new model for the dual
    dual = Model("%s_dual" % irrev_model.id)

    # keep track of dual variables for the new
    dual_vars_for_met = defaultdict(lambda: {"LE": None, "GE": None})

    # add variables and objective coefficients
    for metabolite in irrev_model.metabolites:
        # add constraints based on metabolite constraint sense
        sense_tuples = []
        if metabolite._constraint_sense in ("L", "E"):
            sense_tuples.append(("LE", 1))
        if metabolite._constraint_sense in ("G", "E"):
            sense_tuples.append(("GE", -1))
        if metabolite._constraint_sense not in ("G", "E", "L"):
            raise Exception("Bad metabolite._constraint_sense = %s" % metabolite._constraint_sense)

        # constraints (2-7) to objective (9)
        for sense, factor in sense_tuples:
            var = Reaction("%s_dual_%s" % (metabolite.id, sense))
            # Without auxiliary variables, the objective coefficient would
            # include integer variables when present (for (6) and (7)). However,
            # we will separate out the integer parts into the auxiliary variable
            # objective coefficients.
            var.objective_coefficient = factor * metabolite._bound
            # [dual_vars] >= 0
            var.lower_bound = 0
            var.upper_bound = dual_maximum
            dual.add_reaction(var)
            # remember
            dual_vars_for_met[metabolite.id][sense] = var

    # add constraints and upper & lower bound variables
    for reaction in irrev_model.reactions:
        if reaction.id in integer_vars_to_maintain:
            # keep these integer variables in the dual, with new transformed
            # constraints
            if (reaction.lower_bound not in [0, 1] or
                    reaction.upper_bound not in [0, 1] or
                    reaction.variable_kind != "integer"):
                raise Exception("Reaction %s from integer_vars_to_maintain is not a Boolean integer variable" % reaction.id)
            decision_var = Reaction(reaction.id)
            decision_var.upper_bound = reaction.upper_bound
            decision_var.lower_bound = reaction.lower_bound
            decision_var.variable_kind = reaction.variable_kind
            decision_var.objective_coefficient = 0
            # constraints
            dual.add_reaction(decision_var)

            # add auxiliary variable z_k for each pair y_i, w_j
            # y_i : reaction (and decision_var)
            # w_j : dual_vars_for_met[met.id]["LE"] and dual_vars_for_met[met.id]["GE"]
            # A^y_i,j : coeff

            # get the associated duals
            for met, coeff in iteritems(reaction.metabolites):
                duals = dual_vars_for_met[met.id]
                for var, factor in (duals["LE"], 1), (duals["GE"], -1):
                    if var is None:
                        continue
                    # make the auxiliary variable
                    aux_var = Reaction("%s_auxiliary" % decision_var.id)
                    aux_var.lower_bound = 0
                    aux_var.upper_bound = dual_maximum
                    aux_var.variable_kind = "continuous"
                    aux_var.objective_coefficient = factor * coeff

                    # add auxiliary constraints (11-13)
                    # (11)       [auxiliary_vars] - dual_maximum * corresponding([decision_vars]) <= 0
                    # (12)     - [auxiliary_vars] + corresponding([dual_vars]) >= 0
                    # (13)       [auxiliary_vars] - (corresponding([dual_vars]) - dual_maximum * (1 - correspond([decision_vars]))) >= 0

                    raise NotImplementedError

        else:
            # other variables become constraints, (1) to (10)
            constr = Metabolite("%s_dual_constraint" % reaction.id)
            constr._constraint_sense = "G"
            constr._bound = reaction.objective_coefficient
            for met, coeff in iteritems(reaction.metabolites):
                duals = dual_vars_for_met[met.id]
                var_le, var_ge = duals["LE"], duals["GE"]
                if var_le is not None:
                    var_le.add_metabolites({ constr: coeff })
                if var_ge is not None:
                    var_ge.add_metabolites({ constr: - coeff })

            for bound_name, factor in ("upper_bound", 1), ("lower_bound", -1):
                # upper/lower bound constraints -> variables
                var_bound = Reaction("%s_dual_%s" % (reaction.id, bound_name))
                var_bound.objective_coefficient = factor * getattr(reaction, bound_name)
                # [dual_vars] >= 0
                var_bound.lower_bound = 0
                var_bound.upper_bound = dual_maximum
                # add bound dual variables to dual constraints
                var_bound.add_metabolites({ constr: factor })
                dual.add_reaction(var_bound)

    return dual
