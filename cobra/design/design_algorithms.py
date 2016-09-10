from ..core import Model, Reaction, Metabolite
from ..manipulation.modify import canonical_form

from six import iteritems
from copy import deepcopy


def _add_decision_variable(model, reaction_id):
    """Add an integer decision variable for the given reaction."""
    reaction = model.reactions.get_by_id(reaction_id)
    # add integer variable
    var = Reaction("%s_decision_var" % reaction_id)
    var.lower_bound = 0
    var.upper_bound = 1
    var.variable_kind = "integer"
    var.decision_reaction_id = reaction_id
    model.add_reaction(var)
    # add constraints
    # v <= ub * y  -->  v - ub * y <= 0
    ub_constr = Metabolite("%s_upper_bound" % var.id)
    ub_constr._constraint_sense = "L"
    # v >= lb * y  -->  - v + lb * y <= 0
    lb_constr = Metabolite("%s_lower_bound" % var.id)
    lb_constr._constraint_sense = "L"
    reaction.add_metabolites({lb_constr: - 1,
                              ub_constr:   1})
    var.add_metabolites({lb_constr:   reaction.lower_bound,
                         ub_constr: - reaction.upper_bound})
    return var


def set_up_optknock(model, chemical_objective, knockable_reactions,
                    biomass_objective=None, n_knockouts=5,
                    n_knockouts_required=True, dual_maximum=1000, copy=True):
    """Set up the OptKnock problem described by Burgard et al., 2003:

    Burgard AP, Pharkya P, Maranas CD. Optknock: a bilevel programming
    framework for identifying gene knockout strategies for microbial strain
    optimization. Biotechnol Bioeng. 2003;84(6):647-57.
    https://doi.org/10.1002/bit.10803.

    Arguments
    ---------

    model: :class:`~cobra.core.Model`
        A COBRA model.

    chemical_objective: str
        The ID of the reaction to maximize in the outer problem.

    knockable_reactions: [str]
        A list of reaction IDs that can be knocked out.

    biomass_objective: str
        The ID of the reaction to maximize in the inner problem. By default,
        this is the existing objective function in the passed model.

    n_knockouts: int
        The number of knockouts allowable.

    n_knockouts_required: bool
        Require exactly the number of knockouts specified by n_knockouts.

    dual_maximum: float or int
        The upper bound for dual variables.

    copy: bool
        Copy the model before making any modifications.


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
            raise Exception("Could not find biomass_objective %s in model" %
                            biomass_objective)

    # dual of inner problem
    inner_dual = dual_problem(inner_problem,
                              integer_vars_to_maintain=decision_variable_ids,
                              already_irreversible=False, copy=False,
                              dual_maximum=dual_maximum)

    # add constraints and variables from inner problem to outer problem
    inner_objectives = {}
    for reaction in inner_dual.reactions:
        inner_objectives[reaction.id] = reaction.objective_coefficient
        reaction.objective_coefficient = 0
        if reaction.id in model.reactions:
            existing_reaction = model.reactions.get_by_id(reaction.id)
            for met, coeff in iteritems(reaction._metabolites):
                if met.id in model.metabolites:
                    existing_reaction.add_metabolites(
                        {model.metabolites.get_by_id(met.id): coeff})
                else:
                    existing_reaction.add_metabolites({deepcopy(met): coeff})
        else:
            model.add_reaction(reaction)

    # constraint to set outer and inner objectives equal, and set chemical
    # objective
    equal_objectives_constr = Metabolite("equal_objectives_constraint")
    equal_objectives_constr._constraint_sense = "E"
    equal_objectives_constr._bound = 0
    for reaction in model.reactions:
        if reaction.objective_coefficient != 0:
            reaction.add_metabolites({equal_objectives_constr:
                                      reaction.objective_coefficient})
        inner_objective = inner_objectives.get(reaction.id, 0)
        if inner_objective:
            reaction.add_metabolites(
                {equal_objectives_constr: - inner_objective})
        # set chemical objective
        reaction.objective_coefficient = 1 \
            if reaction.id == chemical_objective else 0

    # add the n_knockouts constraint
    n_knockouts_constr = Metabolite("n_knockouts_constraint")
    n_knockouts_constr._constraint_sense = "E" if n_knockouts_required else "G"
    n_knockouts_constr._bound = len(decision_variable_ids) - n_knockouts
    for r_id in decision_variable_ids:
        reaction = model.reactions.get_by_id(r_id)
        reaction.add_metabolites({n_knockouts_constr: 1})

    return model


def run_optknock(optknock_problem, solver=None, tolerance_integer=1e-9,
                 **kwargs):
    """Run the OptKnock problem created with set_up_optknock.

    Arguments
    ---------

    optknock_problem: :class:`~cobra.core.Model`
        The problem generated by set_up_optknock.

    solver: str
        The name of the preferred solver.

    tolerance_integer: float
        The integer tolerance for the MILP.

    **kwargs
        Keyword arguments are passed to Model.optimize().


    Zachary King 2015

    """
    solution = optknock_problem.optimize(solver=solver,
                                         tolerance_integer=tolerance_integer,
                                         **kwargs)
    solution.knockouts = []
    for reaction in optknock_problem.reactions:
        if solution.x_dict.get(reaction.id, None) == 0:
            r_id = getattr(reaction, "decision_reaction_id", None)
            if r_id is not None:
                solution.knockouts.append(r_id)
    return solution


# This function will generalize the set_up_optknock code to other MILPs:
# def dual_embed(outer_model, inner_model, ..., objective_sense="maximize",
#                integer_vars_to_maintain=[], already_irreversible=False,
#                copy=True, dual_maximum=1000):
#     """Embed the dual of the inner model within the outer model"""


def dual_problem(model, objective_sense="maximize",
                 integer_vars_to_maintain=[],
                 already_irreversible=False, copy=True, dual_maximum=1000):
    """Return a new model representing the dual of the model.

    Make the problem irreversible, then take the dual. Convert the problem:

    .. code-block:: none

        Maximize (c^T)x subject to Ax <= b, x >= 0

    which is something like this in COBRApy:

    .. code-block:: none

        Maximize sum(objective_coefficient_j * reaction_j for all j)
            s.t.
            sum(coefficient_i_j * reaction_j for all j) <= metabolite_bound_i
            reaction_j <= upper_bound_j
            reaction_j >= 0

    to the problem:

    .. code-block:: none

        Minimize (b^T)w subject to (A^T)w >= c, w >= 0

    which is something like this in COBRApy (S matrix is m x n):

    .. code-block:: none

        Minimize sum( metabolite_bound_i * dual_i   for all i ) +
                 sum( upper_bound_j *      dual_m+j for all j ) +
            s.t.
             sum( coefficient_i_j * dual_i for all i ) +
             sum( dual_2m+j' for all j' ) >= objective_coefficient_j
            dual_k >= 0


    Arguments
    ---------

    model : :class:`~cobra.core.Model`
        The COBRA model.

    objective_sense: str
        The objective sense of the starting problem, either 'maximize' or
        'minimize'. A minimization problems will be converted to a maximization
        before taking the dual. This function always returns a minimization
        problem.

    iteger_vars_to_maintain: [str]
        A list of IDs for Boolean integer variables to be maintained in the
        dual problem. See 'Maintaining integer variables' below for more
        details.

    already_irreversible: bool
        If True, then do not convert the model to irreversible.

    copy: bool
        If True, then make a copy of the model before modifying it. This is not
        necessary if already_irreversible is True.

    dual_maximum: float or int
        The upper bound for dual variables.


    **Maintaining integer variables**

    The argument ``integer_vars_to_maintain`` can be used to specify certin
    Boolean integer variables that will be maintained in the dual problem. This
    makes it possible to join outer and inner problems in a bi-level MILP. The
    method for maintaining integer variables is described by Tepper and Shlomi,
    2010:

    Tepper N, Shlomi T. Predicting metabolic engineering knockout strategies
    for chemical production: accounting for competing pathways. Bioinformatics.
    2010;26(4):536-43. https://doi.org/10.1093/bioinformatics/btp704.

    In COBRApy, this roughly translates to transforming (decision variables p,
    integer constraints o):

    .. code-block:: none

        Maximize (c^T)x subject to (A_x)x + (A_y)y <= b, x >= 0

        (1) Maximize sum(objective_coefficient_j * reaction_j for all j)
                s.t.
        (2)     sum(coeff_i_j * reaction_j for all j) +
                sum(decision_coeff_i_j * decision_var_j for all j)
                <= metabolite_bound_i
        (3)     reaction_j <= upper_bound_j
        (4)     reaction_j >= 0

    to the problem:

    .. code-block:: none

        Minimize (b - (A_y)y)^T w subject to (A_x^T)w >= c, w >= 0

    which linearizes to (with auxiliary variables z):

    .. code-block:: none

        Minimize (b^T)w - { ((A_y)y)^T w with yw --> z }
        subject to (A_x^T)w >= c, linearization constraints, w >= 0
          Linearization constraints: z <= w_max * y, z <= w,
                                     z >= w - w_max * (1 - y), z >= 0

        (5) Minimize sum( metabolite_bound_i *  dual_i            for all i ) +
                      sum( upper_bound_j *      dual_m+j          for all j ) +
                    - sum( decision_coeff_i_j * auxiliary_var_i_j
                          for all combinations i, j )
                s.t.
        (6)   - sum( coefficient_i_j * dual_i for all i ) - dual_m+j
              <= - objective_coefficient_j
        (7)     auxiliary_var_i_j - dual_maximum * decision_var_j          <= 0
        (8)     auxiliary_var_i_j - dual_i                                 <= 0
        (9)   - auxiliary_var_i_j + dual_i + dual_maximum * decision_var_j
              <= dual_maximum
       (10)     dual_maximum >= dual_i            >= 0
       (11)     dual_maximum >= dual_m+j          >= 0
       (12)     dual_maximum >= auxiliary_var_i_j >= 0
       (13)                1 >= decision_var_j    >= 0


    Zachary King 2015

    """

    # convert to canonical form and copy
    model = canonical_form(model, objective_sense=objective_sense,
                           already_irreversible=already_irreversible,
                           copy=copy)

    # new model for the dual
    dual = Model("%s_dual" % model.id)

    # keep track of dual_i
    dual_var_for_met = {}

    # add dual variables for constraints. (2) --> dual_i
    for metabolite in model.metabolites:
        # add constraints based on metabolite constraint sense
        if metabolite._constraint_sense != "L":
            raise Exception("Not a less than or equal constraint: %s"
                            % metabolite.id)

        var = Reaction("%s__dual" % metabolite.id)
        # Without auxiliary variables, the objective coefficient would include
        # integer variables when present. However, we will separate out the
        # integer parts into objective coefficients for auxiliary variables.
        var.objective_coefficient = metabolite._bound  # (5)
        # [dual_vars] >= 0
        var.lower_bound = 0
        var.upper_bound = dual_maximum
        dual.add_reaction(var)
        # remember
        dual_var_for_met[metabolite.id] = var

    # keep track of decision variables (integer_vars_to_maintain) as tuples:
    # (reaction in dual problem, reaction in original problem)
    integer_vars_added = []

    # add constraints and upper bound variables
    for reaction in model.reactions:
        # integer vars to maintain
        if reaction.id in integer_vars_to_maintain:
            # keep these integer variables in the dual, with new transformed
            # constraints
            if (reaction.lower_bound not in [0, 1] or
                    reaction.upper_bound not in [0, 1] or
                    reaction.variable_kind != "integer"):
                raise Exception("Reaction %s from integer_vars_to_maintain is "
                                "not a Boolean integer variable" % reaction.id)
            integer_var = Reaction(reaction.id)
            integer_var.upper_bound = reaction.upper_bound
            integer_var.lower_bound = reaction.lower_bound
            integer_var.variable_kind = reaction.variable_kind
            integer_var.objective_coefficient = 0
            # constraints
            dual.add_reaction(integer_var)
            integer_vars_added.append((integer_var, reaction))

        # other vars
        else:
            # other variables become constraints, (1) to (6)
            constr = Metabolite("%s__dual_constrained_by_c" %
                                reaction.id)  # (6)
            constr._constraint_sense = "L"
            constr._bound = - reaction.objective_coefficient
            for met, coeff in iteritems(reaction._metabolites):
                dual_var = dual_var_for_met[met.id]
                dual_var.add_metabolites({constr: - coeff})

            # upper bound constraints -> variables (3) to (5) and (6)
            var_bound = Reaction("%s__dual_for_upper_bound_constraint" %
                                 reaction.id)  # dual_m+j
            var_bound.objective_coefficient = reaction.upper_bound  # (5)
            # [dual_vars] >= 0
            var_bound.lower_bound = 0
            var_bound.upper_bound = dual_maximum
            # add bound dual variables to dual constraints
            var_bound.add_metabolites({constr: -1})  # (6)
            dual.add_reaction(var_bound)

    # add auxiliary variables
    for integer_var, original_reaction in integer_vars_added:
        for metabolite, coeff in iteritems(original_reaction._metabolites):
            dual_var = dual_var_for_met[metabolite.id]
            # create an auxiliary variable
            aux_var = Reaction("%s__auxiliary__%s" % (integer_var.id,
                                                      dual_var.id))
            aux_var.lower_bound = 0
            aux_var.upper_bound = dual_maximum
            aux_var.objective_coefficient = - coeff
            dual.add_reaction(aux_var)

            # add linearization constraints
            # (7)     auxiliary_var_i_j - dual_maximum * decision_var_j    <= 0
            le_decision_constr = Metabolite("%s__le_decision" % aux_var.id)
            le_decision_constr._constraint_sense = "L"
            le_decision_constr._bound = 0
            aux_var.add_metabolites({le_decision_constr: 1})
            integer_var.add_metabolites({le_decision_constr: - dual_maximum})

            # (8)     auxiliary_var_i_j - dual_i                           <= 0
            le_dual_constr = Metabolite("%s__le_dual" % aux_var.id)
            le_dual_constr._constraint_sense = "L"
            le_dual_constr._bound = 0
            aux_var.add_metabolites({le_dual_constr: 1})
            dual_var.add_metabolites({le_dual_constr: -1})

            # (9)   - auxiliary_var_i_j + dual_i +
            #         dual_maximum * decision_var_j <= dual_maximum
            g_constr = Metabolite("%s__g_dual" % aux_var.id)
            g_constr._constraint_sense = "L"
            g_constr._bound = dual_maximum
            aux_var.add_metabolites({g_constr: -1})
            dual_var.add_metabolites({g_constr: 1})
            integer_var.add_metabolites({g_constr: dual_maximum})

    return dual
