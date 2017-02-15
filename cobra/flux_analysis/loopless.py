# -*- coding: utf-8 -*-

from __future__ import absolute_import

from cobra.core import Metabolite, Reaction
from cobra.manipulation.modify import convert_to_irreversible
import numpy as np
from six import iteritems
from cobra.util import nullspace, add_to_solver
from sympy.core.singleton import S


def add_loopless(model, zero_cutoff=1e-12):
    """Add variables and constraints to make a model thermodynamically
    feasible. This removes flux loops. The used formulation is described
    in [1]_.

    Parameters
    ----------
    model : a cobra model
        The model to which to add the constraints.
    zero_cutoff : positive float, optional
        Cutoff used for null space. Coefficients with an absolute value smaller
        than `zero_cutoff` are considered to be zero.

    Returns
    -------
    Nothing

    References
    ----------
    .. [1] Elimination of thermodynamically infeasible loops in steady-state
       metabolic models.
       Schellenberger J, Lewis NE, Palsson BO.
       Biophys J. 2011 Feb 2;100(3):544-53.
       doi: 10.1016/j.bpj.2010.12.3707.
       Erratum in: Biophys J. 2011 Mar 2;100(5):1381.
    """
    internal = [i for i, r in enumerate(model.reactions) if not r.boundary]
    Sint = model.S[:, np.array(internal)]
    Nint = nullspace(Sint).T
    max_bound = max(max(abs(b) for b in r.bounds) for r in model.reactions)
    prob = model.solver.interface

    # Add indicator variables and new constraints
    to_add = []
    for i in internal:
        rxn = model.reactions[i]
        # indicator variable a_i
        indicator = prob.Variable("indicator_" + rxn.id, type="binary")
        # -M*(1 - a_i) <= v_i <= M*a_i
        on_off_constraint = prob.Constraint(
            rxn.flux_expression - max_bound * indicator,
            lb=-max_bound, ub=0, name="on_off_" + rxn.id)
        # -(max_bound + 1) * a_i + 1 <= G_i <= -(max_bound + 1) * a_i + 1000
        delta_g = prob.Variable("delta_g_" + rxn.id)
        delta_g_range = prob.Constraint(
            delta_g + (max_bound + 1) * indicator,
            lb=1, ub=max_bound, name="delta_g_range_" + rxn.id)
        to_add.extend([indicator, on_off_constraint, delta_g, delta_g_range])

    add_to_solver(model, to_add)

    # Add nullspace constraints for G_i
    for i, row in enumerate(Nint):
        name = "nullspace_constraint_" + str(i)
        nullspace_constraint = prob.Constraint(S.Zero, lb=0, ub=0, name=name)
        add_to_solver(model, [nullspace_constraint])
        coefs = {model.solver.variables[
            "delta_g_" + model.reactions[ridx].id]: row[i]
            for i, ridx in enumerate(internal) if abs(row[i]) > zero_cutoff}
        model.solver.constraints[name].set_linear_coefficients(coefs)


def construct_loopless_model(cobra_model):
    """construct a loopless model

    This adds MILP constraints to prevent flux from proceeding in a loop, as
    done in http://dx.doi.org/10.1016/j.bpj.2010.12.3707
    Please see the documentation for an explanation of the algorithm.

    This must be solved with an MILP capable solver.

    """
    # copy the model and make it irreversible
    model = cobra_model.copy()
    convert_to_irreversible(model)
    max_ub = max(model.reactions.list_attr("upper_bound"))
    # a dict for storing S^T
    thermo_stoic = {"thermo_var_" + metabolite.id: {}
                    for metabolite in model.metabolites}
    # Slice operator is so that we don't get newly added metabolites
    original_metabolites = model.metabolites[:]
    for reaction in model.reactions[:]:
        # Boundary reactions are not subjected to these constraints
        if len(reaction._metabolites) == 1:
            continue
        # populate the S^T dict
        bound_id = "thermo_bound_" + reaction.id
        for met, stoic in iteritems(reaction._metabolites):
            thermo_stoic["thermo_var_" + met.id][bound_id] = stoic
        # I * 1000 > v --> I * 1000 - v > 0
        reaction_ind = Reaction(reaction.id + "_indicator")
        reaction_ind.variable_kind = "integer"
        reaction_ind.upper_bound = 1
        reaction_ub = Metabolite(reaction.id + "_ind_ub")
        reaction_ub._constraint_sense = "G"
        reaction.add_metabolites({reaction_ub: -1})
        reaction_ind.add_metabolites({reaction_ub: max_ub})
        # This adds a compensating term for 0 flux reactions, so we get
        # S^T x - (1 - I) * 1001 < -1 which becomes
        # S^T x < 1000 for 0 flux reactions and
        # S^T x < -1 for reactions with nonzero flux.
        reaction_bound = Metabolite(bound_id)
        reaction_bound._constraint_sense = "L"
        reaction_bound._bound = max_ub
        reaction_ind.add_metabolites({reaction_bound: max_ub + 1})
        model.add_reaction(reaction_ind)
    for metabolite in original_metabolites:
        metabolite_var = Reaction("thermo_var_" + metabolite.id)
        metabolite_var.lower_bound = -max_ub
        model.add_reaction(metabolite_var)
        metabolite_var.add_metabolites(
            {model.metabolites.get_by_id(k): v
             for k, v in iteritems(thermo_stoic[metabolite_var.id])})
    return model
