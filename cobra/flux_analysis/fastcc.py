# -*- coding: utf-8 -*-

"""Provide an implementation of FASTCC."""

from __future__ import absolute_import

import logging

from optlang.symbolics import Zero

from cobra.flux_analysis.helpers import normalize_cutoff


LOGGER = logging.getLogger(__name__)


def fastcc(model, flux_threshold=1.0, zero_cutoff=None):
    r"""
    Check consistency of a metabolic network using FASTCC [1]_.

    FASTCC (Fast Consistency Check) is an algorithm for rapid and
    efficient consistency check in metabolic networks. FASTCC is
    a pure LP implementation and is low on computation resource
    demand. FASTCC also circumvents the problem associated with
    reversible reactions for the purpose. Given a global model,
    it will generate a consistent global model i.e., remove
    blocked reactions. For more details on FASTCC, please
    check [1]_.

    Parameters
    ----------
    model: cobra.Model
        The constraint-based model to operate on.
    flux_threshold: float, optional (default 1.0)
        The flux threshold to consider.
    zero_cutoff: float, optional
        The cutoff to consider for zero flux (default model.tolerance).

    Returns
    -------
    cobra.Model
        The consistent constraint-based model.

    Notes
    -----
    The LP used for FASTCC is like so:
    maximize: \sum_{i \in J} z_i
    s.t.    : z_i \in [0, \varepsilon] \forall i \in J, z_i \in \mathbb{R}_+
              v_i \ge z_i \forall i \in J
              Sv = 0 v \in B

    References
    ----------
    .. [1] Vlassis N, Pacheco MP, Sauter T (2014)
           Fast Reconstruction of Compact Context-Specific Metabolic Network
           Models.
           PLoS Comput Biol 10(1): e1003424. doi:10.1371/journal.pcbi.1003424

    """
    zero_cutoff = normalize_cutoff(model, zero_cutoff)

    with model:
        obj_vars = []
        vars_and_cons = []
        prob = model.problem

        for rxn in model.reactions:
            var = prob.Variable("auxiliary_{}".format(rxn.id),
                                lb=0.0, ub=flux_threshold)
            const = prob.Constraint(rxn.forward_variable +
                                    rxn.reverse_variable -
                                    var, name="constraint_{}".format(rxn.id),
                                    lb=0.0)
            vars_and_cons.extend([var, const])
            obj_vars.append(var)
        model.add_cons_vars(vars_and_cons)

        model.objective = prob.Objective(Zero, sloppy=True, direction="max")
        model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})

        sol = model.optimize()

    rxns_to_remove = sol.fluxes[sol.fluxes.abs() < zero_cutoff].index
    consistent_model = model.copy()
    consistent_model.remove_reactions(rxns_to_remove, remove_orphans=True)
    return consistent_model
