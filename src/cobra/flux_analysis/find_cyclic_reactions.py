"""Provides a function to find cyclic reactions in a metabolic model."""

from typing import TYPE_CHECKING, List, Optional, Tuple

import numpy as np
import optlang
from optlang.symbolics import Zero

from ..util import create_stoichiometric_matrix
from ..util import solver as sutil
from .helpers import normalize_cutoff


if TYPE_CHECKING:
    from cobra import Model


def _create_find_cyclic_reactions_problem(
    solver: "optlang.interface",
    s_int: np.ndarray,
    directions_int: np.ndarray,
    zero_cutoff: float,
    bound: float,
) -> Tuple["optlang.interface.Model", List["optlang.interface.Variable"]]:
    """Create an optimization problem to find cyclic reactions.

    This optimization problem includes only internal reactions, assumes
    steady-state constraints, and enforces reaction directions.
    """
    model = solver.Model()

    q_list = []
    for i in range(s_int.shape[1]):
        q = solver.Variable(
            name=f"q_{i}",
            lb=bound * min(0, directions_int[i, 0]),
            ub=bound * max(0, directions_int[i, 1]),
            problem=model,
        )
        q_list.append(q)
    model.add(q_list)

    for idx, row in enumerate(s_int):
        nnz_list = np.flatnonzero(np.abs(row) > zero_cutoff)
        if len(nnz_list) == 0:
            continue
        constraint = solver.Constraint(Zero, lb=0, ub=0, name=f"row_{idx}")
        model.add([constraint], sloppy=True)
        model.constraints[f"row_{idx}"].set_linear_coefficients(
            {q_list[i]: row[i] for i in nnz_list}
        )

    model.update()

    model.objective = solver.Objective(Zero)

    return model, q_list


def find_cyclic_reactions(
    model: "Model",
    zero_cutoff: Optional[float] = None,
    bound: float = 1e3,
    method: str = "optimized",
    required_stop_checks_num: int = 3,
) -> Tuple[List[str], List[Tuple[bool, bool]]]:
    """Find reactions, that can be in a loop in a steady state flux distribution.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to analyze.
    zero_cutoff : float, optional
        The cutoff value to consider a flux as zero.
        The default uses the `model.tolerance` (default None).
    bound : float, optional
        The bound for the reaction fluxes in the optimization problem.
        (default is 1e3).
    method : str, optional
        The method to use for finding cyclic reactions.
        Options are "optimized" (default) or "basic".
        See notes for details.
    required_stop_checks_num : int, optional
        This parameter is used only for the "optimized" method.
        The number of random checks to pass to prove that all cyclic
        reactions were found. (default is 3).

    Returns
    -------
    A tuple containing two lists:
        - A list of reaction IDs that can be part of a loop.
        - A list of tuples indicating the possible directions of
          reactions from the first list in the loop.
          Each tuple contains two boolean values: (can_be_negative, can_be_positive).

    Notes
    -----
    This function considers only the stoichiometric matrix and reaction directions.
    Any other constraints present in the model are ignored.
    * If a reaction is identified as cyclic, there may still be no feasible loop
      when taking into account all bounds and additional constraints.
    * However, if a reaction is identified as non-cyclic, it cannot participate
      in any loop in a steady-state flux distribution.

    The "basic" method for each reaction and direction checks if it can be a part of
    a loop by optimizing linear programming problem.

    The "optimized" method uses a faster randomized approach to firstly find all
    reactions that can be part of a loop and then checks their directions. This method
    usually works at least 2 times faster than the "basic" method.
    The `required_stop_checks_num` parameter is used to descrease the probability
    of missing some cyclic reactions.
    """

    if method not in ["optimized", "basic"]:
        raise ValueError(
            "The `method` parameter must be either 'optimized' or 'basic'."
        )

    if required_stop_checks_num < 1:
        raise ValueError(
            "The `required_stop_checks_num` parameter must be greater than 0."
        )

    zero_cutoff = normalize_cutoff(model, zero_cutoff)

    internal = [i for i, r in enumerate(model.reactions) if not r.boundary]
    s_int = create_stoichiometric_matrix(model)[:, np.array(internal)]
    n = s_int.shape[1]

    bounds_int = np.array([model.reactions[i].bounds for i in internal])
    directions_int = np.sign(bounds_int)

    lp, q_list = _create_find_cyclic_reactions_problem(
        model.problem, s_int, directions_int, zero_cutoff, bound
    )

    candidate_reactions = list(range(n))
    can_positive = [False] * n
    can_negative = [False] * n

    if method == "optimized":
        is_cyclic = [False] * n

        def set_reaction_weights():
            signs = 2 * np.random.randint(low=0, high=2, size=n) - 1
            weights = np.random.uniform(0.5, 1.0, size=n) * signs
            lp.objective.set_linear_coefficients(
                {q_list[i]: weights[i] for i in range(n) if not is_cyclic[i]}
            )

        set_reaction_weights()
        dir_order = ["min", "max"]
        stop_checks_num = 0
        cyclic_reactions_num = 0

        while cyclic_reactions_num < n and stop_checks_num < required_stop_checks_num:
            found_cyclic = False
            reverse_dir_order = False

            for dir in dir_order:
                lp.objective.direction = dir
                lp.optimize()

                sutil.check_solver_status(lp.status)

                if abs(lp.objective.value) > zero_cutoff:
                    remove_coef = {}
                    for i in range(n):
                        if not is_cyclic[i] and abs(q_list[i].primal) > zero_cutoff:
                            is_cyclic[i] = True
                            cyclic_reactions_num += 1
                            remove_coef[q_list[i]] = 0
                            if q_list[i].primal > zero_cutoff:
                                can_positive[i] = True
                            else:
                                can_negative[i] = True

                    found_cyclic = True
                    lp.objective.set_linear_coefficients(remove_coef)
                    stop_checks_num = 0
                    break

                reverse_dir_order = True

            if not found_cyclic:
                set_reaction_weights()
                stop_checks_num += 1

            if reverse_dir_order:
                dir_order.reverse()

        candidate_reactions = [i for i in range(n) if is_cyclic[i]]

    lp.objective = model.problem.Objective(Zero, direction="max")
    for i in candidate_reactions:
        for dir in ["min", "max"]:
            if (dir == "min" and can_negative[i]) or (dir == "max" and can_positive[i]):
                continue

            if directions_int[i, int(dir == "max")] == 0:
                continue

            lp.objective.set_linear_coefficients(
                {q_list[i]: 1.0 if dir == "max" else -1.0}
            )

            lp.optimize()
            sutil.check_solver_status(lp.status)

            if lp.objective.value > zero_cutoff:
                if dir == "max":
                    can_positive[i] = True
                else:
                    can_negative[i] = True

        lp.objective.set_linear_coefficients({q_list[i]: 0.0})

    cyclic_reactions = []
    cyclic_reactions_directions = []
    for i in range(n):
        if can_positive[i] or can_negative[i]:
            cyclic_reactions.append(model.reactions[internal[i]].id)
            cyclic_reactions_directions.append((can_negative[i], can_positive[i]))

    return cyclic_reactions, cyclic_reactions_directions
