"""Provide functions to remove thermodynamically infeasible loops."""

from typing import TYPE_CHECKING, Dict, List, Optional, Union
from warnings import warn

import numpy as np
from optlang.symbolics import Zero

from ..core import get_solution
from ..util import create_stoichiometric_matrix, nullspace
from .fast_snp import nullspace_fast_snp
from .find_cyclic_reactions import find_cyclic_reactions
from .helpers import normalize_cutoff


if TYPE_CHECKING:
    from cobra import Model, Reaction, Solution


def _add_loopless_with_nullspace(
    model: "Model",
    n_int: np.ndarray,
    reactions_to_constrain: List[int],
    max_bound: float,
    zero_cutoff: float,
):
    """Add nullspace-based loopless constraints."""
    prob = model.problem

    # Add indicator variables and new constraints
    to_add = []
    for i, ridx in enumerate(reactions_to_constrain):
        if not (np.abs(n_int[:, i]) > zero_cutoff).any():
            continue

        rxn = model.reactions[ridx]
        # indicator variable a_i
        indicator = prob.Variable(f"indicator_{rxn.id}", type="binary")
        # -M*(1 - a_i) <= v_i <= M*a_i
        on_off_constraint = prob.Constraint(
            rxn.flux_expression - max_bound * indicator,
            lb=-max_bound,
            ub=0,
            name=f"on_off_{rxn.id}",
        )
        # -(max_bound + 1) * a_i + 1 <= G_i <= -(max_bound + 1) * a_i + max_bound
        delta_g = prob.Variable(f"delta_g_{rxn.id}")
        delta_g_range = prob.Constraint(
            delta_g + (max_bound + 1) * indicator,
            lb=1,
            ub=max_bound,
            name=f"delta_g_range_{rxn.id}",
        )
        to_add.extend([indicator, on_off_constraint, delta_g, delta_g_range])

    model.add_cons_vars(to_add)

    # Add nullspace constraints for G_i
    for i, row in enumerate(n_int):
        name = f"nullspace_constraint_{i}"
        nullspace_constraint = prob.Constraint(Zero, lb=0, ub=0, name=name)
        model.add_cons_vars([nullspace_constraint])
        coefs = {
            model.variables[f"delta_g_{model.reactions[ridx].id}"]: row[i]
            for i, ridx in enumerate(reactions_to_constrain)
            if abs(row[i]) > zero_cutoff
        }
        model.constraints[name].set_linear_coefficients(coefs)


def _add_loopless_with_nullspace_directional(
    model: "Model",
    n_int: np.ndarray,
    reactions_to_constrain: List[int],
    flux_threshold: float,
    max_bound: float,
    zero_cutoff: float,
):
    """Add directional nullspace-based loopless constraints."""
    prob = model.problem

    # Add indicator variables and new constraints
    to_add = []
    for ridx in reactions_to_constrain:
        rxn = model.reactions[ridx]

        indicator_maximum = prob.Variable(f"indicator_maximum_{rxn.id}", type="binary")
        indicator_minimum = prob.Variable(f"indicator_minimum_{rxn.id}", type="binary")

        one_direction_constraint = prob.Constraint(
            indicator_maximum + indicator_minimum,
            ub=1,
            name=f"single_nonzero_direction_{rxn.id}",
        )

        # eps * a_i+ <= v_i+ <= M * a_i+
        on_off_constraint_maximum1 = prob.Constraint(
            rxn.forward_variable - max_bound * indicator_maximum,
            ub=0,
            name=f"on_off_maximum1_{rxn.id}",
        )
        on_off_constraint_maximum2 = prob.Constraint(
            rxn.forward_variable - flux_threshold * indicator_maximum,
            lb=0,
            name=f"on_off_maximum2_{rxn.id}",
        )

        # eps * a_i- <= v_i- <= M * a_i-
        on_off_constraint_minimum1 = prob.Constraint(
            rxn.reverse_variable - max_bound * indicator_minimum,
            ub=0,
            name=f"on_off_minimum1_{rxn.id}",
        )
        on_off_constraint_minimum2 = prob.Constraint(
            rxn.reverse_variable - flux_threshold * indicator_minimum,
            lb=0,
            name=f"on_off_minimum2_{rxn.id}",
        )

        to_add.extend(
            [
                indicator_maximum,
                indicator_minimum,
                one_direction_constraint,
                on_off_constraint_maximum1,
                on_off_constraint_maximum2,
                on_off_constraint_minimum1,
                on_off_constraint_minimum2,
            ]
        )

        # a_i+ -> G_i <= -1: G_i <= -(max_bound + 1) * a_i+ + max_bound
        delta_g = prob.Variable(f"delta_g_{rxn.id}")
        delta_g_range_maximum = prob.Constraint(
            delta_g + (max_bound + 1) * indicator_maximum,
            ub=max_bound,
            name=f"delta_g_range_maximum_{rxn.id}",
        )

        # a_i- -> G_i >= 1: G_i >= (max_bound + 1) * a_i- - max_bound
        delta_g_range_minimum = prob.Constraint(
            delta_g - (max_bound + 1) * indicator_minimum,
            lb=-max_bound,
            name=f"delta_g_range_minimum_{rxn.id}",
        )

        to_add.extend([delta_g, delta_g_range_maximum, delta_g_range_minimum])

    model.add_cons_vars(to_add)

    # Add nullspace constraints for G_i
    for i, row in enumerate(n_int):
        name = f"nullspace_constraint_{i}"
        nullspace_constraint = prob.Constraint(Zero, lb=0, ub=0, name=name)
        model.add_cons_vars([nullspace_constraint])
        coefs = {
            model.variables[f"delta_g_{model.reactions[ridx].id}"]: row[i]
            for i, ridx in enumerate(reactions_to_constrain)
            if abs(row[i]) > zero_cutoff
        }
        model.constraints[name].set_linear_coefficients(coefs)


def _add_loopless_with_potentials(
    model: "Model",
    s_int: np.ndarray,
    reactions_to_constrain: List[int],
    max_bound: float,
    zero_cutoff: float,
):
    """Add loopless constraints using metabolite potential variables."""
    prob = model.problem

    # Add indicator variables and new constraints
    to_add = []
    for i in range(s_int.shape[0]):
        to_add.append(prob.Variable(f"potential_{i}"))

    for ridx in reactions_to_constrain:
        rxn = model.reactions[ridx]
        # indicator variable a_i
        indicator = prob.Variable(f"indicator_{rxn.id}", type="binary")
        # -M*(1 - a_i) <= v_i <= M*a_i
        on_off_constraint = prob.Constraint(
            rxn.flux_expression - max_bound * indicator,
            lb=-max_bound,
            ub=0,
            name=f"on_off_{rxn.id}",
        )
        to_add.extend([indicator, on_off_constraint])

    model.add_cons_vars(to_add)

    for i, ridx in enumerate(reactions_to_constrain):
        rxn = model.reactions[ridx]
        col = s_int[:, i]

        name = f"delta_g_range_{i}"
        delta_g_range = prob.Constraint(Zero, lb=1, ub=max_bound, name=name)
        model.add_cons_vars([delta_g_range])

        coefs = {
            model.variables[f"potential_{i}"]: col[i]
            for i in range(s_int.shape[0])
            if abs(col[i]) > zero_cutoff
        }
        coefs[model.variables[f"indicator_{rxn.id}"]] = max_bound + 1
        model.constraints[name].set_linear_coefficients(coefs)


def _add_loopless_with_potentials_directional(
    model: "Model",
    s_int: np.ndarray,
    reactions_to_constrain: List[int],
    flux_threshold: float,
    max_bound: float,
    zero_cutoff: float,
):
    """Add directional loopless constraints using metabolite potential variables."""
    prob = model.problem

    # Add potential and indicator variables
    to_add = []
    for i in range(s_int.shape[0]):
        to_add.append(prob.Variable(f"potential_{i}"))

    for ridx in reactions_to_constrain:
        rxn = model.reactions[ridx]

        indicator_maximum = prob.Variable(f"indicator_maximum_{rxn.id}", type="binary")
        indicator_minimum = prob.Variable(f"indicator_minimum_{rxn.id}", type="binary")

        one_direction_constraint = prob.Constraint(
            indicator_maximum + indicator_minimum,
            ub=1,
            name=f"single_nonzero_direction_{rxn.id}",
        )

        # eps * a_i+ <= v_i+ <= M * a_i+
        on_off_constraint_maximum1 = prob.Constraint(
            rxn.forward_variable - max_bound * indicator_maximum,
            ub=0,
            name=f"on_off_maximum1_{rxn.id}",
        )
        on_off_constraint_maximum2 = prob.Constraint(
            rxn.forward_variable - flux_threshold * indicator_maximum,
            lb=0,
            name=f"on_off_maximum2_{rxn.id}",
        )

        # eps * a_i- <= v_i- <= M * a_i-
        on_off_constraint_minimum1 = prob.Constraint(
            rxn.reverse_variable - max_bound * indicator_minimum,
            ub=0,
            name=f"on_off_minimum1_{rxn.id}",
        )
        on_off_constraint_minimum2 = prob.Constraint(
            rxn.reverse_variable - flux_threshold * indicator_minimum,
            lb=0,
            name=f"on_off_minimum2_{rxn.id}",
        )

        to_add.extend(
            [
                indicator_maximum,
                indicator_minimum,
                one_direction_constraint,
                on_off_constraint_maximum1,
                on_off_constraint_maximum2,
                on_off_constraint_minimum1,
                on_off_constraint_minimum2,
            ]
        )

    model.add_cons_vars(to_add)

    for i, ridx in enumerate(reactions_to_constrain):
        rxn = model.reactions[ridx]
        col = s_int[:, i]

        # a_i+ -> G_i <= -1: G_i <= -(max_bound + 1) * a_i+ + max_bound
        name_maximum = f"delta_g_range_maximum_{rxn.id}"
        delta_g_range_maximum = prob.Constraint(
            Zero,
            ub=max_bound,
            name=name_maximum,
        )
        # a_i- -> G_i >= 1: G_i >= (max_bound + 1) * a_i- - max_bound
        name_minimum = f"delta_g_range_minimum_{rxn.id}"
        delta_g_range_minimum = prob.Constraint(
            Zero,
            lb=-max_bound,
            name=name_minimum,
        )
        model.add_cons_vars([delta_g_range_maximum, delta_g_range_minimum])

        coefs = {
            model.variables[f"potential_{j}"]: col[j]
            for j in range(s_int.shape[0])
            if abs(col[j]) > zero_cutoff
        }
        coefs[model.variables[f"indicator_maximum_{rxn.id}"]] = max_bound + 1
        model.constraints[name_maximum].set_linear_coefficients(coefs)

        coefs.pop(model.variables[f"indicator_maximum_{rxn.id}"])
        coefs[model.variables[f"indicator_minimum_{rxn.id}"]] = -(max_bound + 1)
        model.constraints[name_minimum].set_linear_coefficients(coefs)


def add_loopless(
    model: "Model",
    zero_cutoff: Optional[float] = None,
    method: str = "fastSNP",
    reactions: Optional[List[str]] = None,
    flux_threshold: Optional[float] = None,
) -> None:
    """Modify a model so all feasible flux distributions are loopless.

    It adds variables and constraints to a model which will disallow flux
    distributions with loops. This function *will* modify your model.

    If `method` is set to "original" or "fastSNP" the used formulation
    is described in [1]_. If `method` is set to "fastSNP", it uses a
    faster implementation based on the Fast-SNP algorithm [2]_.

    If `method` is set to "potentials", it uses metabolite potential
    variables instead of nullspace-based constraints.

    In most cases you probably want to use the much faster
    `loopless_solution`. May be used in cases where you want to add complex
    constraints and objectives (for instance quadratic objectives) to the
    model afterwards or use an approximation of Gibbs free energy directions
    in your model.

    Parameters
    ----------
    model : cobra.Model
        The model to which to add the constraints.
    zero_cutoff : positive float, optional
        Cutoff used for null space. Coefficients with an absolute value
        smaller than `zero_cutoff` are considered to be zero. The default
        uses the `model.tolerance` (default None).
    method : str, "original", "fastSNP", or "potentials", optional
        The method used to add loopless constraints. The "original" method
        uses the original nullspace formulation from [1]_, while "fastSNP"
        uses a faster nullspace implementation based on the FastSNP
        algorithm. The "potentials" method adds constraints based on
        metabolite potential variables. The "fastSNP" and "potentials"
        methods are much faster in most cases, with relative performance
        depending on the model and optimization problem.
    reactions : list of str, optional
        The list of reaction IDs to constrain. All cycles within these
        reactions will be removed. If `None`, all reactions will be constrained.
    flux_threshold : float, optional
        Minimum flux required when a directional indicator variable is active.
        If provided, separate forward and reverse indicator variables are
        added for each constrained reaction. This is intended for analyses
        that need to distinguish feasible loopless directions.

    Notes
    -----
    When `flux_threshold` is provided, the directional formulation strongly
    relies on binary indicator variables. If the product of the maximum model
    bound and the solver integrality tolerance is close to `flux_threshold`,
    numerical tolerances can make inactive directions appear feasible.
    Increasing `flux_threshold`, decreasing the maximum model bound,
    lowering the solver integrality tolerance, or using a different MILP
    solver can reduce this risk.

    References
    ----------
    .. [1] Elimination of thermodynamically infeasible loops in steady-state
       metabolic models. Schellenberger J, Lewis NE, Palsson BO. Biophys J.
       2011 Feb 2;100(3):544-53. doi: 10.1016/j.bpj.2010.12.3707. Erratum
       in: Biophys J. 2011 Mar 2;100(5):1381.
    .. [2] Fast-SNP: a fast matrix pre-processing algorithm for efficient
       loopless flux optimization of metabolic models. Saa PA, Nielsen LK.
       Bioinformatics. 2016 Dec;32(24):3807–3814. doi: 10.1093/bioinformatics/btw555.
    """
    if method not in ["original", "fastSNP", "potentials"]:
        raise ValueError(f"unsupported method: {method}")

    zero_cutoff = normalize_cutoff(model, zero_cutoff)

    if reactions is None and method != "original":
        reactions = find_cyclic_reactions(model, zero_cutoff=zero_cutoff)[0]

    reactions_to_constrain = [
        i for i, r in enumerate(model.reactions) if not r.boundary
    ]
    if reactions is not None:
        reactions_set = set(reactions)
        reactions_to_constrain = [
            i
            for i, r in enumerate(model.reactions)
            if not r.boundary and r.id in reactions_set
        ]

    s_int = create_stoichiometric_matrix(model)[:, np.array(reactions_to_constrain)]
    s_int = s_int[np.sum(np.abs(s_int) > zero_cutoff, -1) > 0, :]

    max_bound = max(1000.0, max(max(abs(b) for b in r.bounds) for r in model.reactions))

    if flux_threshold is not None:
        try:
            integrality_tolerance = model.solver.configuration.tolerances.integrality
        except AttributeError:
            integrality_tolerance = None

        if (
            integrality_tolerance is None
            or max_bound * integrality_tolerance >= flux_threshold * 0.5
        ):
            warn(
                "Loopless constraints may not work properly "
                f"with the provided `flux_threshold`={flux_threshold}, "
                f"maximum model bound={max_bound}, and solver integrality "
                f"tolerance={integrality_tolerance}. "
                "This can happen due to numerical instability. "
                "Possible remedies are increasing `flux_threshold`, "
                "decreasing the maximum model bound, "
                "switching to a different solver, "
                "or decreasing the solver `integrality` tolerance. "
                "Please carefully read the note on numerical instability "
                "in the `add_loopless` function documentation.",
                UserWarning,
            )

    if method == "potentials":
        if flux_threshold is not None:
            _add_loopless_with_potentials_directional(
                model,
                s_int,
                reactions_to_constrain,
                flux_threshold,
                max_bound,
                zero_cutoff,
            )
        else:
            _add_loopless_with_potentials(
                model,
                s_int,
                reactions_to_constrain,
                max_bound,
                zero_cutoff,
            )

    elif method in ["original", "fastSNP"]:
        if method == "original":
            n_int = nullspace(s_int).T
        else:
            bounds_int = np.array(
                [model.reactions[i].bounds for i in reactions_to_constrain]
            )
            directions_int = np.sign(bounds_int)
            n_int = nullspace_fast_snp(
                model.problem,
                s_int,
                directions_int,
                zero_cutoff=zero_cutoff,
            ).T

        if flux_threshold is not None:
            _add_loopless_with_nullspace_directional(
                model,
                n_int,
                reactions_to_constrain,
                flux_threshold,
                max_bound,
                zero_cutoff,
            )
        else:
            _add_loopless_with_nullspace(
                model,
                n_int,
                reactions_to_constrain,
                max_bound,
                zero_cutoff,
            )

    else:
        raise ValueError(f"unsupported method: {method}")


def _add_cycle_free(model: "Model", fluxes: Dict[str, float]) -> None:
    """Add constraints for CycleFreeFlux.

    Parameters
    ----------
    model : cobra.Model
        The model to operate on.
    fluxes : dict of {str: float}
        A dictionary having keys as reaction IDs and values as their flux
        values.

    """
    model.objective = model.solver.interface.Objective(
        Zero, direction="min", sloppy=True
    )
    objective_vars = []
    for rxn in model.reactions:
        flux = fluxes[rxn.id]
        if rxn.boundary:
            rxn.bounds = (flux, flux)
            continue
        if flux >= 0:
            rxn.bounds = max(0, rxn.lower_bound), min(flux, rxn.upper_bound)
            objective_vars.append(rxn.forward_variable)
        else:
            rxn.bounds = max(flux, rxn.lower_bound), min(0, rxn.upper_bound)
            objective_vars.append(rxn.reverse_variable)

    model.objective.set_linear_coefficients({v: 1.0 for v in objective_vars})


def loopless_solution(
    model: "Model", fluxes: Optional[Dict[str, float]] = None
) -> "Solution":
    """Convert an existing solution to a loopless one.

    Removes as many loops as possible (see Notes).

    Uses the method from CycleFreeFlux [1]_ and is much faster than
    `add_loopless` and should therefore be the preferred option to get
    loopless flux distributions.

    Parameters
    ----------
    model : cobra.Model
        The model to which to add the constraints.
    fluxes : dict of {str, float}, optional
        A dictionary having keys as reaction IDs and values as their flux
        values. If not None will use the provided flux values to obtain a
        close loopless solution (default None).

    Returns
    -------
    cobra.Solution
        A solution object containing the fluxes with the least amount of
        loops possible or None if the optimization failed (usually happening
        if the flux distribution in `fluxes` is infeasible).

    Notes
    -----
    The returned flux solution has the following properties:

    - It contains the minimal number of loops possible and no loops at all
      if all flux bounds include zero and the objective is not in a cycle.
    - It has the same objective value as the original flux solution and assumes
      that the objective does not participate in a cycle
      (which is usually true since it consumes metabolites).
    - It has the same exact exchange fluxes as the previous solution.
    - All fluxes have the same sign (flow in the same direction) as the
      previous solution.

    When providing fluxes to the method, please note that those have to come from the
    exact same model that you provided, meaning that no bounds or coefficients have
    been changed, and the optimum has remained the same.

    References
    ----------
    .. [1] CycleFreeFlux: efficient removal of thermodynamically infeasible
       loops from flux distributions. Desouki AA, Jarre F, Gelius-Dietrich
       G, Lercher MJ. Bioinformatics. 2015 Jul 1;31(13):2159-65. doi:
       10.1093/bioinformatics/btv096.

    """
    # Need to reoptimize otherwise spurious solution artifacts can cause
    # all kinds of havoc
    # TODO: check solution status
    if fluxes is None:
        sol = model.optimize(objective_sense=None)
        fluxes = sol.fluxes
        opt = sol.objective_value
    else:
        opt = model.slim_optimize()

    with model:
        prob = model.problem
        # Fix the objective
        loopless_obj_constraint = prob.Constraint(
            model.objective.expression,
            lb=opt,
            name="loopless_obj_constraint",
        )
        model.add_cons_vars([loopless_obj_constraint])
        _add_cycle_free(model, fluxes)
        solution = model.optimize(objective_sense=None)
        solution.objective_value = loopless_obj_constraint.primal

    return solution


def loopless_fva_iter(
    model: "Model",
    reaction: "Reaction",
    solution: bool = False,
    zero_cutoff: Optional[float] = None,
) -> Union[float, Dict[str, float]]:
    """Plugin to get a loopless FVA solution from single FVA iteration.

    Assumes the following about `model` and `reaction`:
    1. The model objective is set to be `reaction`.
    2. The model has been optimized and contains the minimum/maximum flux
       for `reaction`.
    3. The model contains an auxiliary variable called "fva_old_objective"
       denoting the previous objective.

    Parameters
    ----------
    model : cobra.Model
        The model to be used.
    reaction : cobra.Reaction
        The reaction currently minimized/maximized.
    solution : bool, optional
        Whether to return the entire solution or only the minimum/maximum
        for `reaction` (default False).
    zero_cutoff : positive float, optional
        Cutoff used for loop removal. Fluxes with an absolute value smaller
        than `zero_cutoff` are considered to be zero. The default is to use
        `model.tolerance` (default None).

    Returns
    -------
    single float or dict of {str: float}
        Returns the minimized/maximized flux through `reaction` if
        `solution` is False. Otherwise, returns a loopless flux
        solution object containing the minimum/maximum flux for `reaction`.

    """
    zero_cutoff = normalize_cutoff(model, zero_cutoff)

    current = model.objective.value
    sol = get_solution(model)
    objective_dir = model.objective.direction

    # Handle a suddenly infeasible solution,
    # usually due to numerical instability
    if current is None:
        return None

    # boundary reactions can not be part of cycles
    if reaction.boundary:
        if solution:
            return sol
        else:
            return current

    with model:
        _add_cycle_free(model, sol.fluxes)
        model.slim_optimize()

        # If the previous optimum is maintained in the loopless solution it was
        # loopless and we are done
        if abs(reaction.flux - current) < zero_cutoff:
            if solution:
                return sol
            return current

        # If previous optimum was not in the loopless solution create a new
        # almost loopless solution containing only loops including the current
        # reaction. Than remove all of those loops.
        ll_sol = get_solution(model).fluxes
        reaction.bounds = (current, current)
        model.slim_optimize()
        almost_ll_sol = get_solution(model).fluxes

    with model:
        # find the reactions with loops using the current reaction and remove
        # the loops
        for rxn in model.reactions:
            rid = rxn.id
            if (abs(ll_sol[rid]) < zero_cutoff) and (
                abs(almost_ll_sol[rid]) > zero_cutoff
            ):
                rxn.bounds = max(0, rxn.lower_bound), min(0, rxn.upper_bound)

        if solution:
            best = model.optimize()
        else:
            model.slim_optimize()
            best = reaction.flux
    model.objective.direction = objective_dir
    return best
