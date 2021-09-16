"""Provide parsimonious FBA implementation."""

from itertools import chain
from typing import TYPE_CHECKING, Callable, Dict, List, Optional, Union
from warnings import warn

from optlang.symbolics import Zero

from ..core.solution import get_solution
from ..util import solver as sutil


if TYPE_CHECKING:
    from optlang.interface import Objective

    from cobra import Model, Reaction, Solution


def optimize_minimal_flux(
    *args, **kwargs
) -> Callable[["Model", float, Union[Dict, "Objective"], List["Reaction"]], "Solution"]:
    """Perform basic pFBA to minimize total flux.

    .. deprecated:: 0.6.0a4
            `optimize_minimal_flux` will be removed in cobrapy 1.0.0, it is
            replaced by `pfba`.

    Parameters
    ----------
    *args: Any
        Non-keyword variable-length arguments.
    **kwargs: Any
        Keyword-only variable-length arguments.

    Returns
    -------
    A function performing the parsimonious FBA.

    """
    warn("optimize_minimal_flux has been renamed to pfba", DeprecationWarning)
    return pfba(*args, **kwargs)


def pfba(
    model: "Model",
    fraction_of_optimum: float = 1.0,
    objective: Union[Dict, "Objective", None] = None,
    reactions: Optional[List["Reaction"]] = None,
) -> "Solution":
    """Perform basic pFBA (parsimonious Enzyme Usage Flux Balance Analysis).

    pFBA [1] adds the minimization of all fluxes the the objective of the
    model. This approach is motivated by the idea that high fluxes have a
    higher enzyme turn-over and that since producing enzymes is costly,
    the cell will try to minimize overall flux while still maximizing the
    original objective function, e.g. the growth rate.

    Parameters
    ----------
    model : cobra.Model
        The model to perform pFBA on.
    fraction_of_optimum : float, optional
        The fraction of optimum which must be maintained. The original
        objective reaction is constrained to be greater than maximal value
        times the `fraction_of_optimum` (default 1.0).
    objective : dict or cobra.Model.objective, optional
        A desired objective to use during optimization in addition to the
        pFBA objective. Dictionaries (reaction as key, coefficient as value)
        can be used for linear objectives (default None).
    reactions : list of cobra.Reaction, optional
        List of cobra.Reaction. Implies `return_frame` to be true. Only
        return fluxes for the given reactions. Faster than fetching all
        fluxes if only a few are needed (default None).

    Returns
    -------
    cobra.Solution
        The solution object to the optimized model with pFBA constraints
        added.

    References
    ----------
    .. [1] Lewis, N. E., Hixson, K. K., Conrad, T. M., Lerman, J. A.,
       Charusanti, P., Polpitiya, A. D., Palsson, B. O. (2010). Omic data
       from evolved E. coli are consistent with computed optimal growth from
       genome-scale models. Molecular Systems Biology, 6,
       390. doi:10.1038/msb.2010.47

    """
    reactions = (
        model.reactions if reactions is None else model.reactions.get_by_any(reactions)
    )
    with model as m:
        add_pfba(m, objective=objective, fraction_of_optimum=fraction_of_optimum)
        m.slim_optimize(error_value=None)
        solution = get_solution(m, reactions=reactions)
    return solution


def add_pfba(
    model: "Model",
    objective: Union[Dict, "Objective", None] = None,
    fraction_of_optimum: float = 1.0,
) -> None:
    """Add pFBA objective to the `model`.

    This adds objective to minimize the summed flux of all reactions to the
    current objective.

    Parameters
    ----------
    model : cobra.Model
        The model to add the objective to.
    objective : dict or cobra.Model.objective, optional
        A desired objective to use during optimization in addition to the
        pFBA objective. Dictionaries (reaction as key, coefficient as value)
        can be used for linear objectives (default None).
    fraction_of_optimum : float, optional
        Fraction of optimum which must be maintained. The original objective
        reaction is constrained to be greater than maximal value times the
        `fraction_of_optimum`.

    See Also
    -------
    pfba

    """
    if objective is not None:
        model.objective = objective
    if model.solver.objective.name == "_pfba_objective":
        raise ValueError("The model already has a pFBA objective.")
    sutil.fix_objective_as_constraint(model, fraction=fraction_of_optimum)
    reaction_variables = (
        (rxn.forward_variable, rxn.reverse_variable) for rxn in model.reactions
    )
    variables = chain(*reaction_variables)
    model.objective = model.problem.Objective(
        Zero, direction="min", sloppy=True, name="_pfba_objective"
    )
    model.objective.set_linear_coefficients({v: 1.0 for v in variables})
