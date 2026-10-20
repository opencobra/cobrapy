"""Provides a function to find blocked reactions in a model."""

import logging
from typing import TYPE_CHECKING, List, Optional, Union

from ..core import Configuration, get_solution
from .helpers import normalize_cutoff
from .variability import flux_variability_analysis


if TYPE_CHECKING:
    from pandas import DataFrame

    from cobra import Model, Reaction


logger = logging.getLogger(__name__)
configuration = Configuration()


class BlockedReactionsResult(list[str]):
    """Result of blocked reactions analysis.

    Behaves as a list of fully blocked reaction identifiers for backward
    compatibility. Fully blocked reactions are blocked in both the forward
    and reverse directions. Direction-specific blocked reactions are also
    available through the ``forward_blocked`` and ``reverse_blocked``
    attributes.

    Attributes
    ----------
    forward_blocked : list of str
        Identifiers of reactions that cannot carry positive flux.
    reverse_blocked : list of str
        Identifiers of reactions that cannot carry negative flux.
    """

    def __init__(self, forward_blocked, reverse_blocked):
        """Initialize blocked reactions result."""
        super().__init__(list(set(forward_blocked) & set(reverse_blocked)))
        self.forward_blocked = forward_blocked
        self.reverse_blocked = reverse_blocked


def _warn_near_zero_cutoff_fluxes(flux_span: "DataFrame", zero_cutoff: float) -> None:
    """Warn if FVA results are close to the blocked reaction cutoff."""
    lower_bound = 0.1 * zero_cutoff
    upper_bound = 10 * zero_cutoff

    closest_below = None
    closest_above = None
    for direction in ("minimum", "maximum"):
        for reaction_id, flux in flux_span[direction].dropna().items():
            dir_flux = flux * (1 if direction == "maximum" else -1)
            if dir_flux < lower_bound or dir_flux > upper_bound:
                continue
            flux_info = (dir_flux, reaction_id, direction, flux)
            if dir_flux < zero_cutoff:
                if closest_below is None or dir_flux > closest_below[0]:
                    closest_below = flux_info
            elif closest_above is None or dir_flux < closest_above[0]:
                closest_above = flux_info

    if closest_below is None and closest_above is None:
        return

    formatted_fluxes = []
    for label, flux_info in (("below", closest_below), ("above", closest_above)):
        if flux_info is None:
            continue
        _, reaction_id, direction, flux = flux_info
        formatted_fluxes.append(
            f"closest {label}: {reaction_id} ({direction}={flux:.3g})"
        )

    fluxes = ", ".join(formatted_fluxes)
    logger.warning(
        "Some reactions have flux bounds close to "
        f"zero_cutoff={zero_cutoff}, which may make "
        f"the blocked reaction result sensitive to numerical tolerance: {fluxes}.",
    )


def find_blocked_reactions(
    model: "Model",
    reaction_list: Optional[List[Union["Reaction", str]]] = None,
    loopless: Optional[str] = None,
    zero_cutoff: Optional[float] = None,
    open_exchanges: bool = False,
    abs_flux_clip: Optional[float] = 0.1,
    processes: Optional[int] = None,
) -> "BlockedReactionsResult":
    """Find reactions that cannot carry any flux.

    The question whether or not a reaction is blocked is highly dependent
    on the current exchange reaction settings for a COBRA model. Hence an
    argument is provided to open all exchange reactions.

    Parameters
    ----------
    model : cobra.Model
        The model to analyze.
    reaction_list : list of cobra.Reaction or str, optional
        List of reactions to consider, the default includes all model
        reactions (default None).
    loopless : str, "fastSNP" or "potentials", optional
        If set, only loopless flux distributions are considered when checking
        whether reactions can carry flux. The value is passed to
        :func:`flux_variability_analysis` as its loopless method. Supported
        values are ``"fastSNP"`` and ``"potentials"`` (default None).
        See :func:`flux_variability_analysis` for more details.
    zero_cutoff : float, optional
        Flux value which is considered to effectively be zero. The default
        is set to use `model.tolerance` (default None).
    open_exchanges : bool, optional
        Whether or not to open all exchange reactions to very high flux
        ranges (default False).
    abs_flux_clip : float, optional
        Maximum absolute flux value passed to
        :func:`flux_variability_analysis` (default 0.1).
    processes : int, optional
        The number of parallel processes to run. Can speed up the
        computations if the number of reactions is large. If not explicitly
        passed, it will be set from the global configuration singleton
        (default None).

    Returns
    -------
    BlockedReactionsResult
        A list of fully blocked reaction identifiers. Also has a
        ``forward_blocked`` attribute with reactions that cannot carry
        positive flux and a ``reverse_blocked`` attribute with reactions
        that cannot carry negative flux.

    Notes
    -----
    Sink and demand reactions are left untouched. Please modify them manually.

    The result can be highly sensitive to ``zero_cutoff``. When exploring an
    appropriate cutoff, call :func:`flux_variability_analysis` directly once
    with the desired ``abs_flux_clip`` and analyze the returned flux bounds
    later.

    """
    if loopless not in (None, "potentials", "fastSNP"):
        raise ValueError(
            "The `loopless` argument must be either None, 'potentials', "
            f"or 'fastSNP', not {loopless!r}."
        )

    zero_cutoff = normalize_cutoff(model, zero_cutoff)

    if abs_flux_clip is not None and abs_flux_clip < zero_cutoff:
        raise ValueError(
            "The `abs_flux_clip` argument must be bigger than `zero_cutoff`."
        )

    with model:
        max_bound = max(
            1000.0, max(max(abs(b) for b in r.bounds) for r in model.reactions)
        )

        if open_exchanges:
            for reaction in model.exchanges:
                reaction.bounds = (
                    min(reaction.lower_bound, -max_bound),
                    max(reaction.upper_bound, max_bound),
                )

        if reaction_list is None:
            reaction_list = model.reactions
        else:
            reaction_list = model.reactions.get_by_any(reaction_list)

        if loopless is None:
            # Limit the search space to reactions which have zero flux.
            # If the reactions already carry flux in this solution,
            # the active direction is known to be feasible and only
            # the opposite direction needs to be checked.
            model.slim_optimize()
            solution = get_solution(model, reactions=reaction_list)
            forward_unknown = solution.fluxes[
                solution.fluxes < 10 * abs_flux_clip
            ].index.tolist()
            reverse_unknown = solution.fluxes[
                solution.fluxes > -10 * abs_flux_clip
            ].index.tolist()
            reaction_list = [
                (reaction_id, "maximum") for reaction_id in forward_unknown
            ] + [(reaction_id, "minimum") for reaction_id in reverse_unknown]

        # Run FVA to find reactions where both the minimal and maximal flux
        # are zero (below the cut off).
        flux_span = flux_variability_analysis(
            model,
            reaction_list=reaction_list,
            loopless=loopless,
            fraction_of_optimum=None,
            abs_flux_clip=abs_flux_clip,
            processes=processes,
        )
        _warn_near_zero_cutoff_fluxes(flux_span, zero_cutoff)

        forward_blocked = flux_span[
            flux_span["maximum"].notna() & (flux_span["maximum"] < zero_cutoff)
        ].index.tolist()
        reverse_blocked = flux_span[
            flux_span["minimum"].notna() & (flux_span["minimum"] > -zero_cutoff)
        ].index.tolist()

        return BlockedReactionsResult(forward_blocked, reverse_blocked)
