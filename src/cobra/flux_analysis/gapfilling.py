"""Provide the base class and utility function for gap filling."""


import logging
from typing import TYPE_CHECKING, Dict, List, Optional, Union

from optlang.interface import OPTIMAL
from optlang.symbolics import Zero

from ..core import Model
from ..util import fix_objective_as_constraint, interface_to_str


if TYPE_CHECKING:
    from cobra import Reaction


logger = logging.getLogger(__name__)


class GapFiller:
    r"""
    The base class for performing gap filling.

    This class implements gap filling based on a mixed-integer approach,
    very similar to that described in [1]_ and the 'no-growth but growth'
    part of [2]_ but with minor adjustments. In short, we add indicator
    variables for using the reactions in the universal model, z_i and then
    solve problem

    minimize: \sum_i c_i * z_i
    s.t.    : Sv = 0
              v_o \ge t
              lb_i \le v_i \le ub_i
              v_i = 0 if z_i = 0

    where lb, ub are respectively the upper, lower flux bounds for reaction i,
    c_i is a cost parameter and the objective v_o is greater than the lower
    bound t. The default costs are 1 for reactions from the universal model,
    100 for exchange (uptake) reactions added and 1 for added demand reactions.

    Note that this is a mixed-integer linear program and as such will be
    expensive to solve for large models. Consider using alternatives [3]_
    such as CORDA instead [4,5]_ .

    Parameters
    ----------
    model : cobra.Model
        The model to perform gap filling on.
    universal : cobra.Model, optional
        A universal model with reactions that can be used to complete the
        `model` (default None).
    lower_bound : float, optional
        The minimally accepted flux for the objective in the filled model
        (default 0.05).
    penalties : dict, optional
        A dictionary with keys being 'universal' (all reactions included in
        the universal model), 'exchange' and 'demand' (all additionally
        added exchange and demand reactions) for the three reaction types.
        Can also have reaction identifiers for reaction specific costs.
        Defaults are 1, 100 and 1 respectively (default None).
    exchange_reactions : bool, optional
        Consider adding exchange (uptake) reactions for all metabolites
        in the model (default False).
    demand_reactions : bool, optional
        Consider adding demand reactions for all metabolites (default True).
    integer_threshold : float, optional
        The threshold at which a value is considered non-zero (aka
        integrality threshold). If gapfilled models fail to validate,
        you may want to lower this value (default 1E-6).

    Attributes
    ----------
    indicators: list of optlang.interface.Variable
        The list of symbolic indicator variables.
    costs: dict of {optlang.interface.Variable: float}
        The dictionary with symbolic variables as keys and their cost as
        values.

    References
    ----------
    .. [1] Reed, Jennifer L., Trina R. Patel, Keri H. Chen, Andrew R. Joyce,
       Margaret K. Applebee, Christopher D. Herring, Olivia T. Bui, Eric M.
       Knight, Stephen S. Fong, and Bernhard O. Palsson. “Systems Approach
       to Refining Genome Annotation.” Proceedings of the National Academy
       of Sciences 103, no. 46 (2006): 17480–17484.

    .. [2] Kumar, Vinay Satish, and Costas D. Maranas. “GrowMatch: An
       Automated Method for Reconciling In Silico/In Vivo Growth
       Predictions.” Edited by Christos A. Ouzounis. PLoS Computational
       Biology 5, no. 3 (March 13, 2009): e1000308.
       doi:10.1371/journal.pcbi.1000308.

    .. [3] http://opencobra.github.io/cobrapy/tags/gapfilling/

    .. [4] Schultz, André, and Amina A. Qutub. “Reconstruction of
       Tissue-Specific Metabolic Networks Using CORDA.” Edited by Costas D.
       Maranas. PLOS Computational Biology 12, no. 3 (March 4, 2016):
       e1004808. doi:10.1371/journal.pcbi.1004808.

    .. [5] Diener, Christian https://github.com/cdiener/corda

    """

    def __init__(
        self,
        model: Model,
        universal: Optional[Model] = None,
        lower_bound: float = 0.05,
        penalties: Optional[Union[Dict[str, int], Dict["Reaction", int]]] = None,
        exchange_reactions: bool = False,
        demand_reactions: bool = True,
        integer_threshold: float = 1e-6,
        **kwargs,
    ) -> None:
        """Initialize a new GapFiller object.

        Other Parameters
        ----------------
        kwargs :
            Further keyword arguments are passed on to the parent class.

        """
        self.original_model = model
        self.lower_bound = lower_bound
        self.model = model.copy()
        tolerances = self.model.solver.configuration.tolerances
        try:
            tolerances.integrality = integer_threshold
        except AttributeError:
            logger.warning(
                f"The current solver interface {interface_to_str(self.model.problem)} "
                f"doesn't support setting the integrality tolerance."
            )
        # TODO (Midnighter): One could debate how useful it is to compare against this
        #  threshold when it is not supported by the chosen solver.
        self.integer_threshold = integer_threshold
        self.universal = universal.copy() if universal else Model("universal")
        self.penalties = {"universal": 1, "exchange": 100, "demand": 1}
        if penalties is not None:
            self.penalties.update(penalties)
        self.indicators = []
        self.costs = {}
        self.extend_model(exchange_reactions, demand_reactions)
        fix_objective_as_constraint(self.model, bound=lower_bound)
        self.add_switches_and_objective()

    def extend_model(
        self, exchange_reactions: bool = False, demand_reactions: bool = True
    ) -> None:
        """Extend gap filling model.

        Add reactions from universal model and optionally exchange and
        demand reactions for all metabolites in the model to perform
        gap filling on.

        Parameters
        ----------
        exchange_reactions : bool, optional
            Consider adding exchange (uptake) reactions for all metabolites
            in the model (default False).
        demand_reactions : bool, optional
            Consider adding demand reactions for all metabolites
            (default True).

        """
        for rxn in self.universal.reactions:
            rxn.gapfilling_type = "universal"
        new_metabolites = self.universal.metabolites.query(
            lambda metabolite: metabolite not in self.model.metabolites
        )
        self.model.add_metabolites(new_metabolites)
        existing_exchanges = []
        for rxn in self.universal.boundary:
            existing_exchanges = existing_exchanges + [
                met.id for met in list(rxn.metabolites)
            ]

        for met in self.model.metabolites:
            if exchange_reactions:
                # check for exchange reaction in model already
                if met.id not in existing_exchanges:
                    rxn = self.universal.add_boundary(
                        met,
                        type="exchange_smiley",
                        lb=-1000,
                        ub=0,
                        reaction_id=f"EX_{met.id}",
                    )
                    rxn.gapfilling_type = "exchange"
            if demand_reactions:
                rxn = self.universal.add_boundary(
                    met,
                    type="demand_smiley",
                    lb=0,
                    ub=1000,
                    reaction_id=f"DM_{met.id}",
                )
                rxn.gapfilling_type = "demand"

        new_reactions = self.universal.reactions.query(
            lambda reaction: reaction not in self.model.reactions
        )
        self.model.add_reactions(new_reactions)

    def update_costs(self) -> None:
        """Update coefficients for the indicator variables in the objective.

        Done incrementally so that second time the function is called,
        active indicators in the current solutions gets higher cost than the
        unused indicators.

        """
        for var in self.indicators:
            if var not in self.costs:
                self.costs[var] = var.cost
            else:
                if var._get_primal() > self.integer_threshold:
                    self.costs[var] += var.cost
        self.model.objective.set_linear_coefficients(self.costs)

    def add_switches_and_objective(self) -> None:
        """Update gap filling model with switches and indicator objective."""
        constraints = []
        big_m = max(max(abs(b) for b in r.bounds) for r in self.model.reactions)
        prob = self.model.problem
        for rxn in self.model.reactions:
            if not hasattr(rxn, "gapfilling_type"):
                continue
            indicator = prob.Variable(
                name=f"indicator_{rxn.id}", lb=0, ub=1, type="binary"
            )
            if rxn.id in self.penalties:
                indicator.cost = self.penalties[rxn.id]
            else:
                indicator.cost = self.penalties[rxn.gapfilling_type]
            indicator.rxn_id = rxn.id
            self.indicators.append(indicator)

            # if z = 1 v_i is allowed non-zero
            # v_i - Mz <= 0   and   v_i + Mz >= 0
            constraint_lb = prob.Constraint(
                rxn.flux_expression - big_m * indicator,
                ub=0,
                name=f"constraint_lb_{rxn.id}",
                sloppy=True,
            )
            constraint_ub = prob.Constraint(
                rxn.flux_expression + big_m * indicator,
                lb=0,
                name=f"constraint_ub_{rxn.id}",
                sloppy=True,
            )

            constraints.extend([constraint_lb, constraint_ub])

        self.model.add_cons_vars(self.indicators)
        self.model.add_cons_vars(constraints, sloppy=True)
        self.model.objective = prob.Objective(Zero, direction="min", sloppy=True)
        self.model.objective.set_linear_coefficients({i: 1 for i in self.indicators})
        self.update_costs()

    def fill(self, iterations: int = 1) -> List[List["Reaction"]]:
        """Perform the gap filling.

        With every iteration, it solves the model, updates the costs and
        records the used reactions.

        Parameters
        ----------
        iterations : int, optional
            The number of rounds of gap filling to perform. For every
            iteration, the penalty for every used reaction increases
            linearly. This way, the algorithm is encouraged to search for
            alternative solutions which may include previously used
            reactions i.e., with enough iterations pathways including 10
            steps will eventually be reported even if the shortest pathway
            is a single reaction (default 1).

        Returns
        -------
        list of list of cobra.Reaction
            A list of lists where each element is a list of reactions that
            were used to gap fill the model.

        Raises
        ------
        RuntimeError
            If the model fails to be validated (i.e. the original model with
            the proposed reactions added, still cannot get the required flux
            through the objective).

        """
        used_reactions = []
        for _ in range(iterations):
            self.model.slim_optimize(
                error_value=None, message="gap filling optimization failed"
            )
            solution = [
                self.model.reactions.get_by_id(ind.rxn_id)
                for ind in self.indicators
                if ind._get_primal() > self.integer_threshold
            ]
            if not self.validate(solution):
                raise RuntimeError(
                    "Failed to validate gap filled model, "
                    "try lowering the integer threshold."
                )
            used_reactions.append(solution)
            self.update_costs()
        return used_reactions

    def validate(self, reactions: List["Reaction"]) -> bool:
        """Validate the model.

        Parameters
        ----------
        reactions: list of cobra.Reaction
            The reactions to add to the model for validation.

        Returns
        -------
        bool
            Whether the model is valid or not.

        """
        with self.original_model as model:
            mets = [x.metabolites for x in reactions]
            all_keys = set().union(*(d.keys() for d in mets))
            model.add_metabolites(all_keys)
            model.add_reactions(reactions)
            model.slim_optimize()
            return (
                model.solver.status == OPTIMAL
                and model.solver.objective.value >= self.lower_bound
            )


def gapfill(
    model: Model,
    universal: Optional[Model] = None,
    lower_bound: float = 0.05,
    penalties: Optional[Dict[str, "Reaction"]] = None,
    demand_reactions: bool = True,
    exchange_reactions: bool = False,
    iterations: int = 1,
):
    """Perform gap filling on a model.

    Parameters
    ----------
    model : cobra.Model
        The model to perform gap filling on.
    universal : cobra.Model, optional
        A universal model with reactions that can be used to complete the
        model. Only gapfill considering demand and exchange reactions if
        left missing (default None).
    lower_bound : float, optional
        The minimally accepted flux for the objective in the filled model.
        (default 0.05).
    penalties : dict, optional
        A dictionary with keys being 'universal' (all reactions included in
        the universal model), 'exchange' and 'demand' (all additionally
        added exchange and demand reactions) for the three reaction types.
        Can also have reaction identifiers for reaction specific costs.
        Defaults are 1, 100 and 1 respectively (default None).
    exchange_reactions : bool, optional
        Consider adding exchange (uptake) reactions for all metabolites
        in the model (default False).
    demand_reactions : bool, optional
        Consider adding demand reactions for all metabolites (default True).
    iterations : int, optional
        The number of rounds of gap filling to perform. For every iteration,
        the penalty for every used reaction increases linearly. This way,
        the algorithm is encouraged to search for alternative solutions
        which may include previously used reactions i.e., with enough
        iterations pathways including 10 steps will eventually be reported
        even if the shortest pathway is a single reaction (default 1).

    Returns
    -------
    list of list of cobra.Reaction
        A list of lists with on set of reactions that completes the model per
        requested iteration.

    Examples
    --------
    >>> from cobra.io import load_model
    >>> from cobra import Model
    >>> from cobra.flux_analysis import gapfill
    >>> model = load_model("iYS1720")
    >>> universal = Model("universal")
    >>> universal.add_reactions([model.reactions.GF6PTA.copy()])
    >>> model.remove_reactions([model.reactions.GF6PTA])
    >>> gapfill(model, universal)
    [[<Reaction GF6PTA at 0x12206a280>]]

    """
    gapfiller = GapFiller(
        model,
        universal=universal,
        lower_bound=lower_bound,
        penalties=penalties,
        demand_reactions=demand_reactions,
        exchange_reactions=exchange_reactions,
    )
    return gapfiller.fill(iterations=iterations)
