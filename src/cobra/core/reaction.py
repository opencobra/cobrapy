"""Define the Reaction class."""

import hashlib
import re
from collections import defaultdict
from copy import copy, deepcopy
from functools import partial
from math import isinf
from operator import attrgetter
from typing import (
    TYPE_CHECKING,
    AnyStr,
    Dict,
    FrozenSet,
    Iterable,
    Iterator,
    List,
    Optional,
    Sequence,
    Set,
    Tuple,
    Union,
)
from warnings import warn


if TYPE_CHECKING:
    from optlang.interface import Variable
    from cobra import Model, Solution
    from cobra.summary import ReactionSummary
    import pandas as pd

from ..exceptions import OptimizationError
from ..manipulation import remove_genes
from ..util.context import get_context, resettable
from ..util.solver import (
    check_solver_status,
    linear_reaction_coefficients,
    set_objective,
)
from ..util.util import format_long_string
from .configuration import Configuration
from .gene import GPR, Gene
from .metabolite import Metabolite
from .object import Object


config = Configuration()

# This regular expression finds any single letter compartment enclosed in
# square brackets at the beginning of the string. For example [c] : foo --> bar
compartment_finder = re.compile(r"^\s*\[([A-Za-z])\]\s*:*")
# Regular expressions to match the arrows
_reversible_arrow_finder = re.compile("<(-+|=+)>")
_forward_arrow_finder = re.compile("(-+|=+)>")
_reverse_arrow_finder = re.compile("<(-+|=+)")


class Reaction(Object):
    """Define the cobra.Reaction class.

    Reaction is a class for holding information regarding
    a biochemical reaction in a cobra.Model object.

    Reactions are by default irreversible with bounds
    `(0.0, cobra.Configuration().upper_bound)`
    if no bounds are provided on creation.
    To create an irreversible reaction use `lower_bound=None`,
    resulting in reaction bounds of
    `(cobra.Configuration().lower_bound, cobra.Configuration().upper_bound)`.

    Parameters
    ----------
    id : str, optional
        The identifier to associate with this reaction (default None).
    name : str, optional
        A human readable name for the reaction (default "").
    subsystem : str, optional
        Subsystem where the reaction is meant to occur (default "").
    lower_bound : float
        The lower flux bound (default 0.0).
    upper_bound : float, optional
        The upper flux bound (default None).
    **kwargs:
        Further keyword arguments are passed on to the parent class.
    """

    # noinspection PyShadowingBuiltins
    def __init__(
        self,
        id: Optional[str] = None,
        name: str = "",
        subsystem: str = "",
        lower_bound: float = 0.0,
        upper_bound: Optional[float] = None,
        **kwargs,
    ) -> None:
        """Initialize the Reaction class."""
        super().__init__(id, name, **kwargs)
        self._gpr = GPR()
        self.subsystem = subsystem

        # The cobra.Genes that are used to catalyze the reaction
        self._genes = set()

        # A dictionary of metabolites and their stoichiometric coefficients in
        # this reaction.
        self._metabolites = {}

        # self.model is None or refers to the cobra.Model that
        # contains self
        self._model = None

        # from cameo ...
        self._lower_bound = (
            lower_bound if lower_bound is not None else config.lower_bound
        )
        self._upper_bound = (
            upper_bound if upper_bound is not None else config.upper_bound
        )

    def _set_id_with_model(self, value: str) -> None:
        """Set Reaction id in model, check that it doesn't already exist.

        The function will rebuild the model reaction index.

        Parameters
        ----------
        value: str
            A string that represents the id.

        Raises
        ------
        ValueError
            If the model already contains a reaction with the id value.
        """
        if value in self.model.reactions:
            raise ValueError(
                f"The model already contains a reaction with the id: {value}"
            )
        forward_variable = self.forward_variable
        reverse_variable = self.reverse_variable
        self._id = value
        self.model.reactions._generate_index()
        forward_variable.name = self.id
        reverse_variable.name = self.reverse_id

    @property
    def reverse_id(self) -> str:
        """Generate the id of reverse_variable from the reaction's id.

        Returns
        -------
        str
            The original id, joined to the word reverse and a partial hash
            of the utf-8 encoded id.
        """
        return "_".join(
            (self.id, "reverse", hashlib.md5(self.id.encode("utf-8")).hexdigest()[0:5])
        )

    @property
    def flux_expression(self) -> Optional["Variable"]:
        """Get Forward flux expression.

        Returns
        -------
        flux_expression: optlang.interface.Variable, optional
            The expression representing the the forward flux (if associated
            with model), otherwise None. Representing the net flux if
            model.reversible_encoding == 'unsplit' or None if reaction is
            not associated with a model
        """
        if self.model is not None:
            return 1.0 * self.forward_variable - 1.0 * self.reverse_variable
        else:
            return None

    @property
    def forward_variable(self) -> Optional["Variable"]:
        """Get an optlang variable representing the forward flux.

        Returns
        -------
        optlang.interface.Variable, optional
            An optlang variable for the forward flux or None if reaction is
            not associated with a model.
        """
        if self.model is not None:
            return self.model.variables[self.id]
        else:
            return None

    @property
    def reverse_variable(self) -> Optional["Variable"]:
        """Get an optlang variable representing the reverse flux.

        Returns
        -------
        optlang.interface.Variable, optional
            An optlang variable for the reverse flux or None if reaction is
            not associated with a model.
        """
        if self.model is not None:
            return self.model.variables[self.reverse_id]
        else:
            return None

    @property
    def objective_coefficient(self) -> float:
        """Get the coefficient for this reaction in a linear objective (float).

        Assuming that the objective of the associated model is summation of
        fluxes from a set of reactions, the coefficient for each reaction
        can be obtained individually using this property. A more general way
        is to use the `model.objective` property directly.

        Returns
        ------
        float
            Linear coefficient if this reaction has any, or 0.0 otherwise.

        Raises
        ------
        AttributeError
            If the model of the reaction is missing (None).
        """
        return linear_reaction_coefficients(self.model, [self]).get(self, 0.0)

    @objective_coefficient.setter
    def objective_coefficient(self, value: float) -> None:
        """Set Objective coefficient.

        Parameters
        ----------
        value: float
            Number to set coefficient to.

        Raises
        ------
        AttributeError
            If reaction does not have model.
        """
        if self.model is None:
            raise AttributeError("cannot assign objective to a missing model")
        if self.flux_expression is not None:
            set_objective(self.model, {self: value}, additive=True)

    def __copy__(self) -> "Reaction":
        """Copy the Reaction.

        Returns
        -------
        Reaction
            A new reaction that is a copy of the original reaction.
        """
        cop = copy(super(Reaction, self))
        return cop

    def __deepcopy__(self, memo: dict) -> "Reaction":
        """Copy the reaction with memo.

        Parameters
        ----------
        memo: dict
            Automatically passed parameter.

        Returns
        -------
        Reaction
            A new reaction that is a deep copy of the original reaction with memo.
        """
        cop = deepcopy(super(Reaction, self), memo)
        return cop

    @staticmethod
    def _check_bounds(lb: float, ub: float) -> None:
        """Check if the lower and upper bounds are valid.

        Parameters
        ----------
        lb: float
            The lower bound.
        ub: float
            The upper bound.

        Raises
        ------
        ValueError
            If the lower bound is higher than upper bound.
        """
        if lb > ub:
            raise ValueError(
                f"The lower bound must be less than or equal to the upper bound "
                f"({lb} <= {ub})."
            )

    def update_variable_bounds(self) -> None:
        """Update and correct variable bounds.

        Sets the forward_variable and reverse_variable bounds based on lower and
        upper bounds. This function corrects for bounds defined as inf or -inf.
        This function will also adjust the associated optlang variables associated
        with the reaction.

        See Also
        -------
        optlang.interface.set_bounds
        """
        if self.model is None:
            return
        # We know that `lb <= ub`.
        if self._lower_bound > 0:
            self.forward_variable.set_bounds(
                lb=None if isinf(self._lower_bound) else self._lower_bound,
                ub=None if isinf(self._upper_bound) else self._upper_bound,
            )
            self.reverse_variable.set_bounds(lb=0, ub=0)
        elif self._upper_bound < 0:
            self.forward_variable.set_bounds(lb=0, ub=0)
            self.reverse_variable.set_bounds(
                lb=None if isinf(self._upper_bound) else -self._upper_bound,
                ub=None if isinf(self._lower_bound) else -self._lower_bound,
            )
        else:
            self.forward_variable.set_bounds(
                lb=0, ub=None if isinf(self._upper_bound) else self._upper_bound
            )
            self.reverse_variable.set_bounds(
                lb=0, ub=None if isinf(self._lower_bound) else -self._lower_bound
            )

    @property
    def lower_bound(self) -> float:
        """Get the lower bound.

        Returns
        -------
        float
            The lower bound of the reaction.
        """
        return self._lower_bound

    @lower_bound.setter
    @resettable
    def lower_bound(self, value: float) -> None:
        """Set the lower bound.

        Parameters
        ----------
        value: float
            The value to set the lower bound.

        Setting the lower bound (float) will also adjust the associated optlang
        variables associated with the reaction.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.

        Raises
        ------
        ValueError
            If lower bound higher than the current upper bound. via _check_bounds.

        See Also
        --------
        _check_bounds
        """
        # Validate bounds before setting them.
        self._check_bounds(value, self._upper_bound)
        self._lower_bound = value
        self.update_variable_bounds()

    @property
    def upper_bound(self) -> float:
        """Get the upper bound.

        Returns
        -------
        float
            The upper bound of the reaction.
        """
        return self._upper_bound

    @upper_bound.setter
    @resettable
    def upper_bound(self, value: float) -> None:
        """Set the upper bound.

        Parameters
        ----------
        value: float
            The value to set the upper bound.

        Setting the upper bound (float) will also adjust the associated optlang
        variables associated with the reaction.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.

        Raises
        ------
        ValueError
            If upper bound lower than the current upper bound. via _check_bounds.

        See Also
        --------
        _check_bounds
        """
        # Validate bounds before setting them.
        self._check_bounds(self._lower_bound, value)
        self._upper_bound = value
        self.update_variable_bounds()

    @property
    def bounds(self) -> Tuple[float, float]:
        """Get or the bounds.

        Returns
        -------
        tuple: lower_bound, upper_bound
            A tuple of floats, representing the lower and upper bound.
        """
        return self.lower_bound, self.upper_bound

    @bounds.setter
    @resettable
    def bounds(self, value: Union[Tuple[float, float], Sequence[float]]) -> None:
        """Set the bounds directly, using a tuple or list.

        Parameters
        ----------
        value: tuple or sequence
            The lower bound and upper bound. Invalid bounds will raise ValueError.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.

        Raises
        ------
        ValueError
            If lower bound higher than upper bound, via _check_bounds.

        """
        lower, upper = value
        # Validate bounds before setting them.
        self._check_bounds(lower, upper)
        self._lower_bound = lower
        self._upper_bound = upper
        self.update_variable_bounds()

    @property
    def flux(self) -> float:
        """
        Get the flux value in the most recent solution.

        Flux is the primal value of the corresponding variable in the model.

        Returns
        -------
        flux: float
            Flux is the primal value of the corresponding variable in the model.

        Warnings
        --------
        * Accessing reaction fluxes through a `Solution` object is the safer,
          preferred, and only guaranteed to be correct way. You can see how to
          do so easily in the examples.
        * Reaction flux is retrieved from the currently defined
          `self._model.solver`. The solver status is checked but there are no
          guarantees that the current solver state is the one you are looking
          for.
        * If you modify the underlying model after an optimization, you will
          retrieve the old optimization values.

        Raises
        ------
        RuntimeError
            If the underlying model was never optimized beforehand or the
            reaction is not part of a model.
        OptimizationError
            If the solver status is anything other than 'optimal'.
        AssertionError
            If the flux value is not within the bounds.

        Examples
        --------
        >>> from cobra.io import load_model
        >>> model = load_model("textbook")
        >>> solution = model.optimize()
        >>> model.reactions.PFK.flux
        7.477381962160283
        >>> solution.fluxes.PFK
        7.4773819621602833
        """
        try:
            check_solver_status(self._model.solver.status)
            return self.forward_variable.primal - self.reverse_variable.primal
        except AttributeError:
            raise RuntimeError(f"reaction '{self.id}' is not part of a model")
        # Due to below all-catch, which sucks, need to reraise these.
        except (RuntimeError, OptimizationError) as err:
            raise err
        # Would love to catch CplexSolverError and GurobiError here.
        except Exception as err:
            raise OptimizationError(
                f"Likely no solution exists. Original solver message: {str(err)}."
            ) from err

    @property
    def reduced_cost(self) -> float:
        """
        Get the reduced cost in the most recent solution.

        Reduced cost is the dual value of the corresponding variable in the
        model.

        Returns
        -------
        reducd_cost: float
            A float representing the reduced cost.

        Warnings
        --------
        * Accessing reduced costs through a `Solution` object is the safer,
          preferred, and only guaranteed to be correct way. You can see how to
          do so easily in the examples.
        * Reduced cost is retrieved from the currently defined
          `self._model.solver`. The solver status is checked but there are no
          guarantees that the current solver state is the one you are looking
          for.
        * If you modify the underlying model after an optimization, you will
          retrieve the old optimization values.

        Raises
        ------
        RuntimeError
            If the underlying model was never optimized beforehand or the
            reaction is not part of a model.
        OptimizationError
            If the solver status is anything other than 'optimal'.

        Examples
        --------
        >>> from cobra.io import load_model
        >>> model = load_model("textbook")
        >>> solution = model.optimize()
        >>> model.reactions.PFK.reduced_cost
        -8.673617379884035e-18
        >>> solution.reduced_costs.PFK
        -8.6736173798840355e-18
        """
        try:
            check_solver_status(self._model.solver.status)
            return self.forward_variable.dual - self.reverse_variable.dual
        except AttributeError:
            raise RuntimeError(f"reaction '{self.id}' is not part of a model")
        # Due to below all-catch, which sucks, need to reraise these.
        except (RuntimeError, OptimizationError) as err:
            raise err
        # Would love to catch CplexSolverError and GurobiError here.
        except Exception as err:
            raise OptimizationError(
                f"Likely no solution exists. Original solver message: {str(err)}."
            ) from err

    # read-only
    @property
    def metabolites(self) -> Dict[Metabolite, float]:
        """Get a dictionary of metabolites and coefficients.

        Returns
        -------
        metaoblites: Dict[Metabolite, float]
            A copy of self._metabolites, which is a dictionary of cobra.Metabolite for
            keys and floats for coeffecieints. Positive coefficient means the reaction
            produces this metabolite, while negative coefficient means the reaction
            consumes this metabolite.
        """
        return self._metabolites.copy()

    @property
    def genes(self) -> FrozenSet:
        """Return the genes of the reaction.

        Returns
        -------
        genes: FrozenSet
        """
        return frozenset(self._genes)

    def update_genes_from_gpr(self) -> None:
        """Update genes of reation based on GPR.

        If the reaction has a model, and new genes appear in the GPR, they will be
        created as Gene() entities and added to the model. If the reaction doesn't have
        a model, genes will be created without a model.

        Genes that no longer appear in the GPR will be removed from the reaction, but
        not the model. If you want to remove them expliclty, use model.remove_genes().
        """
        context = get_context(self)
        if self._gpr.body is not None:
            new_gene_names = self._gpr.genes
        else:
            new_gene_names = set()
        old_genes = self._genes.copy()
        new_genes = set()
        if self._model is None:
            self._genes = {Gene(i) for i in new_gene_names}
        else:
            model_genes = self._model.genes
            self._genes = set()
            for g_id in new_gene_names:
                if not model_genes.has_id(g_id):
                    new_gene = Gene(g_id)
                    new_gene._model = self._model
                    model_genes.append(new_gene)
                    if context:
                        # Remove the gene later
                        context(
                            partial(
                                remove_genes,
                                model=self._model,
                                gene_list=[model_genes.get_by_id(g_id)],
                                remove_reactions=False,
                            )
                        )
                        context(partial(setattr, new_gene, "_model", None))
                        # Maybe should be
                        # context(partial(self._model.genes.__isub__, [new_gene]))
                new_gene = model_genes.get_by_id(g_id)
                self._genes.add(new_gene)
                new_genes.add(new_gene)

        # Make the genes aware that it is involved in this reaction
        for g in self._genes:
            self._associate_gene(g)
            if context:
                context(partial(self._dissociate_gene, g))

        # make the old genes aware they are no longer involved in this reaction
        for g in old_genes.difference(new_genes):
            try:
                self._dissociate_gene(g)
                if context:
                    context(partial(self._associate_gene, g))
            except KeyError:
                warn(f"could not remove old gene {g.id} from reaction {self.id}")
            if g in self._genes:  # if an old gene is still a new gene
                raise Exception("something wrong with sets. Shouldn't happen.")

    @property
    def gene_reaction_rule(self) -> str:
        """See gene reaction rule as str.

        Uses the to_string() method of the GPR class

        Returns
        -------
        str

        """
        return self._gpr.to_string()

    @gene_reaction_rule.setter
    @resettable
    def gene_reaction_rule(self, new_rule: str) -> None:
        """Set a new GPR for the reaction, using a str expression.

        Will use the new GPR to update reaction genes.

        Parameters
        ----------
        new_rule : str
            which will be parsed by the string parser in GPR, GPR.from_string(new_rule).
            It makes a new GPR, and does not modify the existing one.

        See Also
        --------
        update_genes_from_gpr()

        """
        self._gpr = GPR.from_string(new_rule)
        self.update_genes_from_gpr()

    @property
    def gene_name_reaction_rule(self):
        """Display gene_reaction_rule with names intead.

        Do NOT use this string for computation. It is intended to give a
        representation of the rule using more familiar gene names instead of
        the often cryptic ids.

        """
        names = {i.id: i.name for i in self._genes}
        return self._gpr.to_string(names=names)

    @property
    def gpr(self) -> GPR:
        """Return the GPR associated with the reaction.

        Returns
        -------
        gpr: GPR
            The GPR class, see cobra.core.gene.GPR() for details.
        """
        return self._gpr

    @gpr.setter
    @resettable
    def gpr(self, value: GPR) -> None:
        """Set a new GPR for the reaction, using GPR() class.

        Also updates the reaction genes based on GPR.

        Parameters
        ----------
        value : GPR() class.

        See Also
        --------
        cobra.core.gene.GPR()
        update_genes_from_gpr()
        """
        self._gpr = value
        self.update_genes_from_gpr()

    @property
    def functional(self) -> bool:
        """All required enzymes for reaction are functional.

        Returns
        -------
        bool
            True if the gene-protein-reaction (GPR) rule is fulfilled for
            this reaction, or if reaction is not associated to a model,
            otherwise False.
        """
        if self._model:
            return self._gpr.eval(
                {gene.id for gene in self.genes if not gene.functional}
            )
        return True

    @property
    def x(self) -> float:
        """Get the flux through the reaction in the most recent solution.

        Flux values are computed from the primal values of the variables in
        the solution.

        Returns
        -------
        flux: float
            Float representing the flux value.

        .. deprecated ::
        Use reaction.flux instead.
        """
        warn("Please use reaction.flux instead.", DeprecationWarning)
        return self.flux

    @property
    def y(self) -> float:
        """Get the reduced cost of the reaction in the most recent solution.

        Reduced costs are computed from the dual values of the variables in
        the solution.

        Returns
        -------
        flux: float
            Float representing the reduced cost value.

        .. deprecated ::
        Use reaction.reduced_cost instead.
        """
        warn("Please use reaction.reduced_cost instead.", DeprecationWarning)
        return self.reduced_cost

    @property
    def reversibility(self) -> bool:
        """Whether the reaction can proceed in both directions (reversible).

        This is computed from the current upper and lower bounds.

        Returns
        -------
        bool:
            True if the reaction is reversible (lower bound lower than 0 and upper
            bound is higher than 0).

        """
        return self._lower_bound < 0 < self._upper_bound

    @reversibility.setter
    def reversibility(self, value: bool) -> None:
        """Set reversiblitly (ignored).

        At this stage this is not set and is ignored. Present for compatiblity reasons.
        If you'd like to change the reversibility of the reaction change
        lower and upper bounds.

        Parameters
        ----------
        value: bool
            This value is ignored.
        """
        warn("Setting reaction reversibility is ignored")

    @property
    def boundary(self) -> bool:
        """Whether or not this reaction is an exchange reaction.

        Returns
        -------
        bool:   `True` if the reaction has either no products or reactants.
        """
        return len(self.metabolites) == 1 and not (self.reactants and self.products)

    @property
    def model(self) -> Optional["Model"]:
        """Return the model the reaction is a part of.

        Returns
        -------
        model: cobra.Model, optional
            The model this reaction belongs to. None if there is no model associated
            with this reaction.
        """
        return self._model

    def _update_awareness(self) -> None:
        """Update awareness for genes and metaoblites of the reaction.

        Make sure all metabolites and genes that are associated with
        this reaction are aware of it.

        """
        for x in self._metabolites:
            x._reaction.add(self)
        for x in self._genes:
            x._reaction.add(self)

    def remove_from_model(self, remove_orphans: bool = False) -> None:
        """Remove the reaction from a model.

        This removes all associations between a reaction the associated
        model, metabolites and genes.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        remove_orphans : bool
            Remove orphaned genes and metabolites from the model as well (default
            False).
        """
        self._model.remove_reactions([self], remove_orphans=remove_orphans)

    def delete(self, remove_orphans: bool = False) -> None:
        """Remove the reaction from a model.

        This removes all associations between a reaction the associated
        model, metabolites and genes.

        The change is reverted upon exit when using the model as a context.

        .. deprecated ::
        use `reaction.remove_from_model` instead.

        Parameters
        ----------
        remove_orphans : bool
            Remove orphaned genes and metabolites from the model as well (default
            False).
        """
        warn(
            "delete is deprecated. Use reaction.remove_from_model instead",
            DeprecationWarning,
        )
        self.remove_from_model(remove_orphans=remove_orphans)

    def __getstate__(self) -> Dict:
        """Get state for reaction.

        This serializes the reaction object. The GPR will be converted to a string
        to avoid unneccessary copies due to interdependencies of used objects.

        Returns
        -------
        dict
            The state/attributes of the reaction in serilized form.

        """
        state = self.__dict__.copy()
        state["_gpr"] = str(self._gpr)
        return state

    def __setstate__(self, state: Dict) -> None:
        """Set state for reaction.

        Probably not necessary to set _model as the cobra.Model that
        contains self sets the _model attribute for all metabolites and genes
        in the reaction.

        However, to increase performance speed we do want to let the metabolite
        and gene know that they are employed in this reaction

        Parameters
        ----------
        state: dict
            A dictionary of state, where keys are attribute names (str).
            Similar to __dict__.
        """
        # These are necessary for old pickles which store attributes
        # which have since been superceded by properties.
        if "reaction" in state:
            state.pop("reaction")
        if "gene_reaction_rule" in state:
            state["_gene_reaction_rule"] = state.pop("gene_reaction_rule")
        if "lower_bound" in state:
            state["_lower_bound"] = state.pop("lower_bound")
        if "upper_bound" in state:
            state["_upper_bound"] = state.pop("upper_bound")

        # Used for efficient storage in newer cobrapy versions
        if "_gpr" not in state:
            state["_gpr"] = state["_gene_reaction_rule"]
        if type(state["_gpr"]) is str:
            state["_gpr"] = GPR.from_string(state["_gpr"])

        self.__dict__.update(state)
        for x in state["_metabolites"]:
            x._model = self._model
            x._reaction.add(self)
        for x in state["_genes"]:
            x._model = self._model
            x._reaction.add(self)

    def copy(self) -> "Reaction":
        """Copy a reaction.

        The referenced metabolites and genes are also copied.

        Returns
        -------
        cobra.Reaction
            A copy of the Reaction.
        """
        # no references to model when copying
        model = self._model
        self._model = None
        for i in self._metabolites:
            i._model = None
        for i in self._genes:
            i._model = None
        # now we can copy
        new_reaction = deepcopy(self)
        # restore the references
        self._model = model
        for i in self._metabolites:
            i._model = model
        for i in self._genes:
            i._model = model
        return new_reaction

    def __add__(self, other: "Reaction") -> "Reaction":
        """Add two reactions and return a new one.

        The stoichiometry will be the combined stoichiometry of the two
        reactions, and the gene reaction rule will be both rules combined by an
        and. All other attributes (i.e. reaction bounds) will match those of
        the first reaction.

        Does not modify in place.

        Parameters
        ----------
        other: cobra.Reaction
            Another reaction to add to the current one.

        Returns
        -------
        Reaction - new reaction with the added properties.

        """
        new_reaction = self.copy()
        if other == 0:
            return new_reaction
        else:
            new_reaction += other

        return new_reaction

    __radd__ = __add__

    def __iadd__(self, other: "Reaction") -> "Reaction":
        """Add two reactions in place and return the modified first one.

        The stoichiometry will be the combined stoichiometry of the two
        reactions, and the gene reaction rule will be both rules combined by an
        and. All other attributes (i.e. reaction bounds) will match those of
        the first reaction.

        Modifies in place.

        Parameters
        ----------
        other: cobra.Reaction
            Another reaction to add to the current one.

        Returns
        -------
        Reaction - original reaction (self) with the added properties.
        """
        self.add_metabolites(other._metabolites, combine=True)
        rule1 = self.gene_reaction_rule.strip()
        rule2 = other.gene_reaction_rule.strip()
        if rule1 != "" and rule2 != "":
            self.gene_reaction_rule = (
                f"({self.gene_reaction_rule}) and " f"({other.gene_reaction_rule})"
            )
        elif rule1 != "" and rule2 == "":
            self.gene_reaction_rule = rule1
        elif rule1 == "" and rule2 != "":
            self.gene_reaction_rule = rule2
        return self

    def __sub__(self, other: "Reaction") -> "Reaction":
        """Subtract two reactions and return a new one.

        The stoichiometry will be the subtracted stoichiometry of the two
        reactions, and the gene_reaction_rule will be the gene_reaction_rule of the
        first reaction. All other attributes (i.e. reaction bounds) will match those of
        the first reaction.

        Does not modify in place. The name will still be that of the first reaction.

        Parameters
        ----------
        other: Reaction
            The reaction to subtract from self.

        Returns
        -------
        Reaction - new reaction with the added properties.
        """
        new = self.copy()
        new -= other
        return new

    def __isub__(self, other: "Reaction") -> "Reaction":
        """Subtract metabolites of one reaction from another in place.

        The stoichiometry will be the metabolites of self minus the metabolites
         of the other. All other attributes including gene_reaction_rule
         (i.e. reaction bounds) will match those of
        the first reaction.

        Modifies in place and changes the original reaction.

        Parameters
        ----------
        other: Reaction
            The reaction to subtract from self.

        Returns
        -------
        Reaction - self with the subtracted metabolites.
        """
        self.subtract_metabolites(other._metabolites, combine=True)
        return self

    def __imul__(self, coefficient: float) -> "Reaction":
        """Scale coefficients in a reaction by a given value in place.

        E.g. A -> B becomes 2A -> 2B.

        If coefficient is less than zero, the reaction is reversed and the
        bounds are swapped.

        Parameters
        ----------
        coefficient: float
            Value to scale coefficients of metabolites by. If less than zero, reverses
            the reaction.

        Returns
        -------
        Reaction
            Returns the same reaction modified in place.
        """
        self._metabolites = {
            met: value * coefficient for met, value in self._metabolites.items()
        }

        if coefficient < 0:
            self.bounds = (-self.upper_bound, -self.lower_bound)

        if self._model:
            self._model._populate_solver([self])

        context = get_context(self)
        if context:
            context(partial(self._model._populate_solver, [self]))
            context(partial(self.__imul__, 1.0 / coefficient))

        return self

    def __mul__(self, coefficient: float) -> "Reaction":
        """Scale coefficients in a reaction by a given value and return new reaction.

        E.g. A -> B becomes 2A -> 2B.

        If coefficient is less than zero, the reaction is reversed and the
        bounds are swapped.

        Parameters
        ----------
        coefficient: float
            Value to scale coefficients of metabolites by. If less than zero, reverses
            the reaction.

        Returns
        -------
        Reaction
            Returns a new reaction, identical to the original except coefficients.
        """
        new = self.copy()
        new *= coefficient
        return new

    @property
    def reactants(self) -> List[Metabolite]:
        """Return a list of reactants for the reaction.

        Returns
        -------
        list
            A list of the metabolites consudmed (coefficient < 0) by the reaction.
        """
        return [k for k, v in self._metabolites.items() if v < 0]

    @property
    def products(self) -> List[Metabolite]:
        """Return a list of products for the reaction.

        Returns
        -------
        list
            A list of the metabolites produced (coefficient > 0) by the reaction.
        """
        return [k for k, v in self._metabolites.items() if v >= 0]

    def get_coefficient(self, metabolite_id: Union[str, Metabolite]) -> float:
        """Return the stoichiometric coefficient of a metabolite.

        Parameters
        ----------
        metabolite_id : str or cobra.Metabolite

        """
        if isinstance(metabolite_id, Metabolite):
            return self._metabolites[metabolite_id]

        _id_to_metabolites = {m.id: m for m in self._metabolites}
        return self._metabolites[_id_to_metabolites[metabolite_id]]

    def get_coefficients(
        self, metabolite_ids: Iterable[Union[str, Metabolite]]
    ) -> Iterator[float]:
        """Return the stoichiometric coefficients for a list of metabolites.

        Parameters
        ----------
        metabolite_ids : iterable
            Containing ``str`` or ``cobra.Metabolite``s.

        Returns
        -------
        map: Iterable
            Returns the result of map function, which is a map object (an Iterable).

        """
        return map(self.get_coefficient, metabolite_ids)

    def add_metabolites(
        self,
        metabolites_to_add: Dict[Metabolite, float],
        combine: bool = True,
        reversibly: bool = True,
    ) -> None:
        """Add metabolites and stoichiometric coefficients to the reaction.

        If the final coefficient for a metabolite is 0 then it is removed
        from the reaction.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        metabolites_to_add : dict
            Dictionary with metabolite objects or metabolite identifiers as
            keys and coefficients as values. If keys are strings (name of a
            metabolite) the reaction must already be part of a model and a
            metabolite with the given name must exist in the model.

        combine : bool
            Describes behavior if a metabolite already exists in the reaction (default
            True).
            True causes the coefficients to be added.
            False causes the coefficient to be replaced.

        reversibly : bool
            Whether to add the change to the context to make the change
            reversibly or not (primarily intended for internal use). Default is True.

        Raises
        ------
        KeyError
            If the metabolite string id is not in the model.
        ValueError
            If the metabolite key in the dictionary is a string, and there is no model
            for the reaction.
        """
        old_coefficients = self.metabolites
        new_metabolites = []
        _id_to_metabolites = dict([(x.id, x) for x in self._metabolites])

        for metabolite, coefficient in metabolites_to_add.items():
            # Make sure metabolites being added belong to the same model, or
            # else copy them.
            if isinstance(metabolite, Metabolite):
                if (metabolite.model is not None) and (
                    metabolite.model is not self._model
                ):
                    metabolite = metabolite.copy()

            met_id = str(metabolite)
            # If a metabolite already exists in the reaction then
            # just add them.
            if met_id in _id_to_metabolites:
                reaction_metabolite = _id_to_metabolites[met_id]
                if combine:
                    self._metabolites[reaction_metabolite] += coefficient
                else:
                    self._metabolites[reaction_metabolite] = coefficient
            else:
                # If the reaction is in a model, ensure we aren't using
                # a duplicate metabolite.
                if self._model:
                    try:
                        metabolite = self._model.metabolites.get_by_id(met_id)
                    except KeyError as e:
                        if isinstance(metabolite, Metabolite):
                            new_metabolites.append(metabolite)
                        else:
                            # do we want to handle creation here?
                            raise e
                elif isinstance(metabolite, str):
                    # if we want to handle creation, this should be changed
                    raise ValueError(
                        f"Reaction '{self.id}' does not belong to a model. "
                        f"Either add the reaction to a model or use Metabolite objects "
                        f"instead of strings as keys."
                    )
                self._metabolites[metabolite] = coefficient
                # make the metabolite aware that it is involved in this
                # reaction
                metabolite._reaction.add(self)

        # from cameo ...
        model = self.model
        if model is not None:
            model.add_metabolites(new_metabolites)

            for metabolite, coefficient in self._metabolites.items():
                model.constraints[metabolite.id].set_linear_coefficients(
                    {
                        self.forward_variable: coefficient,
                        self.reverse_variable: -coefficient,
                    }
                )

        for metabolite, the_coefficient in list(self._metabolites.items()):
            if the_coefficient == 0:
                # make the metabolite aware that it no longer participates
                # in this reaction
                metabolite._reaction.remove(self)
                self._metabolites.pop(metabolite)

        context = get_context(self)
        if context and reversibly:
            if combine:
                # Just subtract the metabolites that were added
                context(
                    partial(
                        self.subtract_metabolites,
                        metabolites_to_add,
                        combine=True,
                        reversibly=False,
                    )
                )
            else:
                # Reset them with add_metabolites
                mets_to_reset = {
                    key: old_coefficients[model.metabolites.get_by_any(key)[0]]
                    for key in metabolites_to_add.keys()
                }

                context(
                    partial(
                        self.add_metabolites,
                        mets_to_reset,
                        combine=False,
                        reversibly=False,
                    )
                )

    def subtract_metabolites(
        self,
        metabolites: Dict[Metabolite, float],
        combine: bool = True,
        reversibly: bool = True,
    ) -> None:
        """Subtract metabolites from a reaction.

        That means add the metabolites with -1*coefficient. If the final
        coefficient for a metabolite is 0 then the metabolite is removed from
        the reaction.

        Notes
        -----
        * A final coefficient < 0 implies a reactant.
        * The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        metabolites : dict
            Dictionary where the keys are of class Metabolite and the values
            are the coefficients. These metabolites will be added to the
            reaction.

        combine : bool
            Describes behavior if a metabolite already exists in the reaction (default
            True).
            True causes the coefficients to be added.
            False causes the coefficient to be replaced.

        reversibly : bool
            Whether to add the change to the context to make the change
            reversibly or not ,primarily intended for internal use (default True).
        """
        self.add_metabolites(
            {k: -v for k, v in metabolites.items()},
            combine=combine,
            reversibly=reversibly,
        )

    @property
    def reaction(self) -> str:
        """Return Human readable reaction str.

        Returns
        -------
        str
            The reaction in a human readble str.
        """
        return self.build_reaction_string()

    @reaction.setter
    def reaction(self, value: str) -> None:
        """Set reaction from Human readable reaction str.

        Parameters
        -------
        value: str
            The reaction in a human readble str.

        See Also
        --------
        build_reaction_from_string()
        """
        self.build_reaction_from_string(value)

    def build_reaction_string(self, use_metabolite_names: bool = False) -> str:
        """Generate a human readable reaction str.

        Parameters
        ----------
        use_metabolite_names: bool
            Whether to use metabolite names (when True) or metabolite ids (when False,
            default).

        Returns
        -------
        str
            A human readable str.
        """

        def format(number: Union[int, float]) -> str:
            return "" if number == 1 else str(number).rstrip(".") + " "

        id_type = "id"
        if use_metabolite_names:
            id_type = "name"
        reactant_bits = []
        product_bits = []
        for met in sorted(self._metabolites, key=attrgetter("id")):
            coefficient = self._metabolites[met]
            name = str(getattr(met, id_type))
            if coefficient >= 0:
                product_bits.append(format(coefficient) + name)
            else:
                reactant_bits.append(format(abs(coefficient)) + name)

        reaction_string = " + ".join(reactant_bits)
        if not self.reversibility:
            if self.lower_bound < 0 and self.upper_bound <= 0:
                reaction_string += " <-- "
            else:
                reaction_string += " --> "
        else:
            reaction_string += " <=> "
        reaction_string += " + ".join(product_bits)
        return reaction_string

    def check_mass_balance(self) -> Dict[str, float]:
        """Compute mass and charge balance for the reaction.

        Returns
        -------
        dict
            a dict of {element: amount} for unbalanced elements.
            "charge" is treated as an element in this dict
            This should be empty for balanced reactions.

        Raises
        ------
        ValueError
            No elements were found in metabolite.
        """
        reaction_element_dict = defaultdict(int)
        for metabolite, coefficient in self._metabolites.items():
            if metabolite.charge is not None:
                reaction_element_dict["charge"] += coefficient * metabolite.charge
            if metabolite.elements is None:
                raise ValueError(f"No elements found in metabolite {metabolite.id}")
            for element, amount in metabolite.elements.items():
                reaction_element_dict[element] += coefficient * amount
        # filter out 0 values
        return {k: v for k, v in reaction_element_dict.items() if v != 0}

    @property
    def compartments(self) -> Set:
        """Return set of compartments the metabolites are in.

        Returns
        -------
        set
            A set of compartments the metabolites are in.
        """
        return {
            met.compartment for met in self._metabolites if met.compartment is not None
        }

    def get_compartments(self) -> list:
        """List compartments the metabolites are in.

        Returns
        -------
        list
            A list of compartments the metabolites are in.
        .. deprecated ::
        Use reaction.compartments() instead.
        """
        warn("use Reaction.compartments instead", DeprecationWarning)
        return list(self.compartments)

    def _associate_gene(self, cobra_gene: Gene) -> None:
        """Associates a cobra.Gene object with a cobra.Reaction.

        Parameters
        ----------
        cobra_gene : cobra.core.Gene.Gene

        """
        self._genes.add(cobra_gene)
        cobra_gene._reaction.add(self)
        cobra_gene._model = self._model

    def _dissociate_gene(self, cobra_gene: Gene) -> None:
        """Dissociates a cobra.Gene object with a cobra.Reaction.

        Parameters
        ----------
        cobra_gene : cobra.core.Gene.Gene

        """
        self._genes.discard(cobra_gene)
        cobra_gene._reaction.discard(self)

    def knock_out(self) -> None:
        """Knockout reaction by setting its bounds to zero."""
        self.bounds = (0, 0)

    def build_reaction_from_string(
        self,
        reaction_str: str,
        verbose: bool = True,
        fwd_arrow: Optional[AnyStr] = None,
        rev_arrow: Optional[AnyStr] = None,
        reversible_arrow: Optional[AnyStr] = None,
        term_split: str = "+",
    ) -> None:
        """Build reaction from reaction equation reaction_str using parser.

        Takes a string and using the specifications supplied in the optional
        arguments infers a set of metabolites, metabolite compartments and
        stoichiometries for the reaction.  It also infers the reversibility
        of the reaction from the reaction arrow.

        Changes to the associated model are reverted upon exit when using
        the model as a context.

        Parameters
        ----------
        reaction_str : str
            a string containing a reaction formula (equation)
        verbose: bool
            setting verbosity of function (default True)
        fwd_arrow : AnyStr, optional
            Str or bytes that encode forward irreversible reaction arrows (default
            None).
        rev_arrow : AnyStr, optional
            Str or bytes that encode backward irreversible reaction arrows (default
            None).
        reversible_arrow : AnyStr, optional
            Str or bytes that encode reversible reaction arrows (default None).
        term_split : str
            dividing individual metabolite entries (default "+")".

        Raises
        ------
        ValueError
            No arrow found in reaction string.
        """
        # set the arrows
        forward_arrow_finder = (
            _forward_arrow_finder
            if fwd_arrow is None
            else re.compile(re.escape(fwd_arrow))
        )
        reverse_arrow_finder = (
            _reverse_arrow_finder
            if rev_arrow is None
            else re.compile(re.escape(rev_arrow))
        )
        reversible_arrow_finder = (
            _reversible_arrow_finder
            if reversible_arrow is None
            else re.compile(re.escape(reversible_arrow))
        )
        if self._model is None:
            warn("no model found")
            model = None
        else:
            model = self._model
        found_compartments = compartment_finder.findall(reaction_str)
        if len(found_compartments) == 1:
            compartment = found_compartments[0]
            reaction_str = compartment_finder.sub("", reaction_str)
        else:
            compartment = None

        # reversible case
        arrow_match = reversible_arrow_finder.search(reaction_str)
        if arrow_match is not None:
            self.bounds = config.lower_bound, config.upper_bound
        else:  # irreversible
            # try forward
            arrow_match = forward_arrow_finder.search(reaction_str)
            if arrow_match is not None:
                self.bounds = 0, config.upper_bound
            else:
                # must be reverse
                arrow_match = reverse_arrow_finder.search(reaction_str)
                if arrow_match is None:
                    raise ValueError(f"no suitable arrow found in '{reaction_str}'")
                else:
                    self.bounds = config.lower_bound, 0
        reactant_str = reaction_str[: arrow_match.start()].strip()
        product_str = reaction_str[arrow_match.end() :].strip()

        self.subtract_metabolites(self.metabolites, combine=True)

        for substr, factor in ((reactant_str, -1), (product_str, 1)):
            if len(substr) == 0:
                continue
            for term in substr.split(term_split):
                term = term.strip()
                if term.lower() == "nothing":
                    continue
                if " " in term:
                    num_str, met_id = term.split()
                    num = float(num_str.lstrip("(").rstrip(")")) * factor
                else:
                    met_id = term
                    num = factor
                if compartment is not None:
                    met_id += f"[{compartment}]"
                try:
                    met = model.metabolites.get_by_id(met_id)
                except KeyError:
                    if verbose:
                        print(f"unknown metabolite '{met_id}' created")
                    met = Metabolite(met_id, compartment=compartment)
                self.add_metabolites({met: num})

    def summary(
        self,
        solution: Optional["Solution"] = None,
        fva: Optional[Union[float, "pd.DataFrame"]] = None,
    ) -> "ReactionSummary":
        """
        Create a summary of the reaction flux.

        Parameters
        ----------
        solution : cobra.Solution, optional
            A previous model solution to use for generating the summary. If
            ``None``, the summary method will generate a parsimonious flux
            distribution (default None).
        fva : pandas.DataFrame or float, optional
            Whether or not to include flux variability analysis in the output.
            If given, `fva` should either be a previous FVA solution matching the
            model or a float between 0 and 1 representing the fraction of the
            optimum objective to be searched (default None).

        Returns
        -------
        cobra.summary.ReactionSummary

        See Also
        --------
        Metabolite.summary
        Model.summary

        """
        from cobra.summary import ReactionSummary

        return ReactionSummary(
            reaction=self,
            model=self._model,
            solution=solution,
            fva=fva,
        )

    def __str__(self) -> str:
        """Return reaction id and reaction as str.

        Returns
        -------
        str
            A string comprised out of reaction id and reaction.
        """
        return f"{self.id}: {self.build_reaction_string()}"

    def _repr_html_(self) -> str:
        """Generate html representation of reaction.

        Returns
        -------
        str
            HTML representation of the reaction.
        """
        return f"""
        <table>
            <tr>
                <td><strong>Reaction identifier</strong></td><td>{format_long_string(
            self.id, 100)}</td>
            </tr><tr>
                <td><strong>Name</strong></td><td>{format_long_string(
            self.name, 100)}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{f"{id(self):#x}"}</td>
            </tr><tr>
                <td><strong>Stoichiometry</strong></td>
                <td>
                    <p style='text-align:right'>{format_long_string(
            self.build_reaction_string(), 200)}</p>
                    <p style='text-align:right'>{format_long_string(
            self.build_reaction_string(True), 200)}</p>
                </td>
            </tr><tr>
                <td><strong>GPR</strong></td><td>{format_long_string(
            self.gene_reaction_rule, 100)}</td>
            </tr><tr>
                <td><strong>Lower bound</strong></td><td>{self.lower_bound}</td>
            </tr><tr>
                <td><strong>Upper bound</strong></td><td>{self.upper_bound}</td>
            </tr>
        </table>
        """
