"""Define the Metabolite class."""

import re
from typing import TYPE_CHECKING, Dict, Optional, Union
from warnings import warn

from ..exceptions import OptimizationError
from ..util.solver import check_solver_status
from ..util.util import format_long_string
from .formula import elements_and_molecular_weights
from .species import Species


if TYPE_CHECKING:
    from optlang.interface import Container
    from pandas import DataFrame

    from cobra.core import Solution
    from cobra.summary.metabolite_summary import MetaboliteSummary

# Numbers are not required because of the |(?=[A-Z])? block. See the
# discussion in https://github.com/opencobra/cobrapy/issues/128 for
# more details.
element_re = re.compile("([A-Z][a-z]?)([0-9.]+[0-9.]?|(?=[A-Z])?)")


class Metabolite(Species):
    """Class for information about metabolite in cobra.Reaction.

    Metabolite is a class for holding information regarding
    a metabolite in a cobra.Reaction object.

    Parameters
    ----------
    id : str
        the identifier to associate with the metabolite
    formula : str
        Chemical formula (e.g. H2O)
    name : str
        A human readable name.
    charge : float
       The charge number of the metabolite
    compartment: str or None
       Compartment of the metabolite.
    """

    # noinspection PyShadowingBuiltins
    def __init__(
        self,
        id: Optional[str] = None,
        formula: Optional[str] = None,
        name: Optional[str] = "",
        charge: Optional[float] = None,
        compartment: Optional[str] = None,
    ) -> None:
        """Initialize Metaboblite cobra Species.

        Parameters
        ----------
        id : str
            the identifier to associate with the metabolite
        formula : str
            Chemical formula (e.g. H2O)
        name : str
            A human readable name.
        charge : float
           The charge number of the metabolite
        compartment: str or None
           Compartment of the metabolite.
        """
        super().__init__(id=id, name=name)
        self.formula = formula
        # because in a Model a metabolite may participate in multiple Reactions
        self.compartment = compartment
        self.charge = charge

        self._bound = 0.0

    def _set_id_with_model(self, value: str) -> None:
        """Set id with value.

        Parameters
        ----------
        value: str
        """
        if value in self.model.metabolites:
            raise ValueError(
                f"The model already contains a metabolite with the id:" f" {value}"
            )
        self.model.constraints[self.id].name = value
        self._id = value
        self.model.metabolites._generate_index()

    @property
    def constraint(self) -> "Container":
        """Get the constraints associated with this metabolite from the solver.

        Returns
        -------
        optlang.<interface>.Containter
            the optlang constraint for this metabolite
        """
        if self.model is not None:
            return self.model.constraints[self.id]

    @property
    def elements(self) -> Optional[Dict[str, Union[int, float]]]:
        """Get dicitonary of elements and counts.

        Dictionary of elements as keys and their count in the metabolite
        as integer. When set, the `formula` property is updated accordingly.

        Returns
        -------
        composition: None or Dict
            A dictionary of elements and counts, where count is int unless it is needed
            to be a float.
            Returns None in case of error.

        """
        tmp_formula = self.formula
        if tmp_formula is None:
            return {}
        # necessary for some old pickles which use the deprecated
        # Formula class
        tmp_formula = str(self.formula)
        # commonly occurring characters in incorrectly constructed formulas
        if "*" in tmp_formula:
            warn(f"invalid character '*' found in formula '{self.formula}'")
            tmp_formula = tmp_formula.replace("*", "")
        if "(" in tmp_formula or ")" in tmp_formula:
            warn(f"invalid formula (has parenthesis) in '{self.formula}'")
            return None
        composition = {}
        parsed = element_re.findall(tmp_formula)
        for element, count in parsed:
            if count == "":
                count = 1
            else:
                try:
                    count = float(count)
                    int_count = int(count)
                    if count == int_count:
                        count = int_count
                    else:
                        warn(f"{count} is not an integer (in formula {self.formula})")
                except ValueError:
                    warn(f"failed to parse {count} (in formula {self.formula})")
                    return None
            if element in composition:
                composition[element] += count
            else:
                composition[element] = count
        return composition

    @elements.setter
    def elements(self, elements_dict: Dict[str, Union[int, float]]) -> None:
        """Update formula based on elements dictionary.

        Parameters
        ----------
        elements_dict: dict
            A dicitonary of elements as keys, count as items.
        """

        def stringify(element, number):
            return element if number == 1 else element + str(number)

        self.formula = "".join(
            stringify(e, n) for e, n in sorted(elements_dict.items())
        )

    @property
    def formula_weight(self) -> Union[int, float]:
        """Calculate the formula weight.

        Returns
        ------
        float, int
            Weight of formula, based on the weight and count of elements. Can be int if
            the formula weight is a whole number, but unlikely.
        """
        try:
            return sum(
                [
                    count * elements_and_molecular_weights[element]
                    for element, count in self.elements.items()
                ]
            )
        except KeyError as e:
            warn(f"The element {e} does not appear in the periodic table")

    @property
    def y(self) -> float:
        """Return the shadow price for the metabolite in the most recent solution.

        Shadow prices are computed from the dual values of the bounds in
        the solution.
        .. deprecated ::
        Use metabolite.shadow_price instead.

        Returns
        -------
        float
            Float representing the shadow price.
        """
        warn("Please use metabolite.shadow_price instead.", DeprecationWarning)
        return self.shadow_price

    @property
    def shadow_price(self) -> float:
        """Return the shadow price for the metabolite in the most recent solution.

        Shadow price is the dual value of the corresponding constraint in the
        model.

        Returns
        -------
        shadow_price: float

        Warnings
        --------
        * Accessing shadow prices through a `Solution` object is the safer,
          preferred, and only guaranteed to be correct way. You can see how to
          do so easily in the examples.
        * Shadow price is retrieved from the currently defined
          `self._model.solver`. The solver status is checked but there are no
          guarantees that the current solver state is the one you are looking
          for.
        * If you modify the underlying model after an optimization, you will
          retrieve the old optimization values.

        Raises
        ------
        RuntimeError
            If the underlying model was never optimized beforehand or the
            metabolite is not part of a model.
        OptimizationError
            If the solver status is anything other than 'optimal'.

        Examples
        --------
        >>> from cobra.io import load_model
        >>> model = load_model("textbook")
        >>> solution = model.optimize()
        >>> model.metabolites.glc__D_e.shadow_price
        -0.09166474637510488
        >>> solution.shadow_prices.glc__D_e
        -0.091664746375104883
        """
        try:
            check_solver_status(self._model.solver.status)
            return self._model.constraints[self.id].dual
        except AttributeError:
            raise RuntimeError(f"metabolite '{self.id}' is not part of a model")
        # Due to below all-catch, which sucks, need to reraise these.
        except (RuntimeError, OptimizationError) as err:
            raise err.with_traceback()
        # Would love to catch CplexSolverError and GurobiError here.
        except Exception as err:
            raise OptimizationError(
                f"Likely no solution exists. Original solver message: {str(err)}."
            ) from err

    def remove_from_model(self, destructive: bool = False) -> None:
        """Remove the association from self.model.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        destructive : bool, default False
            If False then the metabolite is removed from all
            associated reactions.  If True then all associated
            reactions are removed from the Model.
        """
        self._model.remove_metabolites(self, destructive)

    def summary(
        self,
        solution: Optional["Solution"] = None,
        fva: Optional[Union[float, "DataFrame"]] = None,
    ) -> "MetaboliteSummary":
        """Create a summary of the producing and consuming fluxes.

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
        cobra.summary.MetaboliteSummary

        See Also
        --------
        Reaction.summary
        Model.summary

        """
        from cobra.summary import MetaboliteSummary

        return MetaboliteSummary(
            metabolite=self,
            model=self._model,
            solution=solution,
            fva=fva,
        )

    def _repr_html_(self) -> str:
        """Return the metabolite as an HTML string."""
        return f"""
        <table>
            <tr>
                <td><strong>Metabolite identifier</strong></td><td>{self.id}</td>
            </tr><tr>
                <td><strong>Name</strong></td><td>{format_long_string(self.name)}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{id(self):#x}</td>
            </tr><tr>
                <td><strong>Formula</strong></td><td>{self.formula}</td>
            </tr><tr>
                <td><strong>Compartment</strong></td><td>{self.compartment}</td>
            </tr><tr>
                <td><strong>In {len(self.reactions)} reaction(s)</strong></td><td>
                    {format_long_string(", ".join(r.id for r in self.reactions), 200)}
                    </td>
            </tr>
        </table>"""
