from __future__ import print_function

from ..core import Model, Reaction, Metabolite
from ..solvers import get_solver_name
from ..manipulation import modify


class SUXModelMILP(Model):
    """Model with additional Universal and Exchange reactions.
    Adds corresponding dummy reactions and dummy metabolites for each added
    reaction which are used to impose MILP constraints to minimize the
    total number of added reactions. See the figure for more
    information on the structure of the matrix.
    """

    def __init__(self, model, Universal=None, threshold=.05,
                 penalties=None, dm_rxns=True, ex_rxns=False):
        Model.__init__(self, "")
        # store parameters
        self.threshold = threshold
        if penalties is None:
            self.penalties = {"Universal": 1, "Exchange": 100, "Demand": 1}
        else:
            self.penalties = penalties
        # want to only operate on a copy of Universal so as not to mess up
        # is this necessary?
        if Universal is None:
            Universal = Model("Universal_Reactions")
        else:
            Universal = Universal.copy()

        modify.convert_to_irreversible(Universal)

        for rxn in Universal.reactions:
            rxn.notes["gapfilling_type"] = "Universal"

        # SUX += Exchange (when exchange generator has been written)
        # For now, adding exchange reactions to Universal - could add to a new
        # model called exchange and allow their addition or not....
        if ex_rxns:
            for m in model.metabolites:
                rxn = Reaction('SMILEY_EX_' + m.id)
                rxn.lower_bound = 0
                rxn.upper_bound = 1000
                rxn.add_metabolites({m: 1.0})
                rxn.notes["gapfilling_type"] = "Exchange"
                Universal.add_reaction(rxn)

        if dm_rxns:
            # ADD DEMAND REACTIONS FOR ALL METABOLITES TO UNIVERSAL MODEL
            for m in model.metabolites:
                rxn = Reaction('SMILEY_DM_' + m.id)
                rxn.lower_bound = 0
                rxn.upper_bound = 1000
                rxn.add_metabolites({m: -1.0})
                rxn.notes["gapfilling_type"] = "Demand"
                Universal.add_reaction(rxn)

        Model.add_reactions(self, model.copy().reactions)
        Model.add_reactions(self, Universal.reactions)

        # all reactions with an index < len(model.reactions) were original
        self.original_reactions = self.reactions[:len(model.reactions)]
        self.added_reactions = self.reactions[len(model.reactions):]

        # Add MILP indicator reactions
        indicators = []
        for reaction in self.added_reactions:
            dummy_metabolite = Metabolite("dummy_met_" + reaction.id)
            dummy_metabolite._constraint_sense = "L"
            reaction.add_metabolites({dummy_metabolite: 1})
            indicator_reaction = Reaction("indicator_" + reaction.id)
            indicator_reaction.add_metabolites(
                {dummy_metabolite: -1 * reaction.upper_bound})
            indicator_reaction.lower_bound = 0
            indicator_reaction.upper_bound = 1
            indicator_reaction.variable_kind = "integer"
            indicator_reaction.objective_coefficient = \
                self.penalties[reaction.notes["gapfilling_type"]]
            indicators.append(indicator_reaction)
        Model.add_reactions(self, indicators)

        # original reaction objectives need to be set to lower bounds
        self._update_objectives()

    def _update_objectives(self, added=True):
        """Update the metabolite which encodes the objective function
        with the objective coefficients for the reaction, and impose
        penalties for added reactions.
        """
        for reaction in self.original_reactions:
            if reaction.objective_coefficient > 0:
                reaction.lower_bound = max(
                    reaction.lower_bound,
                    reaction.objective_coefficient * self.threshold)
            reaction.objective_coefficient = 0

    def add_reactions(self, reactions):
        Model.add_reactions(self, reactions)
        self.original_reactions.extend(reactions)
        self._update_objectives()

    def solve(self, solver=None, iterations=1, debug=False, time_limit=100,
              **solver_parameters):
        """solve the MILP problem"""
        if solver is None:
            solver = get_solver_name(mip=True)
        used_reactions = [None] * iterations
        numeric_error_cutoff = 0.0001
        self._update_objectives()
        for i in range(iterations):
            used_reactions[i] = []
            self.optimize(objective_sense="minimize",
                          solver=solver, **solver_parameters)
            if debug:
                print("Iteration %d: Status is %s" % (i, self.solution.status))
            for reaction in self.added_reactions:
                # The dummy reaction should have a flux of either 0 or 1.
                # If it is 1 (nonzero), then the reaction was used in
                # the solution.
                ind = self.reactions.get_by_id("indicator_" + reaction.id)
                if ind.x > numeric_error_cutoff:
                    used_reactions[i].append(reaction)
                    ind.objective_coefficient += \
                        self.penalties[reaction.notes["gapfilling_type"]]
                    if debug:
                        print('    ', reaction, reaction.objective_coefficient)

        return used_reactions


def growMatch(model, Universal, dm_rxns=False, ex_rxns=False,
              penalties=None, **solver_parameters):
    """runs growMatch"""
    SUX = SUXModelMILP(model, Universal, dm_rxns=dm_rxns, ex_rxns=ex_rxns,
                       penalties=penalties)
    return SUX.solve(**solver_parameters)


def SMILEY(model, metabolite_id, Universal,
           dm_rxns=False, ex_rxns=False, penalties=None, **solver_parameters):
    """
    runs the SMILEY algorithm to determine which gaps should be
    filled in order for the model to create the metabolite with the
    given metabolite_id.

    This function is good for running the algorithm once. For more fine-
    grained control, create a SUXModelMILP object, add a demand reaction
    for the given metabolite_id, and call the solve function on the
    SUXModelMILP object.
    """
    SUX = SUXModelMILP(model, Universal, dm_rxns=dm_rxns, ex_rxns=ex_rxns,
                       penalties=penalties)
    # change the objective to be the metabolite
    for reaction in SUX.original_reactions:
        reaction.objective_coefficient = 0
    demand_name = "SMILEY_DM_" + metabolite_id
    if demand_name not in SUX.reactions:
        demand_reaction = Reaction(demand_name)
        demand_reaction.add_metabolites(
            {SUX.metabolites.get_by_id(metabolite_id): -1})
        SUX.add_reaction(demand_reaction)
    else:
        demand_reaction = SUX.reactions.get_by_id(demand_name)
    demand_reaction.lower_bound = SUX.threshold
    return SUX.solve(**solver_parameters)
