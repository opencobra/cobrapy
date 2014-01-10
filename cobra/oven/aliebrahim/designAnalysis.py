import numpy
import pylab

#from ... import solvers
from cobra import solvers
from itertools import combinations


def plot_production_envelope(model, target_id, n_points=20, plot=True,
        solver_name="glpk"):
    """Plot the production envelope for the model given a target

    Parameters
    ----------
    model : cobra model
        The cobra model should already have the uptake rates se
    target_id : str
        The id of the exchange reaction for the target compound
    n_points : int
        The number of points to calculate for the production envolope
    plot : bool, optional
        Whether or not a plot should be made of the production envelope

    Returns
    -------
    growth_rates : :class:`numpy.ndarray`
        An array of growth rates
    production_rates : :class:`numpy.ndarray`
        An array of the corresponding maximum production rate at the
        given growth rate.

    """
    solver = solvers.solver_dict[solver_name]
    target_id = str(target_id)
    target_reaction = model.reactions.get_by_id(target_id)
    original_target_bounds = (target_reaction.lower_bound,
                              target_reaction.upper_bound)
    lp = solver.create_problem(model)
    if solver.solve_problem(lp) != "optimal":
        return ([0], [0])
    solution = solver.format_solution(lp, model)
    max_growth_rate = solution.f
    max_growth_production = solution.x_dict[target_reaction.id]
    #growth_coupled = max_growth_production > 0
    # extract the current objective so it can be changed
    original_objectives = {}
    for reaction in model.reactions:
        if reaction.objective_coefficient != 0:
            original_objectives[reaction] = reaction.objective_coefficient
            reaction.objective_coefficient = 0
    # calculate the maximum possible production rate
    target_reaction.objective_coefficient = 1
    model.optimize(objective_sense="minimize")
    min_production_rate = model.solution.f
    model.optimize(objective_sense="maximize")
    max_production_rate = model.solution.f
    production_rates = numpy.linspace(min_production_rate,
        max_production_rate, n_points)
    # ensure the point of production at maximum growth is included
    production_rates[
        numpy.abs(production_rates - max_growth_production).argmin()] = \
        max_growth_production
    # if the 0 point was overwritten in the last operation
    if production_rates[0] != 0:
        production_rates[1] = production_rates[0]
        production_rates[0] = 0
    growth_rates = production_rates * 0
    # make the objective coefficient what it was before
    target_reaction.objective_coefficient = 0
    for reaction, coefficient in original_objectives.iteritems():
        reaction.objective_coefficient = coefficient
    # calculate the maximum growth rate at each production rate
    for i in range(n_points):
        target_reaction.lower_bound = production_rates[i]
        target_reaction.upper_bound = production_rates[i]
        solver.update_problem(lp, model)
        if solver.solve_problem(lp) == "optimal":
            growth_rates[i] = solver.get_objective_value(lp)
        else:
            growth_rates[i] = 0
    # reset the bounds on the target reaction
    target_reaction.lower_bound = original_target_bounds[0]
    target_reaction.upper_bound = original_target_bounds[1]
    if plot:
        pylab.plot(growth_rates, production_rates)
        pylab.title("Production envelope for %s" % (target_id))
        pylab.xlabel("Growth rate")
        pylab.ylabel("Production rate")
        pylab.xlim(xmin=0)
        pylab.ylim(ymin=0)
    return (growth_rates, production_rates)


def analyze_growth_coupled_num_knockouts(model, knockout_reaction, target_name="EX_etoh_e"):
    None
    

def analyze_growth_coupled_design_subset(model, knockout_reactions, knockout_count, target_name="EX_etoh_e"):
    lp = model.optimize()
    best_score = 0
    best = []
    lb = [None] * k  # store lower bounds when reactions are knocked out
    ub = [None] * k  # store upper bounds when reactions are knocked out
    for subset in combinations(knockout_reactions, knockout_reactions):
        # knockout reactions
        for i, reaction_name in enumerate(subset):
            reaction = model.reactions.get_by_id(str(reaction_name))
            (lb[i], ub[i]) = (reaction.lower_bound, reaction.upper_bound)
            (reaction.lower_bound, reaction.upper_bound) = (0.0, 0.0)
        model.optimize()
        production = model.solution.x_dict[target_name]
        # identical performance
        if abs(production - best_score) < 0.001:
            best.append(subset)
        # better performance
        elif production > best_score:
            best_score = model.solution.x_dict[target_name]
            best = [subset]
        print model.solution.f, model.solution.x_dict[target_name]
        # reset reactions
        for i, reaction_name in enumerate(subset):
            (reaction.lower_bound, reaction.upper_bound) = (lb[i], ub[i])
    return best_score, best

if __name__ == "__main__":
    from cobra.test import ecoli_pickle, create_test_model
    from time import time

    model = create_test_model(ecoli_pickle)
    #from IPython import embed; embed()
    model.reactions.get_by_id("EX_o2_e").lower_bound = 0
    #analyze_strain_design(model, ["ABTA", "ACALD", "ACKr", "ATPS4rpp", "F6PA",
    #          "GLUDy", "LDH_D", "MGSA", "PFL", "TPI"])

    for i in ["ABTA", "ACALD", "ACKr", "ATPS4rpp", "F6PA",
              "GLUDy", "LDH_D", "MGSA", "PFL", "TPI"]:
        model.reactions.get_by_id(i).lower_bound = 0
        model.reactions.get_by_id(i).upper_bound = 0
    start = time()
    plot_production_envelope(model, "EX_etoh_e", solver_name="glpk", n_points=40, plot=True)
    print "ran in %.2f seconds" % (time() - start)
    pylab.show()
    # calculates in approx 1 seconds on 3.4 GHz i7
