import numpy
import pylab


def plot_production_envelope(model, target_id, n_points=20, plot=True):
    target_id = str(target_id)
    target_reaction = model.reactions.get_by_id(target_id)
    original_target_bounds = (target_reaction.lower_bound,
                              target_reaction.upper_bound)
    hot_start = model.optimize()
    max_growth_rate = model.solution.f
    max_growth_production = target_reaction.x
    # TODO - ensure the model solves
    growth_coupled = False
    if max_growth_production > 0:
        growth_coupled = True
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
    growth_rates = production_rates * 0
    # make the objective coefficient what it was before
    target_reaction.objective_coefficient = 0
    for reaction, coefficient in original_objectives.iteritems():
        reaction.objective_coefficient = coefficient
    # calculate the maximum growth rate at each production rate
    for i in range(n_points):
        target_reaction.lower_bound = production_rates[i]
        target_reaction.upper_bound = production_rates[i]
        hot_start = model.optimize(the_problem=hot_start)
        growth_rates[i] = model.solution.f
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


if __name__ == "__main__":
    from cobra.oven.legacyIO import load_pickle
    from cobra.test import ecoli_pickle
    from time import time
    import pylab
    model = load_pickle(ecoli_pickle)
    model.reactions.get_by_id("EX_o2(e)").lower_bound = 0
    for i in ["ABTA", "ACALD", "ACKr", "ATPS4rpp", "F6PA",
              "GLUDy", "LDH_D", "MGSA", "PFL", "TPI"]:
        model.reactions.get_by_id(i).lower_bound = 0
        model.reactions.get_by_id(i).upper_bound = 0
    start = time()
    plot_production_envelope(model, "EX_etoh(e)")
    print "ran in %.2f seconds" % (time() - start)
    pylab.show()
    # calculates in approx 1.2 seconds on 3.4 GHz i7
