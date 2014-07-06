
## Simulating with FBA

# This example is available as an IPython [notebook](http://nbviewer.ipython.or
# g/github/opencobra/cobrapy/blob/master/documentation_builder/simulating.ipynb
# ).
# 
# Simulations using flux balance analysis can be solved using Model.optimize().
# This will maximize or minimize (maximizing is the default) flux through the
# objective reactions.
# 

import cobra.test
model = cobra.test.create_test_model()


### Running FBA

model.optimize()
# Output:
# <Solution 0.38 at 0x660d990>

# The Model.optimize() function will return a Solution object, which will also
# be stored at model.solution. A solution object has several attributes:
#  - f: the objective value
#  - status: the status from the linear programming solver
#  - x_dict: a dictionary of {reaction_id: flux_value} (also called "primal")
#  - x: a list for x_dict
#  - y_dict: a dictionary of {metabolite_id: dual_value}.
#  - y: a list for y_dict

# For example, after the last call to model.optimize(), the status should be
# 'optimal' if the solver returned no errors, and f should be the objective
# value

model.solution.status
# Output:
# 'optimal'

model.solution.f
# Output:
# 0.38000797227551136

### Changing the Objectives

# The objective function is determined from the objective_coefficient attribute
# of the objective reaction(s). Currently in the model, there is only one
# objective reaction, with an objective coefficient of 1.

{reaction: reaction.objective_coefficient for reaction in model.reactions
        if reaction.objective_coefficient > 0}
# Output:
# {<Reaction biomass_iRR1083_metals at 0x660d350>: 1.0}

# The objective function can be changed by using the function
# Model.change_objective, which will take either a reaction object or just its
# name.

# change the objective to ATPM
# the upper bound should be 1000 so we get the actual optimal value
model.reactions.get_by_id("ATPM").upper_bound = 1000.
model.change_objective("ATPM")
{reaction: reaction.objective_coefficient for reaction in model.reactions
        if reaction.objective_coefficient > 0}
# Output:
# {<Reaction ATPM at 0x52cb190>: 1.0}

model.optimize()
# Output:
# <Solution 119.67 at 0x4c93110>

# The objective function can also be changed by setting
# Reaction.objective_coefficient directly.

model.reactions.get_by_id("ATPM").objective_coefficient = 0.
model.reactions.get_by_id("biomass_iRR1083_metals").objective_coefficient = 1.
{reaction: reaction.objective_coefficient for reaction in model.reactions
        if reaction.objective_coefficient > 0}
# Output:
# {<Reaction biomass_iRR1083_metals at 0x660d350>: 1.0}
