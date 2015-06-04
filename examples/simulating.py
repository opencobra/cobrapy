
# coding: utf-8

# # Simulating with FBA

# This example is available as an IPython [notebook](http://nbviewer.ipython.org/github/opencobra/cobrapy/blob/master/documentation_builder/simulating.ipynb).
# 
# Simulations using flux balance analysis can be solved using Model.optimize(). This will maximize or minimize (maximizing is the default) flux through the objective reactions.
# 

# In[1]:

import cobra.test
model = cobra.test.create_test_model()


# ## Running FBA

# In[2]:

model.optimize()


# The Model.optimize() function will return a Solution object, which will also be stored at model.solution. A solution object has several attributes:
# 
#  - f: the objective value
#  - status: the status from the linear programming solver
#  - x_dict: a dictionary of {reaction_id: flux_value} (also called "primal")
#  - x: a list for x_dict
#  - y_dict: a dictionary of {metabolite_id: dual_value}.
#  - y: a list for y_dict

# For example, after the last call to model.optimize(), the status should be 'optimal' if the solver returned no errors, and f should be the objective value

# In[3]:

model.solution.status


# In[4]:

model.solution.f


# ## Changing the Objectives

# The objective function is determined from the objective_coefficient attribute of the objective reaction(s). Currently in the model, there is only one objective reaction, with an objective coefficient of 1.

# In[5]:

model.objective


# The objective function can be changed by assigning Model.objective, which can be a reaction object (or just it's name), or a dict of {Reaction: objective_coefficient}.

# In[6]:

# change the objective to ATPM
# the upper bound should be 1000 so we get the actual optimal value
model.reactions.get_by_id("ATPM").upper_bound = 1000.
model.objective = "ATPM"
model.objective


# In[7]:

model.optimize()


# The objective function can also be changed by setting Reaction.objective_coefficient directly.

# In[8]:

model.reactions.get_by_id("ATPM").objective_coefficient = 0.
model.reactions.get_by_id("biomass_iRR1083_metals").objective_coefficient = 1.
model.objective

