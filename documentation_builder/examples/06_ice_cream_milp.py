# Advanced user example showing how to set up and solve an MILP
#

from cobra import Model, Metabolite, Reaction
#solver = 'cplex'  #With libglpk-java there is an untraced memory bug.

cone_selling_price = 7.
cone_production_cost = 3.
popsicle_selling_price = 2.
popsicle_production_cost = 1.
starting_budget = 100.

print()
print(('I can sell cones for $%1.2f.'%(cone_selling_price)))
print(('Cones cost me $%1.2f to produce.'%(cone_production_cost)))
print(('I can sell popsicles for $%1.2f.'%(popsicle_selling_price)))
print(('Popsicles cost me $%1.2f to produce.'%(popsicle_production_cost)))
print(('My total budget was capped at $%1.2f today.'%(starting_budget)))

# problem is:
# max profit
# s.t.
# cone_production_cost*cone_production + popsidle_production_cost*popsicle_production <= starting_budget
# number of cones and popsicles has to be integer...

# first, we'll solve the continuous case just to make sure everything is
# working (it should only make cones)...then we'll tighten the constraints to
# integer... and it should make popsicles.

cobra_model = Model('MILP_implementation_test')
cone_out = Metabolite(id='cone_out', compartment='c')
cone_in = Metabolite(id='cone_in', compartment='c')
cone_consumed = Metabolite(id='cone_consumed', compartment='c')

popsicle_out = Metabolite(id='popsicle_out', compartment='c')
popsicle_in = Metabolite(id='popsicle_in', compartment='c')
popsicle_consumed = Metabolite(id='popsicle_consumed', compartment='c')

the_reactions = []

# SOURCE
Cone_source = Reaction(name='Cone_source')
temp_metabolite_dict = {cone_out: 1}
Cone_source.add_metabolites(temp_metabolite_dict)
the_reactions.append(Cone_source)

Popsicle_source = Reaction(name='Popsicle_source')
temp_metabolite_dict = {popsicle_out: 1}
Popsicle_source.add_metabolites(temp_metabolite_dict)
the_reactions.append(Popsicle_source)


## PRODUCTION
Cone_production = Reaction(name='Cone_production')
temp_metabolite_dict = {cone_out: -1,
                        cone_in: 1}
Cone_production.add_metabolites(temp_metabolite_dict)
the_reactions.append(Cone_production)


Popsicle_production = Reaction(name='Popsicle_production')
temp_metabolite_dict = {popsicle_out: -1,
                        popsicle_in: 1}
Popsicle_production.add_metabolites(temp_metabolite_dict)
the_reactions.append(Popsicle_production)

## CONSUMPTION
Cone_consumption = Reaction(name='Cone_consumption')
temp_metabolite_dict = {cone_in: -1,
                        cone_consumed: 1}
Cone_consumption.add_metabolites(temp_metabolite_dict)
the_reactions.append(Cone_consumption)

Popsicle_consumption = Reaction(name='Popsicle_consumption')
temp_metabolite_dict = {popsicle_in: -1,
                        popsicle_consumed: 1}
Popsicle_consumption.add_metabolites(temp_metabolite_dict)
the_reactions.append(Popsicle_consumption)

# SINK
Cone_consumed_sink = Reaction(name='Cone_consumed_sink')
temp_metabolite_dict = {cone_consumed: -1}
Cone_consumed_sink.add_metabolites(temp_metabolite_dict)
the_reactions.append(Cone_consumed_sink)

Popsicle_consumed_sink = Reaction(name='Popsicle_consumed_sink')
temp_metabolite_dict = {popsicle_consumed: -1}
Popsicle_consumed_sink.add_metabolites(temp_metabolite_dict)
the_reactions.append(Popsicle_consumed_sink)

## add all reactions
cobra_model.add_reactions(the_reactions)

# set objective coefficients
Cone_consumption.objective_coefficient = cone_selling_price
Popsicle_consumption.objective_coefficient = popsicle_selling_price

Cone_production.objective_coefficient = -1*cone_production_cost
Popsicle_production.objective_coefficient = -1*popsicle_production_cost


production_capacity_constraint = Metabolite(id='production_capacity_constraint')
production_capacity_constraint._constraint_sense = 'L'
production_capacity_constraint._bound = starting_budget;

Cone_production.add_metabolites({production_capacity_constraint: cone_production_cost })
Popsicle_production.add_metabolites({production_capacity_constraint: popsicle_production_cost })

print()
print('Here is what happens in the continuous (LP) case...')

the_program = cobra_model.optimize(objective_sense='maximize')
print()
print(('Status is: %s'%cobra_model.solution.status))
print(('Objective value is: %1.2f'%cobra_model.solution.f))

for the_reaction, the_value in cobra_model.solution.x_dict.items():
    print('%s: %1.2f'%(the_reaction, the_value))

print()
print('Who wants 1/3 of a cone?  Cones and popsicles are units aka integers, reformulate as MILP')
Cone_production.variable_kind = 'integer'
Popsicle_production.variable_kind = 'integer'

the_program = cobra_model.optimize(objective_sense='maximize')
print()
print(('Status is: %s'%cobra_model.solution.status))
print(('Objective value is: %1.2f'%cobra_model.solution.f))

for the_reaction, the_value in cobra_model.solution.x_dict.items():
    print('%s: %1.2f'%(the_reaction, the_value))

print()
print('We now make full items')
