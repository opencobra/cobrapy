# This example shows how to modify an existing cobra model and change the
# objective function. We'll use the '3 oxoacyl acyl carrier protein synthase n C140'
# reaction from the STM_1.0 model which currently has the ID of 'my_new_reaction':
#  
# 1.0 malACP[c] + 1.0 h[c] + 1.0 ddcaACP[c] -> 1.0 co2[c] + 1.0 ACP[c] + 1.0 3omrsACP[c]
#

from time import time

from cobra.flux_analysis import single_deletion
from cobra import Reaction
from cobra.manipulation import initialize_growth_medium
from cobra.test import create_test_model, salmonella_pickle  # test filename

cobra_model = create_test_model(salmonella_pickle)

reaction = Reaction('my_new_reaction')
reaction.name = '3 oxoacyl acyl carrier protein synthase n C140'
reaction.subsystem = 'Cell Envelope Biosynthesis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.objective_coefficient = 0.  # This is the default

# Adding metabolites to a reaction requires using a dictionary of the 
# metabolites and their stoichiometric coefficients.  A group of metabolites 
# can be added all at once or they can be added one at a time.

# Create the metabolites
from cobra import Metabolite
ACP_c = Metabolite('ACP_c', formula='C11H21N2O7PRS',
    name='acyl-carrier-protein', compartment='c')
omrsACP_c = Metabolite('3omrsACP_c', formula='C25H45N2O9PRS',
    name='3-Oxotetradecanoyl-acyl-carrier-protein', compartment='c')
co2_c = Metabolite('co2_c', formula='CO2', name='CO2', compartment='c')
malACP_c = Metabolite('malACP_c', formula='C14H22N2O10PRS',
    name='Malonyl-acyl-carrier-protein', compartment='c')
h_c = Metabolite('h_c', formula='H',
    name='H', compartment='c')
# If the metabolites are already in Model.metabolites then it is permissible,
# but not necessary, to use them instead of creating new ones
ddcaACP_c = cobra_model.metabolites.get_by_id('ddcaACP_c')

#add the metabolites to the reaction:
reaction.add_metabolites({malACP_c: -1.0,
                          h_c: -1.0,
                          ddcaACP_c: -1.0,
                          co2_c: 1.0,
                          ACP_c: 1.0,
                          omrsACP_c: 1.0})

#Print the reaction to make sure it worked:
print((reaction.reaction))

print('%i reactions in original model' % len(cobra_model.reactions))
#Add the reaction to the model
cobra_model.add_reaction(reaction)
print(('%i reaction in updated model' % len(cobra_model.reactions)))

cobra_model.optimize()
print('solution for old objective (should be approximately 0.380):')
print((cobra_model.solution.f))

print(('Changing objective to newly added reaction: %s' % reaction.id))
cobra_model.change_objective(reaction)
cobra_model.optimize()
print(('solution for %s (should be approximately 0.845):' % reaction.id))
print((cobra_model.solution.f))
