
## Building a Model

# This simple example (available as an IPython [notebook](http://nbviewer.ipyth
# on.org/github/opencobra/cobrapy/blob/master/documentation_builder/building_mo
# del.ipynb)) demonstrates how to create a model, create a reaction, and then
# add the reaction to the model.
# 
# We'll use the '3OAS140' reaction from the STM_1.0 model:
# 
# 1.0 malACP[c] + 1.0 h[c] + 1.0 ddcaACP[c] $\rightarrow$ 1.0 co2[c] + 1.0
# ACP[c] + 1.0 3omrsACP[c]

# First, create the model and reaction.

from cobra import Model, Reaction, Metabolite
# Best practise: SBML compliant IDs
cobra_model = Model('example_cobra_model')

reaction = Reaction('3OAS140')
reaction.name = '3 oxoacyl acyl carrier protein synthase n C140 '
reaction.subsystem = 'Cell Envelope Biosynthesis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.objective_coefficient = 0. # this is the default


# We need to create metabolites as well. If we were using an existing model, we
# could use get_by_id to get the apporpriate Metabolite objects instead.

ACP_c = Metabolite('ACP_c', formula='C11H21N2O7PRS',
    name='acyl-carrier-protein', compartment='c')
omrsACP_c = Metabolite('3omrsACP_c', formula='C25H45N2O9PRS',
    name='3-Oxotetradecanoyl-acyl-carrier-protein', compartment='c')
co2_c = Metabolite('co2_c', formula='CO2', name='CO2', compartment='c')
malACP_c = Metabolite('malACP_c', formula='C14H22N2O10PRS',
    name='Malonyl-acyl-carrier-protein', compartment='c')
h_c = Metabolite('h_c', formula='H', name='H', compartment='c')
ddcaACP_c = Metabolite('ddcaACP_c', formula='C23H43N2O8PRS',
    name='Dodecanoyl-ACP-n-C120ACP', compartment='c')


# Adding metabolites to a reaction requires using a dictionary of the
# metabolites and their stoichiometric coefficients. A group of metabolites can
# be added all at once, or they can be added one at a time.

reaction.add_metabolites({malACP_c: -1.0,
                          h_c: -1.0,
                          ddcaACP_c: -1.0,
                          co2_c: 1.0,
                          ACP_c: 1.0,
                          omrsACP_c: 1.0})


reaction.reaction  # This gives a string representation of the reaction
# Output:
# 'malACP_c + h_c + ddcaACP_c --> co2_c + 3omrsACP_c + ACP_c'

# The gene_reaction_rule is a boolean representation of the gene requirements
# for this reaction to be active as described in [Schellenberger et al 2011
# Nature Protocols
# 6(9):1290-307](http://dx.doi.org/doi:10.1038/nprot.2011.308). We will assign
# the gene reaction rule string, which will automatically create the
# corresponding gene objects.

reaction.gene_reaction_rule = '( STM2378  or STM1197 )'
reaction.genes
# Output:
# frozenset({<Gene STM2378 at 0x3739b10>, <Gene STM1197 at 0x3739b50>})

# At this point in time, the model is still empty

print('%i reactions in initial model' % len(cobra_model.reactions))
print('%i metabolites in initial model' % len(cobra_model.metabolites))
print('%i genes in initial model' % len(cobra_model.genes))
# Prints:
# 0 reactions in initial model
# 0 metabolites in initial model
# 0 genes in initial model

# We will add the reaction to the model, which will also add all associated
# metabolites and genes

cobra_model.add_reaction(reaction)

# Now there are things in the model
print('%i reaction in model' % len(cobra_model.reactions))
print('%i metabolites in model' % len(cobra_model.metabolites))
print('%i genes in model' % len(cobra_model.genes))
# Prints:
# 1 reaction in model
# 6 metabolites in model
# 2 genes in model

# We can iterate through the model objects to observe the contents

# Iterate through the the objects in the model
print("Reactions")
print("---------")
for x in cobra_model.reactions:
    print("%s : %s" % (repr(x), x.reaction))
print("Metabolites")
print("-----------")
for x in cobra_model.metabolites:
    print('%s : %s' % (repr(x), x.formula))
print("Genes")
print("-----")
for x in cobra_model.genes:
    reactions_list_str = ", ".join((repr(i) for i in x.reactions))
    print("%s is associated with reactions: %s" % (repr(x), reactions_list_str))
# Prints:
# Reactions
# ---------
# <Reaction 3OAS140 at 0x4a18b90> : malACP_c + h_c + ddcaACP_c --> co2_c + 3omrsACP_c + ACP_c
# Metabolites
# -----------
# <Metabolite co2_c at 0x594ba10> : CO2
# <Metabolite malACP_c at 0x594ba90> : C14H22N2O10PRS
# <Metabolite h_c at 0x594bb10> : H
# <Metabolite 3omrsACP_c at 0x594b950> : C25H45N2O9PRS
# <Metabolite ACP_c at 0x594b990> : C11H21N2O7PRS
# <Metabolite ddcaACP_c at 0x594bb90> : C23H43N2O8PRS
# Genes
# -----
# <Gene STM2378 at 0x3739b10> is associated with reactions: <Reaction 3OAS140 at 0x4a18b90>
# <Gene STM1197 at 0x3739b50> is associated with reactions: <Reaction 3OAS140 at 0x4a18b90>
