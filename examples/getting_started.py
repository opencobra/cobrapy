
## Getting Started

# This example is available as an IPython [notebook](http://nbviewer.ipython.or
# g/github/opencobra/cobrapy/blob/master/documentation_builder/getting_started.
# ipynb).

# To begin with, cobrapy comes with two bundled models for _Salmonella_ and _E.
# coli_. To load a test model, type

from __future__ import print_function
import cobra.test
model = cobra.test.create_test_model(cobra.test.salmonella_pickle)


# The reactions, metabolites, and genes attributes of the cobrapy model are are
# a special type of list called a DictList, and each one is made up of
# Reaction, Metabolite and Gene objects respectively.

print(len(model.reactions))
print(len(model.metabolites))
print(len(model.genes))
# Prints:
# 2546
# 1802
# 1264

# Just like a regular list, objects in the DictList can be retrived by index.
# For example, to get the 30th reaction in the model (at index 29 because of
# [0-indexing](https://en.wikipedia.org/wiki/Z ero-based_numbering)):

model.reactions[29]
# Output:
# <Reaction 2AGPA180tipp at 0x6613f90>

# Addictionally, items can be retrived by their id using the get_by_id()
# function. For example, to get the cytosolic atp metabolite object (the id is
# “atp_c”), we can do the following:

model.metabolites.get_by_id("atp_c")
# Output:
# <Metabolite atp_c at 0x65a03d0>

# As an added bonus, users with an interactive shell such as IPython will be
# able to tab-complete to list elements inside a list. While this is not
# recommended behavior for most code because of the possibility for characters
# like "-" inside ids, this is very useful while in an interactive prompt:

model.reactions.EX_glc__D_e.lower_bound
# Output:
# 0.0

### Reactions

# We will consider the reaction glucose 6-phosphate isomerase, which
# interconverts glucose 6-phosphate and fructose 6-phosphate. The reaction id
# for this reaction in our test model is PGI.

pgi = model.reactions.get_by_id("PGI")
pgi
# Output:
# <Reaction PGI at 0x783e350>

# We can view the full name and reaction catalyzed as strings

print(pgi.name)
print(pgi.reaction)
# Prints:
# glucose 6 phosphate isomerase
# g6p_c <=> f6p_c

# We can also view reaction upper and lower bounds. Because the pgi.lower_bound
# < 0, and pgi.upper_bound > 0, pgi is reversible

print(pgi.lower_bound, "< pgi <", pgi.upper_bound)
print(pgi.reversibility)
# Prints:
# -1000.0 < pgi < 1000.0
# True

# We can also ensure the reaction is mass balanced. This function will return
# elements which violate mass balance. If it comes back empty, then the
# reaction is mass balanced.

pgi.check_mass_balance()
# Output:
# []

# In order to add a metabolite, we pass in a dict with the metabolite object
# and its coefficient

pgi.add_metabolites({model.metabolites.get_by_id("h_c"): -1})
pgi.reaction
# Output:
# 'g6p_c + h_c <=> f6p_c'

# The reaction is no longer mass balanced

pgi.check_mass_balance()
# Output:
# ['PGI', {'C': 0.0, 'H': -1.0, 'O': 0.0, 'P': 0.0}]

# We can remove the metabolite, and the reaction will be balanced once again.

pgi.pop(model.metabolites.get_by_id("h_c"))
print(pgi.reaction)
print(pgi.check_mass_balance())
# Prints:
# g6p_c <=> f6p_c
# []

### Metabolites

# We will consider cytosolic atp as our metabolite, which has the id atp_c in
# our test model.

atp = model.metabolites.get_by_id("atp_c")
atp
# Output:
# <Metabolite atp_c at 0x65a03d0>

# We can print out the metabolite name and compartment (cytosol in this case).

print(atp.name)
print(atp.compartment)
# Prints:
# ATP
# c

# We can see that ATP is a charged molecule in our model.

atp.charge
# Output:
# -4

# We can see the chemical formula for the metabolite as well.

print(atp.formula)
# Prints:
# C10H12N5O13P3

# The reactions attribute gives a frozenset of all reactions using the given
# metabolite. We can use this to count the number of reactions which use atp.

len(atp.reactions)
# Output:
# 348

# A metabolite like glucose 6-phosphate will participate in fewer reactions.

model.metabolites.get_by_id("g6p_c").reactions
# Output:
# frozenset({<Reaction G6PDH2r at 0x72c99d0>,
#            <Reaction G6PP at 0x72c9bd0>,
#            <Reaction G6Pt6_2pp at 0x72c9c90>,
#            <Reaction GLCptspp at 0x72e1890>,
#            <Reaction HEX1 at 0x74a9c10>,
#            <Reaction PGI at 0x783e350>,
#            <Reaction PGMT at 0x783e8d0>,
#            <Reaction TRE6PH at 0x7bc7290>,
#            <Reaction TRE6PS at 0x7bc7550>,
#            <Reaction AB6PGH at 0x7f79bd0>})

### Genes

# The gene_reaction_rule is a boolean representation of the gene requirements
# for this reaction to be active as described in [Schellenberger et al 2011
# Nature Protocols
# 6(9):1290-307](http://dx.doi.org/doi:10.1038/nprot.2011.308).
# 
# The GPR is stored as the gene_reaction_rule for a Reaction object as a
# string.

gpr = pgi.gene_reaction_rule
gpr
# Output:
# 'STM4221'

# Corresponding gene objects also exist. These objects are tracked by the
# reactions itself, as well as by the model

pgi.genes
# Output:
# frozenset({<Gene STM4221 at 0x783e3d0>})

pgi_gene = model.genes.get_by_id("STM4221")
pgi_gene
# Output:
# <Gene STM4221 at 0x783e3d0>

# Each gene keeps track of the reactions it catalyzes

pgi_gene.reactions
# Output:
# frozenset({<Reaction PGI at 0x783e350>})

# Altering the gene_reaction_rule will create new gene objects if necessary and
# update all relationships.

pgi.gene_reaction_rule = "(spam or eggs)"
pgi.genes
# Output:
# frozenset({<Gene spam at 0x7f80b10>, <Gene eggs at 0x7f80cd0>})

pgi_gene.reactions
# Output:
# frozenset()

# Newly created genes are also added to the model

model.genes.get_by_id("spam")
# Output:
# <Gene spam at 0x7f80b10>

# The delete_model_genes function will evaluate the gpr and set the upper and
# lower bounds to 0 if the reaction is knocked out. This function can preserve
# existing deletions or reset them using the cumulative_deletions flag.

cobra.manipulation.delete_model_genes(model, ["spam"], cumulative_deletions=True)
print(pgi.lower_bound, "< pgi <", pgi.upper_bound)
cobra.manipulation.delete_model_genes(model, ["eggs"], cumulative_deletions=True)
print(pgi.lower_bound, "< pgi <", pgi.upper_bound)
# Prints:
# -1000.0 < pgi < 1000.0
# 0.0 < pgi < 0.0

# The undelete_model_genes can be used to reset a gene deletion

cobra.manipulation.undelete_model_genes(model)
print(pgi.lower_bound, "< pgi <", pgi.upper_bound)
# Prints:
# -1000.0 < pgi < 1000.0
