
Building a Model
================

This simple example demonstrates how to create a model, create a
reaction, and then add the reaction to the model.

We'll use the '3OAS140' reaction from the STM\_1.0 model:

1.0 malACP[c] + 1.0 h[c] + 1.0 ddcaACP[c] :math:`\rightarrow` 1.0 co2[c]
+ 1.0 ACP[c] + 1.0 3omrsACP[c]

First, create the model and reaction.

.. code:: python

    from cobra import Model, Reaction, Metabolite
    # Best practise: SBML compliant IDs
    cobra_model = Model('example_cobra_model')
    
    reaction = Reaction('3OAS140')
    reaction.name = '3 oxoacyl acyl carrier protein synthase n C140 '
    reaction.subsystem = 'Cell Envelope Biosynthesis'
    reaction.lower_bound = 0.  # This is the default
    reaction.upper_bound = 1000.  # This is the default
    reaction.objective_coefficient = 0. # this is the default

We need to create metabolites as well. If we were using an existing
model, we could use get\_by\_id to get the apporpriate Metabolite
objects instead.

.. code:: python

    ACP_c = Metabolite('ACP_c',
                       formula='C11H21N2O7PRS',
                       name='acyl-carrier-protein',
                       compartment='c')
    omrsACP_c = Metabolite('3omrsACP_c',
                           formula='C25H45N2O9PRS',
                           name='3-Oxotetradecanoyl-acyl-carrier-protein',
                           compartment='c')
    co2_c = Metabolite('co2_c',
                       formula='CO2',
                       name='CO2',
                       compartment='c')
    malACP_c = Metabolite('malACP_c',
                          formula='C14H22N2O10PRS',
                          name='Malonyl-acyl-carrier-protein',
                          compartment='c')
    h_c = Metabolite('h_c',
                     formula='H',
                     name='H',
                     compartment='c')
    ddcaACP_c = Metabolite('ddcaACP_c',
                           formula='C23H43N2O8PRS',
                           name='Dodecanoyl-ACP-n-C120ACP',
                           compartment='c')

Adding metabolites to a reaction requires using a dictionary of the
metabolites and their stoichiometric coefficients. A group of
metabolites can be added all at once, or they can be added one at a
time.

.. code:: python

    reaction.add_metabolites({malACP_c: -1.0,
                              h_c: -1.0,
                              ddcaACP_c: -1.0,
                              co2_c: 1.0,
                              ACP_c: 1.0,
                              omrsACP_c: 1.0})
    
    
    reaction.reaction  # This gives a string representation of the reaction




.. parsed-literal::

    'malACP_c + h_c + ddcaACP_c --> 3omrsACP_c + ACP_c + co2_c'



The gene\_reaction\_rule is a boolean representation of the gene
requirements for this reaction to be active as described in
`Schellenberger et al 2011 Nature Protocols
6(9):1290-307 <http://dx.doi.org/doi:10.1038/nprot.2011.308>`__. We will
assign the gene reaction rule string, which will automatically create
the corresponding gene objects.

.. code:: python

    reaction.gene_reaction_rule = '( STM2378 or STM1197 )'
    reaction.genes




.. parsed-literal::

    frozenset({<Gene STM1197 at 0x7feea0ae9850>, <Gene STM2378 at 0x7feea0ae9b10>})



At this point in time, the model is still empty

.. code:: python

    print('%i reactions in initial model' % len(cobra_model.reactions))
    print('%i metabolites in initial model' % len(cobra_model.metabolites))
    print('%i genes in initial model' % len(cobra_model.genes))


.. parsed-literal::

    0 reactions in initial model
    0 metabolites in initial model
    0 genes in initial model


We will add the reaction to the model, which will also add all
associated metabolites and genes

.. code:: python

    cobra_model.add_reaction(reaction)
    
    # Now there are things in the model
    print('%i reaction in model' % len(cobra_model.reactions))
    print('%i metabolites in model' % len(cobra_model.metabolites))
    print('%i genes in model' % len(cobra_model.genes))


.. parsed-literal::

    1 reaction in model
    6 metabolites in model
    2 genes in model


We can iterate through the model objects to observe the contents

.. code:: python

    # Iterate through the the objects in the model
    print("Reactions")
    print("---------")
    for x in cobra_model.reactions:
        print("%s : %s" % (x.id, x.reaction))
    print("Metabolites")
    print("-----------")
    for x in cobra_model.metabolites:
        print('%s : %s' % (x.id, x.formula))
    print("Genes")
    print("-----")
    for x in cobra_model.genes:
        reactions_list_str = "{" + ", ".join((i.id for i in x.reactions)) + "}"
        print("%s is associated with reactions: %s" % (x.id, reactions_list_str))


.. parsed-literal::

    Reactions
    ---------
    3OAS140 : malACP_c + h_c + ddcaACP_c --> 3omrsACP_c + ACP_c + co2_c
    Metabolites
    -----------
    3omrsACP_c : C25H45N2O9PRS
    ACP_c : C11H21N2O7PRS
    co2_c : CO2
    malACP_c : C14H22N2O10PRS
    h_c : H
    ddcaACP_c : C23H43N2O8PRS
    Genes
    -----
    STM2378 is associated with reactions: {3OAS140}
    STM1197 is associated with reactions: {3OAS140}

