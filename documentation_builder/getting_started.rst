
Getting Started
===============

To begin with, cobrapy comes with bundled models for *Salmonella* and
*E. coli*, as well as a "textbook" model of *E. coli* core metabolism.
To load a test model, type

.. code:: python

    from __future__ import print_function
    import cobra.test
    
    model = cobra.test.create_test_model("textbook")  # "ecoli" and "salmonella" are also valid arguments

The reactions, metabolites, and genes attributes of the cobrapy model
are a special type of list called a DictList, and each one is made up of
Reaction, Metabolite and Gene objects respectively.

.. code:: python

    print(len(model.reactions))
    print(len(model.metabolites))
    print(len(model.genes))


.. parsed-literal::

    95
    72
    137


Just like a regular list, objects in the DictList can be retrived by
index. For example, to get the 30th reaction in the model (at index 29
because of
`0-indexing <https://en.wikipedia.org/wiki/Z%20ero-based_numbering>`__):

.. code:: python

    model.reactions[29]




.. parsed-literal::

    <Reaction EX_glu__L_e at 0x7fbbe05e5590>



Addictionally, items can be retrived by their id using the get\_by\_id()
function. For example, to get the cytosolic atp metabolite object (the
id is "atp\_c"), we can do the following:

.. code:: python

    model.metabolites.get_by_id("atp_c")




.. parsed-literal::

    <Metabolite atp_c at 0x7fbbe0617350>



As an added bonus, users with an interactive shell such as IPython will
be able to tab-complete to list elements inside a list. While this is
not recommended behavior for most code because of the possibility for
characters like "-" inside ids, this is very useful while in an
interactive prompt:

.. code:: python

    model.reactions.EX_glc__D_e.lower_bound




.. parsed-literal::

    -10.0



Reactions
---------

We will consider the reaction glucose 6-phosphate isomerase, which
interconverts glucose 6-phosphate and fructose 6-phosphate. The reaction
id for this reaction in our test model is PGI.

.. code:: python

    pgi = model.reactions.get_by_id("PGI")
    pgi




.. parsed-literal::

    <Reaction PGI at 0x7fbbe0611790>



We can view the full name and reaction catalyzed as strings

.. code:: python

    print(pgi.name)
    print(pgi.reaction)


.. parsed-literal::

    glucose-6-phosphate isomerase
    g6p_c <=> f6p_c


We can also view reaction upper and lower bounds. Because the
pgi.lower\_bound < 0, and pgi.upper\_bound > 0, pgi is reversible

.. code:: python

    print(pgi.lower_bound, "< pgi <", pgi.upper_bound)
    print(pgi.reversibility)


.. parsed-literal::

    -1000.0 < pgi < 1000.0
    True


We can also ensure the reaction is mass balanced. This function will
return elements which violate mass balance. If it comes back empty, then
the reaction is mass balanced.

.. code:: python

    pgi.check_mass_balance()




.. parsed-literal::

    {}



In order to add a metabolite, we pass in a dict with the metabolite
object and its coefficient

.. code:: python

    pgi.add_metabolites({model.metabolites.get_by_id("h_c"): -1})
    pgi.reaction




.. parsed-literal::

    'g6p_c + h_c <=> f6p_c'



The reaction is no longer mass balanced

.. code:: python

    pgi.check_mass_balance()




.. parsed-literal::

    {'H': -1.0}



We can remove the metabolite, and the reaction will be balanced once
again.

.. code:: python

    pgi.pop(model.metabolites.get_by_id("h_c"))
    print(pgi.reaction)
    print(pgi.check_mass_balance())


.. parsed-literal::

    g6p_c <=> f6p_c
    {}


It is also possible to build the reaction from a string. However, care
must be taken when doing this to ensure reaction id's match those in the
model. The direction of the arrow is also used to update the upper and
lower bounds.

.. code:: python

    pgi.reaction = "g6p_c --> f6p_c + h_c + green_eggs + ham"


.. parsed-literal::

    unknown metabolite 'green_eggs' created
    unknown metabolite 'ham' created


.. code:: python

    pgi.reaction




.. parsed-literal::

    'g6p_c --> green_eggs + ham + h_c + f6p_c'



.. code:: python

    pgi.reaction = "g6p_c <=> f6p_c"
    pgi.reaction




.. parsed-literal::

    'g6p_c <=> f6p_c'



Metabolites
-----------

We will consider cytosolic atp as our metabolite, which has the id
atp\_c in our test model.

.. code:: python

    atp = model.metabolites.get_by_id("atp_c")
    atp




.. parsed-literal::

    <Metabolite atp_c at 0x7fbbe0617350>



We can print out the metabolite name and compartment (cytosol in this
case).

.. code:: python

    print(atp.name)
    print(atp.compartment)


.. parsed-literal::

    ATP
    c


We can see that ATP is a charged molecule in our model.

.. code:: python

    atp.charge




.. parsed-literal::

    -4



We can see the chemical formula for the metabolite as well.

.. code:: python

    print(atp.formula)


.. parsed-literal::

    C10H12N5O13P3


The reactions attribute gives a frozenset of all reactions using the
given metabolite. We can use this to count the number of reactions which
use atp.

.. code:: python

    len(atp.reactions)




.. parsed-literal::

    13



A metabolite like glucose 6-phosphate will participate in fewer
reactions.

.. code:: python

    model.metabolites.get_by_id("g6p_c").reactions




.. parsed-literal::

    frozenset({<Reaction G6PDH2r at 0x7fbbe05fd050>,
               <Reaction GLCpts at 0x7fbbe05fd150>,
               <Reaction PGI at 0x7fbbe0611790>,
               <Reaction Biomass_Ecoli_core at 0x7fbbe0650ed0>})



Genes
-----

The gene\_reaction\_rule is a boolean representation of the gene
requirements for this reaction to be active as described in
`Schellenberger et al 2011 Nature Protocols
6(9):1290-307 <http://dx.doi.org/doi:10.1038/nprot.2011.308>`__.

The GPR is stored as the gene\_reaction\_rule for a Reaction object as a
string.

.. code:: python

    gpr = pgi.gene_reaction_rule
    gpr




.. parsed-literal::

    'b4025'



Corresponding gene objects also exist. These objects are tracked by the
reactions itself, as well as by the model

.. code:: python

    pgi.genes




.. parsed-literal::

    frozenset({<Gene b4025 at 0x7fbbe063dc90>})



.. code:: python

    pgi_gene = model.genes.get_by_id("b4025")
    pgi_gene




.. parsed-literal::

    <Gene b4025 at 0x7fbbe063dc90>



Each gene keeps track of the reactions it catalyzes

.. code:: python

    pgi_gene.reactions




.. parsed-literal::

    frozenset({<Reaction PGI at 0x7fbbe0611790>})



Altering the gene\_reaction\_rule will create new gene objects if
necessary and update all relationships.

.. code:: python

    pgi.gene_reaction_rule = "(spam or eggs)"
    pgi.genes




.. parsed-literal::

    frozenset({<Gene eggs at 0x7fbbe0611b50>, <Gene spam at 0x7fbbe0611e90>})



.. code:: python

    pgi_gene.reactions




.. parsed-literal::

    frozenset()



Newly created genes are also added to the model

.. code:: python

    model.genes.get_by_id("spam")




.. parsed-literal::

    <Gene spam at 0x7fbbe0611e90>



The delete\_model\_genes function will evaluate the gpr and set the
upper and lower bounds to 0 if the reaction is knocked out. This
function can preserve existing deletions or reset them using the
cumulative\_deletions flag.

.. code:: python

    cobra.manipulation.delete_model_genes(model, ["spam"], cumulative_deletions=True)
    print(pgi.lower_bound, "< pgi <", pgi.upper_bound)
    cobra.manipulation.delete_model_genes(model, ["eggs"], cumulative_deletions=True)
    print(pgi.lower_bound, "< pgi <", pgi.upper_bound)


.. parsed-literal::

    -1000 < pgi < 1000
    0.0 < pgi < 0.0


The undelete\_model\_genes can be used to reset a gene deletion

.. code:: python

    cobra.manipulation.undelete_model_genes(model)
    print(pgi.lower_bound, "< pgi <", pgi.upper_bound)


.. parsed-literal::

    -1000 < pgi < 1000

