
Metabolites in cobrapy
======================

This example is available as an IPython
`notebook <http://nbviewer.ipython.org/github/opencobra/cobrapy/blob/master/documentation_builder/metabolites.ipynb>`__.

We will consider cytosolic atp as our metabolite, which has the id
atp\_c in our test model.

.. code:: python

    import cobra.test
    model = cobra.test.create_test_model()
    atp = model.metabolites.get_by_id("atp_c")
    atp



.. parsed-literal::

    <Metabolite atp_c at 0x5b91c90>



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

    348



A metabolite like glucose 6-phosphate will participate in fewer
reactions.

.. code:: python

    model.metabolites.get_by_id("g6p_c").reactions



.. parsed-literal::

    frozenset({<Reaction G6PDH2r at 0x6841090>,
               <Reaction G6PP at 0x6841290>,
               <Reaction G6Pt6_2pp at 0x6841350>,
               <Reaction GLCptspp at 0x6853f90>,
               <Reaction HEX1 at 0x6a18390>,
               <Reaction PGI at 0x6f50ad0>,
               <Reaction PGMT at 0x6f55090>,
               <Reaction TRE6PH at 0x7129a90>,
               <Reaction TRE6PS at 0x7129d50>,
               <Reaction AB6PGH at 0x74f0450>})


