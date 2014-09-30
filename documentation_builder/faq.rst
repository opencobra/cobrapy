
FAQ
===

This document will address frequently asked questions not addressed in
other pages of the documentation.

How do I install cobrapy?
~~~~~~~~~~~~~~~~~~~~~~~~~

Please see the
`INSTALL.md <https://github.com/opencobra/cobrapy/blob/master/INSTALL.md>`__
file.

How do I rename reactions or metabolites?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TL;DR Use Model.repair afterwards

When renaming metabolites or reactions, there are issues because cobra
indexes based off of ID's, which can cause errors. For example:

.. code:: python

    from __future__ import print_function
    import cobra.test
    model = cobra.test.create_test_model()
    
    for metabolite in model.metabolites:
        metabolite.id = "test_" + metabolite.id
    
    try:
        model.metabolites.get_by_id(model.metabolites[0].id)
    except KeyError as e:
        print("KeyError:", e)

.. parsed-literal::

    KeyError: 'test_dcaACP_c'


The Model.repair function will rebuild the necessary indexes

.. code:: python

    model.repair()
    model.metabolites.get_by_id(model.metabolites[0].id)



.. parsed-literal::

    <Metabolite test_dcaACP_c at 0x5f68ed0>



How do I delete a gene?
~~~~~~~~~~~~~~~~~~~~~~~

That depends on what precisely you mean by delete a gene.

If you want to simulate the model with a gene knockout, use the
cobra.maniupulation.delete\_model\_genes function. The effects of this
function are reversed by cobra.manipulation.undelete\_model\_genes.

.. code:: python

    model = cobra.test.create_test_model()
    PGI = model.reactions.get_by_id("PGI")
    print("bounds before knockout:", (PGI.lower_bound, PGI.upper_bound))
    cobra.manipulation.delete_model_genes(model, ["STM4221"])
    print("bounds after knockouts", (PGI.lower_bound, PGI.upper_bound))

.. parsed-literal::

    bounds before knockout: (-1000.0, 1000.0)
    bounds after knockouts (0.0, 0.0)


If you want to actually remove all traces of a gene from a model, this
is more difficult because this will require changing all the
gene\_reaction\_rule strings for reactions involving the gene.
