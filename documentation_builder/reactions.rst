
Reactions in cobrapy
====================

This example is available as an IPython
`notebook <http://nbviewer.ipython.org/github/opencobra/cobrapy/blob/master/documentation_builder/reactions.ipynb>`__.

We will consider the reaction glucose 6-phosphate isomerase, which
interconverts glucose 6-phosphate and fructose 6-phosphate. The reaction
id for this reaction in our test model is PGI.

.. code:: python

    import cobra.test
    model = cobra.test.create_test_model()
    pgi = model.reactions.get_by_id("PGI")
    pgi



.. parsed-literal::

    <Reaction PGI at 0x655fc10>



We can view the full name and reaction catalyzed as strings

.. code:: python

    print(pgi.name)
    print(pgi.reaction)

.. parsed-literal::

    1.0 g6p_c <=> 1.0 f6p_c


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

    ['PGI', {'C': -6.0, 'H': -11.0, 'O': -9.0, 'P': -1.0}]



In order to add a metabolite, we pass in a dict with the metabolite
object and its coefficient

.. code:: python

    pgi.add_metabolites({model.metabolites.get_by_id("h_c"): -1})
    pgi.reaction



.. parsed-literal::

    'g6p_c + h_c <=> '



The reaction is no longer mass balanced

.. code:: python

    pgi.check_mass_balance()



.. parsed-literal::

    ['PGI', {'C': 0.0, 'H': -1.0, 'O': 0.0, 'P': 0.0}]



We can remove the metabolite, and the reaction will be balanced once
again.

.. code:: python

    pgi.pop(model.metabolites.get_by_id("h_c"))
    print pgi.reaction
    print pgi.check_mass_balance()

.. parsed-literal::

    1.0 g6p_c <=> 1.0 f6p_c
    []

