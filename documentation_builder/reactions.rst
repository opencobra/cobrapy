
Reactions in cobrapy
====================

We will consider the reaction glucose 6-phosphate isomerase, which
interconverts glucose 6-phosphate and fructose 6-phosphate. The reaction
id for this reaction in our test model is PGI.

.. code:: python

>>> import cobra.test
>>> model = cobra.test.create_test_model()
>>> pgi = model.reactions.get_by_id("PGI")
>>> print pgi.name
    glucose 6 phosphate isomerase


We can view the reaction catalyzed by PGI as a reaction string

.. code:: python

>>> print pgi.reaction
    1.0 g6p_c <=> 1.0 f6p_c


We can also view reaction upper and lower bounds. Because the
pgi.lower\_bound < 0, and pgi.upper\_bound > 0, pgi is reversible

.. code:: python

>>> print pgi.lower_bound, "< pgi <", pgi.upper_bound
    -1000.0 < pgi < 1000.0
>>> print pgi.reversibility
    True


We can also ensure the reaction is mass balanced

.. code:: python

>>> pgi.check_mass_balance()
    []



In order to add a metabolite, we pass in a dict with the metabolite
object and its coefficient

.. code:: python

>>> pgi.add_metabolites({model.metabolites.get_by_id("h_c"): -1})
>>> print pgi.reaction
    1.0 g6p_c + 1 h_c <=> 1.0 f6p_c


The reaction is no longer mass balanced

.. code:: python

>>> pgi.check_mass_balance()
    ['PGI', {'C': 0.0, 'H': -1.0, 'O': 0.0, 'P': 0.0}]



We can remove the metabolite, and the reaction will be balanced once
again.

.. code:: python

>>> pgi.pop(model.metabolites.get_by_id("h_c"))
>>> print pgi.reaction
    1.0 g6p_c <=> 1.0 f6p_c
>>> print pgi.check_mass_balance()
    []

View the IPython notebook_.

.. _notebook: http://nbviewer.ipython.org/github/opencobra/cobrapy/blob/master/documentation_builder/reactions.ipynb
