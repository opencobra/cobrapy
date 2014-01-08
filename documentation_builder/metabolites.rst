
Metabolites in cobrapy
======================

We will consider cytosolic atp as our metabolite, which has the id
atp\_c in our test model.

.. code:: python

>>> import cobra.test
>>> model = cobra.test.create_test_model()
>>> atp = model.metabolites.get_by_id("atp_c")
>>> print atp.name
    ATP


The metabolite is labeled as being in the cytosol

.. code:: python

>>> atp.compartment
    'c'



We can see that ATP is a charged molecule in our model.

.. code:: python

>>> atp.charge
    -4



We can see the chemical formula for the metabolite as well.

.. code:: python

>>> atp.formula
    'C10H12N5O13P3'
