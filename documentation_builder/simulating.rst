
Simulating with FBA
===================

This example is available as an IPython
`notebook <http://nbviewer.ipython.org/github/opencobra/cobrapy/blob/master/documentation_builder/simulating.ipynb>`__.

Simulations using flux balance analysis can be solved using
Model.optimize(). This will maximize or minimize (maximizing is the
default) flux through the objective reactions.

.. code:: python

    import cobra.test
    model = cobra.test.create_test_model()

Running FBA
-----------

.. code:: python

    model.optimize()



.. parsed-literal::

    <Solution 0.38 at 0x660d990>



The Model.optimize() function will return a Solution object, which will
also be stored at model.solution. A solution object has several
attributes: - f: the objective value - status: the status from the
linear programming solver - x\_dict: a dictionary of {reaction\_id:
flux\_value} (also called "primal") - x: a list for x\_dict - y\_dict: a
dictionary of {metabolite\_id: dual\_value}. - y: a list for y\_dict

For example, after the last call to model.optimize(), the status should
be 'optimal' if the solver returned no errors, and f should be the
objective value

.. code:: python

    model.solution.status



.. parsed-literal::

    'optimal'



.. code:: python

    model.solution.f



.. parsed-literal::

    0.38000797227551136



Changing the Objectives
-----------------------

The objective function is determined from the objective\_coefficient
attribute of the objective reaction(s). Currently in the model, there is
only one objective reaction, with an objective coefficient of 1.

.. code:: python

    {reaction: reaction.objective_coefficient for reaction in model.reactions
            if reaction.objective_coefficient > 0}



.. parsed-literal::

    {<Reaction biomass_iRR1083_metals at 0x660d350>: 1.0}



The objective function can be changed by using the function
Model.change\_objective, which will take either a reaction object or
just its name.

.. code:: python

    # change the objective to ATPM
    # the upper bound should be 1000 so we get the actual optimal value
    model.reactions.get_by_id("ATPM").upper_bound = 1000.
    model.change_objective("ATPM")
    {reaction: reaction.objective_coefficient for reaction in model.reactions
            if reaction.objective_coefficient > 0}



.. parsed-literal::

    {<Reaction ATPM at 0x52cb190>: 1.0}



.. code:: python

    model.optimize()



.. parsed-literal::

    <Solution 119.67 at 0x4c93110>



The objective function can also be changed by setting
Reaction.objective\_coefficient directly.

.. code:: python

    model.reactions.get_by_id("ATPM").objective_coefficient = 0.
    model.reactions.get_by_id("biomass_iRR1083_metals").objective_coefficient = 1.
    {reaction: reaction.objective_coefficient for reaction in model.reactions
            if reaction.objective_coefficient > 0}



.. parsed-literal::

    {<Reaction biomass_iRR1083_metals at 0x660d350>: 1.0}


