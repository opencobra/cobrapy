
Mixed-Integer Linear Programming
================================

This example was originally contributed by Joshua Lerman. It is
available as an IPython
`notebook <http://nbviewer.ipython.org/github/opencobra/cobrapy/blob/master/documentation_builder/milp.ipynb>`__.

An ice cream stand sells cones and popsicles. It wants to maximize its
profit, but is subject to a budget.

We can write this problem as a linear program:

    **max** cone :math:`\cdot` cone\_margin + popsicle :math:`\cdot`
    popsicle margin

    *subject to*

    cone :math:`\cdot` cone\_cost + popsicle :math:`\cdot`
    popsicle\_cost :math:`\le` budget

.. code:: python

    cone_selling_price = 7.
    cone_production_cost = 3.
    popsicle_selling_price = 2.
    popsicle_production_cost = 1.
    starting_budget = 100.

This problem can be written as a cobra.Model

.. code:: python

    from cobra import Model, Metabolite, Reaction
    
    cone = Reaction("cone")
    popsicle = Reaction("popsicle")
    
    # constrainted to a budget
    budget = Metabolite("budget")
    budget._constraint_sense = "L"
    budget._bound = starting_budget
    cone.add_metabolites({budget: cone_production_cost})
    popsicle.add_metabolites({budget: popsicle_production_cost})
    
    # objective coefficient is the profit to be made from each unit
    cone.objective_coefficient = cone_selling_price - cone_production_cost
    popsicle.objective_coefficient = popsicle_selling_price - popsicle_production_cost
    
    m = Model("lerman_ice_cream_co")
    m.add_reactions((cone, popsicle))
    
    m.optimize().x_dict



.. parsed-literal::

    {'cone': 33.333333333333336, 'popsicle': 0.0}



In reality, cones and popsicles can only be sold in integer amounts. We
can use the variable kind attribute of a cobra.Reaction to enforce this.

.. code:: python

    cone.variable_kind = "integer"
    popsicle.variable_kind = "integer"
    m.optimize().x_dict



.. parsed-literal::

    {'cone': 33.0, 'popsicle': 1.0}



Now the model makes both popsicles and cones.
