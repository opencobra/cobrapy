
Mixed-Integer Linear Programming
================================

Ice Cream
---------

This example was originally contributed by Joshua Lerman.

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
    popsicle.objective_coefficient = popsicle_selling_price - \
                                     popsicle_production_cost
    
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

Restaurant Order
----------------

To tackle the less immediately obvious problem from the following `XKCD
comic <http://xkcd.com/287/>`__:

.. code:: python

    from IPython.display import Image
    Image(url=r"http://imgs.xkcd.com/comics/np_complete.png")




.. raw:: html

    <img src="http://imgs.xkcd.com/comics/np_complete.png"/>



We want a solution satisfying the following constraints:

:math:`\left(\begin{matrix}2.15&2.75&3.35&3.55&4.20&5.80\end{matrix}\right) \cdot \vec v = 15.05`

:math:`\vec v_i \ge 0`

:math:`\vec v_i \in \mathbb{Z}`

This problem can be written as a COBRA model as well.

.. code:: python

    total_cost = Metabolite("constraint")
    total_cost._bound = 15.05
    
    costs = {"mixed_fruit": 2.15, "french_fries": 2.75, "side_salad": 3.35,
             "hot_wings": 3.55, "mozarella_sticks": 4.20, "sampler_plate": 5.80}
    
    m = Model("appetizers")
    
    for item, cost in costs.items():
        r = Reaction(item)
        r.add_metabolites({total_cost: cost})
        r.variable_kind = "integer"
        m.add_reaction(r)
    
    # To add to the problem, suppose we don't want to eat all mixed fruit.
    m.reactions.mixed_fruit.objective_coefficient = 1
        
    m.optimize(objective_sense="minimize").x_dict




.. parsed-literal::

    {'french_fries': 0.0,
     'hot_wings': 2.0,
     'mixed_fruit': 1.0,
     'mozarella_sticks': 0.0,
     'sampler_plate': 1.0,
     'side_salad': 0.0}



There is another solution to this problem, which would have been
obtained if we had maximized for mixed fruit instead of minimizing.

.. code:: python

    m.optimize(objective_sense="maximize").x_dict




.. parsed-literal::

    {'french_fries': 0.0,
     'hot_wings': 0.0,
     'mixed_fruit': 7.0,
     'mozarella_sticks': 0.0,
     'sampler_plate': 0.0,
     'side_salad': 0.0}



Boolean Indicators
------------------

To give a COBRA-related example, we can create boolean variables as
integers, which can serve as indicators for a reaction being active in a
model. For a reaction flux :math:`v` with lower bound -1000 and upper
bound 1000, we can create a binary variable :math:`b` with the following
constraints:

:math:`b \in \{0, 1\}`

:math:`-1000 \cdot b \le v \le 1000 \cdot b`

To introduce the above constraints into a cobra model, we can rewrite
them as follows

:math:`v \le b \cdot 1000 \Rightarrow v- 1000\cdot b \le 0`

:math:`-1000 \cdot b \le v \Rightarrow v + 1000\cdot b \ge 0`

.. code:: python

    import cobra.test
    model = cobra.test.create_test_model("textbook")
    
    # an indicator for pgi
    pgi = model.reactions.get_by_id("PGI")
    # make a boolean variable
    pgi_indicator = Reaction("indicator_PGI")
    pgi_indicator.lower_bound = 0
    pgi_indicator.upper_bound = 1
    pgi_indicator.variable_kind = "integer"
    # create constraint for v - 1000 b <= 0
    pgi_plus = Metabolite("PGI_plus")
    pgi_plus._constraint_sense = "L"
    # create constraint for v + 1000 b >= 0
    pgi_minus = Metabolite("PGI_minus")
    pgi_minus._constraint_sense = "G"
    
    pgi_indicator.add_metabolites({pgi_plus: -1000, pgi_minus: 1000})
    pgi.add_metabolites({pgi_plus: 1, pgi_minus: 1})
    model.add_reaction(pgi_indicator)
    
    
    # an indicator for zwf
    zwf = model.reactions.get_by_id("G6PDH2r")
    zwf_indicator = Reaction("indicator_ZWF")
    zwf_indicator.lower_bound = 0
    zwf_indicator.upper_bound = 1
    zwf_indicator.variable_kind = "integer"
    # create constraint for v - 1000 b <= 0
    zwf_plus = Metabolite("ZWF_plus")
    zwf_plus._constraint_sense = "L"
    # create constraint for v + 1000 b >= 0
    zwf_minus = Metabolite("ZWF_minus")
    zwf_minus._constraint_sense = "G"
    
    zwf_indicator.add_metabolites({zwf_plus: -1000, zwf_minus: 1000})
    zwf.add_metabolites({zwf_plus: 1, zwf_minus: 1})
    
    # add the indicator reactions to the model
    model.add_reaction(zwf_indicator)


In a model with both these reactions active, the indicators will also be
active

.. code:: python

    solution = model.optimize()
    print("PGI indicator = %d" % solution.x_dict["indicator_PGI"])
    print("ZWF indicator = %d" % solution.x_dict["indicator_ZWF"])
    print("PGI flux = %.2f" % solution.x_dict["PGI"])
    print("ZWF flux = %.2f" % solution.x_dict["G6PDH2r"])


.. parsed-literal::

    PGI indicator = 1
    ZWF indicator = 1
    PGI flux = 4.86
    ZWF flux = 4.96


Because these boolean indicators are in the model, additional
constraints can be applied on them. For example, we can prevent both
reactions from being active at the same time by adding the following
constraint:

:math:`b_\text{pgi} + b_\text{zwf} = 1`

.. code:: python

    or_constraint = Metabolite("or")
    or_constraint._bound = 1
    zwf_indicator.add_metabolites({or_constraint: 1})
    pgi_indicator.add_metabolites({or_constraint: 1})
    
    solution = model.optimize()
    print("PGI indicator = %d" % solution.x_dict["indicator_PGI"])
    print("ZWF indicator = %d" % solution.x_dict["indicator_ZWF"])
    print("PGI flux = %.2f" % solution.x_dict["PGI"])
    print("ZWF flux = %.2f" % solution.x_dict["G6PDH2r"])


.. parsed-literal::

    PGI indicator = 1
    ZWF indicator = 0
    PGI flux = 9.82
    ZWF flux = 0.00

