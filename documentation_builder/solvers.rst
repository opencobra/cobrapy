
Solver Interface
================

Each cobrapy solver must expose the following API. The solvers all will
have their own distinct LP object types, but each can be manipulated by
these functions. This API can be used directly when implementing
algorithms efficiently on linear programs because it has 2 primary
benefits:

1. Avoid the overhead of creating and destroying LP's for each operation

2. Many solver objects preserve the basis between subsequent LP's,
   making each subsequent LP solve faster

We will walk though the API with the cglpk solver, which links the
cobrapy solver API with `GLPK <http://www.gnu.org/software/glpk/>`__'s C
API.

.. code:: python

    import cobra.test
    
    model = cobra.test.create_test_model("textbook")
    solver = cobra.solvers.cglpk

Attributes and functions
------------------------

Each solver has some attributes:

solver\_name
~~~~~~~~~~~~

The name of the solver. This is the name which will be used to select
the solver in cobrapy functions.

.. code:: python

    solver.solver_name




.. parsed-literal::

    'cglpk'



.. code:: python

    model.optimize(solver="cglpk")




.. parsed-literal::

    <Solution 0.87 at 0x7f9148bb2250>



\_SUPPORTS\_MILP
~~~~~~~~~~~~~~~~

The presence of this attribute tells cobrapy that the solver supports
mixed-integer linear programming

.. code:: python

    solver._SUPPORTS_MILP




.. parsed-literal::

    True



solve
~~~~~

Model.optimize is a wrapper for each solver's solve function. It takes
in a cobra model and returns a solution

.. code:: python

    solver.solve(model)




.. parsed-literal::

    <Solution 0.87 at 0x7f917d50ed50>



create\_problem
~~~~~~~~~~~~~~~

This creates the LP object for the solver.

.. code:: python

    lp = solver.create_problem(model, objective_sense="maximize")
    lp




.. parsed-literal::

    <cobra.solvers.cglpk.GLP at 0x46e8aa0>



solve\_problem
~~~~~~~~~~~~~~

Solve the LP object and return the solution status

.. code:: python

    solver.solve_problem(lp)




.. parsed-literal::

    'optimal'



format\_solution
~~~~~~~~~~~~~~~~

Extract a cobra.Solution object from a solved LP object

.. code:: python

    solver.format_solution(lp, model)




.. parsed-literal::

    <Solution 0.87 at 0x7f917d50e9d0>



get\_objective\_value
~~~~~~~~~~~~~~~~~~~~~

Extract the objective value from a solved LP object

.. code:: python

    solver.get_objective_value(lp)




.. parsed-literal::

    0.8739215069684305



get\_status
~~~~~~~~~~~

Get the solution status of a solved LP object

.. code:: python

    solver.get_status(lp)




.. parsed-literal::

    'optimal'



change\_variable\_objective
~~~~~~~~~~~~~~~~~~~~~~~~~~~

change the objective coefficient a reaction at a particular index. This
does not change any of the other objectives which have already been set.
This example will double and then revert the biomass coefficient.

.. code:: python

    model.reactions.index("Biomass_Ecoli_core")




.. parsed-literal::

    12



.. code:: python

    solver.change_variable_objective(lp, 12, 2)
    solver.solve_problem(lp)
    solver.get_objective_value(lp)




.. parsed-literal::

    1.747843013936861



.. code:: python

    solver.change_variable_objective(lp, 12, 1)
    solver.solve_problem(lp)
    solver.get_objective_value(lp)




.. parsed-literal::

    0.8739215069684305



change variable\_bounds
~~~~~~~~~~~~~~~~~~~~~~~

change the lower and upper bounds of a reaction at a particular index.
This example will set the lower bound of the biomass to an infeasible
value, then revert it.

.. code:: python

    solver.change_variable_bounds(lp, 12, 1000, 1000)
    solver.solve_problem(lp)




.. parsed-literal::

    'infeasible'



.. code:: python

    solver.change_variable_bounds(lp, 12, 0, 1000)
    solver.solve_problem(lp)




.. parsed-literal::

    'optimal'



change\_coefficient
~~~~~~~~~~~~~~~~~~~

Change a coefficient in the stoichiometric matrix. In this example, we
will set the entry for ADP in the ATMP reaction to in infeasible value,
then reset it.

.. code:: python

    model.metabolites.index("atp_c")




.. parsed-literal::

    16



.. code:: python

    model.reactions.index("ATPM")




.. parsed-literal::

    10



.. code:: python

    solver.change_coefficient(lp, 16, 10, -10)
    solver.solve_problem(lp)




.. parsed-literal::

    'infeasible'



.. code:: python

    solver.change_coefficient(lp, 16, 10, -1)
    solver.solve_problem(lp)




.. parsed-literal::

    'optimal'



set\_parameter
~~~~~~~~~~~~~~

Set a solver parameter. Each solver will have its own particular set of
unique paramters. However, some have unified names. For example, all
solvers should accept "tolerance\_feasibility."

.. code:: python

    solver.set_parameter(lp, "tolerance_feasibility", 1e-9)

.. code:: python

    solver.set_parameter(lp, "objective_sense", "minimize")
    solver.solve_problem(lp)
    solver.get_objective_value(lp)




.. parsed-literal::

    0.0



.. code:: python

    solver.set_parameter(lp, "objective_sense", "maximize")
    solver.solve_problem(lp)
    solver.get_objective_value(lp)




.. parsed-literal::

    0.8739215069684304



Example with FVA
----------------

Consider flux variability analysis (FVA), which requires maximizing and
minimizing every reaction with the original biomass value fixed at its
optimal value. If we used the cobra Model API in a naive implementation,
we would do the following:

.. code:: python

    %%time
    # work on a copy of the model so the original is not changed
    fva_model = model.copy()
    
    # set the lower bound on the objective to be the optimal value
    f = fva_model.optimize().f
    for objective_reaction, coefficient in fva_model.objective.items():
        objective_reaction.lower_bound = coefficient * f
    
    # now maximize and minimze every reaction to find its bounds
    fva_result = {}
    for r in fva_model.reactions:
        fva_model.change_objective(r)
        fva_result[r.id] = {}
        fva_result[r.id]["maximum"] = fva_model.optimize(objective_sense="maximize").f
        fva_result[r.id]["minimum"] = fva_model.optimize(objective_sense="minimize").f


.. parsed-literal::

    CPU times: user 144 ms, sys: 667 µs, total: 145 ms
    Wall time: 141 ms


Instead, we could use the solver API to do this more efficiently. This
is roughly how cobrapy implementes FVA. It keeps uses the same LP object
and repeatedly maximizes and minimizes it. This allows the solver to
preserve the basis, and is much faster. The speed increase is even more
noticeable the larger the model gets.

.. code:: python

    %%time
    # create the LP object
    lp = solver.create_problem(model)
    
    # set the lower bound on the objective to be the optimal value
    solver.solve_problem(lp)
    f = solver.get_objective_value(lp)
    for objective_reaction, coefficient in model.objective.items():
        objective_index = model.reactions.index(objective_reaction)
        # old objective is no longer the objective
        solver.change_variable_objective(lp, objective_index, 0.)
        solver.change_variable_bounds(lp, objective_index, f * coefficient, objective_reaction.upper_bound)
    
    # now maximize and minimze every reaction to find its bounds
    fva_result = {}
    for index, r in enumerate(model.reactions):
        solver.change_variable_objective(lp, index, 1.)
        fva_result[r.id] = {}
        solver.solve_problem(lp, objective_sense="maximize")
        fva_result[r.id]["maximum"] = solver.get_objective_value(lp)
        solver.solve_problem(lp, objective_sense="minimize")
        fva_result[r.id]["minimum"] = solver.get_objective_value(lp)
        solver.change_variable_objective(lp, index, 0.)


.. parsed-literal::

    CPU times: user 9.85 ms, sys: 251 µs, total: 10.1 ms
    Wall time: 9.94 ms

