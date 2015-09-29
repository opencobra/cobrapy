
FAQ
===

This document will address frequently asked questions not addressed in
other pages of the documentation.

How do I install cobrapy?
~~~~~~~~~~~~~~~~~~~~~~~~~

Please see the
`INSTALL.md <https://github.com/opencobra/cobrapy/blob/master/INSTALL.md>`__
file.

How do I cite cobrapy?
~~~~~~~~~~~~~~~~~~~~~~

Please cite the 2013 publication:
`10.1186/1752-0509-7-74 <http://dx.doi.org/doi:10.1186/1752-0509-7-74>`__

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
        print(repr(e))


.. parsed-literal::

    KeyError('test_dcaACP_c',)


The Model.repair function will rebuild the necessary indexes

.. code:: python

    model.repair()
    model.metabolites.get_by_id(model.metabolites[0].id)




.. parsed-literal::

    <Metabolite test_dcaACP_c at 0x688b450>



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

How do I change the reversibility of a Reaction?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Reaction.reversibility is a property in cobra which is computed when it
is requested from the lower and upper bounds.

.. code:: python

    model = cobra.test.create_test_model()
    model.reactions.get_by_id("PGI").reversibility




.. parsed-literal::

    True



Trying to set it directly will result in an error:

.. code:: python

    try:
        model.reactions.get_by_id("PGI").reversibility = False
    except Exception as e:
        print(repr(e))


.. parsed-literal::

    AttributeError("can't set attribute",)


The way to change the reversibility is to change the bounds to make the
reaction irreversible.

.. code:: python

    model.reactions.get_by_id("PGI").lower_bound = 10
    model.reactions.get_by_id("PGI").reversibility




.. parsed-literal::

    False



How do I generate an LP file from a COBRA model?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While the cobrapy does not include python code to support this feature
directly, many of the bundled solvers have this capability. Create the
problem with one of these solvers, and use its appropriate function.

Please note that unlike the LP file format, the MPS file format does not
specify objective direction and is always a minimzation. Some (but not
all) solvers will rewrite the maximization as a minimzation.

.. code:: python

    model = cobra.test.create_test_model()
    # glpk through cglpk
    glp = cobra.solvers.cglpk.create_problem(model)
    glp.write("test.lp")
    glp.write("test.mps")  # will not rewrite objective
    # gurobi
    gurobi_problem = cobra.solvers.gurobi_solver.create_problem(model)
    gurobi_problem.write("test.lp")
    gurobi_problem.write("test.mps")  # rewrites objective
    # cplex
    cplex_problem = cobra.solvers.cplex_solver.create_problem(model)
    cplex_problem.write("test.lp")
    cplex_problem.write("test.mps")  # rewrites objective

How do I visualize my flux solutions?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cobrapy works well with the `escher <https://escher.github.io/>`__
package, which is well suited to this purpose. Consult the `escher
documentation <https://escher.readthedocs.org/en/latest/>`__ for
examples.
