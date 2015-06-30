
Gapfillling
===========

GrowMatch and SMILEY are gap-filling algorithms, which try to to make
the minimal number of changes to a model and allow it to simulate
growth. For more information, see `Kumar et
al. <http://dx.doi.org/10.1371/journal.pcbi.1000308>`__. Please note
that these algorithms are Mixed-Integer Linear Programs, which need
solvers such as gurobi or cplex to function correctly.

.. code:: python

    import cobra.test
    
    model = cobra.test.create_test_model("salmonella")

In this model D-Fructose-6-phosphate is an essential metabolite. We will
remove all the reactions using it, and at them to a separate model.

.. code:: python

    # remove some reactions and add them to the universal reactions
    Universal = cobra.Model("Universal_Reactions")
    for i in [i.id for i in model.metabolites.f6p_c.reactions]:
        reaction = model.reactions.get_by_id(i)
        Universal.add_reaction(reaction.copy())
        reaction.remove_from_model()

Now, because of these gaps, the model won't grow.

.. code:: python

    model.optimize().f




.. parsed-literal::

    3.067723590211908e-08



We will use GrowMatch to add back the minimal number of reactions from
this set of "universal" reactions (in this case just the ones we
removed) to allow it to grow.

.. code:: python

    cobra.flux_analysis.growMatch(model, Universal)




.. parsed-literal::

    [[<Reaction GF6PTA at 0x7fcecddbc390>,
      <Reaction MAN6PI_reverse at 0x7fcecddbc450>,
      <Reaction F6PA_reverse at 0x7fcecddbc490>,
      <Reaction PGI_reverse at 0x7fcecddbc510>,
      <Reaction TKT2_reverse at 0x7fcecddbc590>]]



We can obtain multiple possible reaction sets by having the algorithm go
through multiple iterations.

.. code:: python

    result = cobra.flux_analysis.growMatch(model, Universal, iterations=4)
    for i, entries in enumerate(result):
        print("---- Run %d ----" % (i + 1))
        for e in entries:
            print(e.id)


.. parsed-literal::

    ---- Run 1 ----
    FBP
    GF6PTA
    MAN6PI_reverse
    PGI_reverse
    TKT2_reverse
    ---- Run 2 ----
    TALA
    F6PP
    GF6PTA
    MAN6PI_reverse
    F6PA_reverse
    ---- Run 3 ----
    F6PP
    GF6PTA
    MAN6PI_reverse
    F6PA_reverse
    TKT2_reverse
    ---- Run 4 ----
    TALA
    FBP
    GF6PTA
    MAN6PI_reverse
    PGI_reverse

