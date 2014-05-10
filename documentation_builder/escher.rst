Escher and cobrapy
==================

cobrapy integrates well with the escher_ package. If it has been
installed, viewing maps is extremely simple. Maps can be
easily downloaded and viewed from a map repository_.

.. code:: python

    from escher import Builder
    import cobra.test
    
    model = cobra.test.create_test_model(cobra.test.ecoli_pickle)
    cobra.flux_analysis.parsimonious.optimize_minimal_flux(model)
    wt_solution = model.solution.x_dict
    
    # mutant flux has PGI knocked out
    model.reactions.PGI.lower_bound = 0
    model.reactions.PGI.upper_bound = 0
    cobra.flux_analysis.parsimonious.optimize_minimal_flux(model)
    mutant_solution = model.solution.x_dict

    wt_map = Builder("iJO1366_central_metabolism", reaction_data=wt_solution)
    wt_map.display_in_browser()


    mutant_map = Builder("iJO1366_central_metabolism", reaction_data=mutant_solution)
    mutant_map.display_in_browser()

Using maps is even nicer when using the IPython_ notebook, as shown in this
example_.

.. _escher: http://zakandrewking.github.io/escher/
.. _repository: https://github.com/zakandrewking/escher/tree/gh-pages/maps
.. _IPython: http://ipython.org/
.. _example: http://nbviewer.ipython.org/github/opencobra/cobrapy/blob/master/documentation_builder/escher.ipynb
