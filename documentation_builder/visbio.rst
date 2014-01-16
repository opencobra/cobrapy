Visbio and cobrapy
======================

cobrapy integrates well with the visbio_ package. If it has been
installed, viewing maps is extremely simple. The default map is the *E. coli*
core map, but other maps can be viewed by passing the correct map_name. The
maps are then downloaded from a map repository_.

.. code:: python
>>> from visbio import Map
>>> import cobra.test
>>> model = cobra.test.create_test_model(cobra.test.ecoli_pickle)
>>> model.optimize()
>>> wt_solution = model.solution.x_dict
>>> # mutant flux has PGI knocked out
>>> model.reactions.PGI.lower_bound = 0
>>> model.reactions.PGI.upper_bound = 0
>>> model.optimize()
>>> mutant_solution = model.solution.x_dict
>>> wt_map = Map(flux=wt_solution)
>>> mutant_map = Map(flux=model.solution.x_dict)
>>> wt_map.view_browser()  # opens in a browser
>>> mutant_map.create_standalone_html("mutant.html")  # saves an html file

Using maps is even nicer when using the IPython_ notebook, as shown in this
example_.

.. _visbio: https://github.com/zakandrewking/visbio
.. _repository: https://github.com/zakandrewking/visbio/tree/gh-pages/maps
.. _IPython: http://ipython.org/
.. _example: http://nbviewer.ipython.org/github/opencobra/cobrapy/blob/master/documentation_builder/visbio.ipynb
