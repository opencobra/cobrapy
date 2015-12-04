
Simulating with FBA
===================

Simulations using flux balance analysis can be solved using
Model.optimize(). This will maximize or minimize (maximizing is the
default) flux through the objective reactions.

.. code:: python

    import pandas
    pandas.options.display.max_rows = 100
    
    import cobra.test
    model = cobra.test.create_test_model("textbook")

Running FBA
-----------

.. code:: python

    model.optimize()




.. parsed-literal::

    <Solution 0.87 at 0x7fe558058b50>



The Model.optimize() function will return a Solution object, which will
also be stored at model.solution. A solution object has several
attributes:

-  f: the objective value
-  status: the status from the linear programming solver
-  x\_dict: a dictionary of {reaction\_id: flux\_value} (also called
   "primal")
-  x: a list for x\_dict
-  y\_dict: a dictionary of {metabolite\_id: dual\_value}.
-  y: a list for y\_dict

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

    0.8739215069684305



Changing the Objectives
-----------------------

The objective function is determined from the objective\_coefficient
attribute of the objective reaction(s). Currently in the model, there is
only one objective reaction, with an objective coefficient of 1.

.. code:: python

    model.objective




.. parsed-literal::

    {<Reaction Biomass_Ecoli_core at 0x7fe526516490>: 1.0}



The objective function can be changed by assigning Model.objective,
which can be a reaction object (or just it's name), or a dict of
{Reaction: objective\_coefficient}.

.. code:: python

    # change the objective to ATPM
    # the upper bound should be 1000 so we get the actual optimal value
    model.reactions.get_by_id("ATPM").upper_bound = 1000.
    model.objective = "ATPM"
    model.objective




.. parsed-literal::

    {<Reaction ATPM at 0x7fe526516210>: 1}



.. code:: python

    model.optimize().f




.. parsed-literal::

    174.99999999999997



The objective function can also be changed by setting
Reaction.objective\_coefficient directly.

.. code:: python

    model.reactions.get_by_id("ATPM").objective_coefficient = 0.
    model.reactions.get_by_id("Biomass_Ecoli_core").objective_coefficient = 1.
    model.objective




.. parsed-literal::

    {<Reaction Biomass_Ecoli_core at 0x7fe526516490>: 1.0}



Running FVA
-----------

FBA will not give always give unique solution, because multiple flux
states can achieve the same optimum. FVA (or flux variability analysis)
finds the ranges of each metabolic flux at the optimum.

.. code:: python

    fva_result = cobra.flux_analysis.flux_variability_analysis(model, model.reactions[:20])
    pandas.DataFrame.from_dict(fva_result).T




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>maximum</th>
          <th>minimum</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>ACALD</th>
          <td>9.466331e-29</td>
          <td>3.720797e-15</td>
        </tr>
        <tr>
          <th>ACALDt</th>
          <td>-6.310887e-29</td>
          <td>3.720797e-15</td>
        </tr>
        <tr>
          <th>ACKr</th>
          <td>-2.524355e-28</td>
          <td>3.933509e-15</td>
        </tr>
        <tr>
          <th>ACONTa</th>
          <td>6.007250e+00</td>
          <td>6.007250e+00</td>
        </tr>
        <tr>
          <th>ACONTb</th>
          <td>6.007250e+00</td>
          <td>6.007250e+00</td>
        </tr>
        <tr>
          <th>ACt2r</th>
          <td>6.121561e-28</td>
          <td>3.933509e-15</td>
        </tr>
        <tr>
          <th>ADK1</th>
          <td>-4.042971e-14</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>AKGDH</th>
          <td>5.064376e+00</td>
          <td>5.064376e+00</td>
        </tr>
        <tr>
          <th>AKGt2r</th>
          <td>0.000000e+00</td>
          <td>7.079399e-15</td>
        </tr>
        <tr>
          <th>ALCD2x</th>
          <td>0.000000e+00</td>
          <td>5.729185e-15</td>
        </tr>
        <tr>
          <th>ATPM</th>
          <td>8.390000e+00</td>
          <td>8.390000e+00</td>
        </tr>
        <tr>
          <th>ATPS4r</th>
          <td>4.551401e+01</td>
          <td>4.551401e+01</td>
        </tr>
        <tr>
          <th>Biomass_Ecoli_core</th>
          <td>8.739215e-01</td>
          <td>8.739215e-01</td>
        </tr>
        <tr>
          <th>CO2t</th>
          <td>-2.280983e+01</td>
          <td>-2.280983e+01</td>
        </tr>
        <tr>
          <th>CS</th>
          <td>6.007250e+00</td>
          <td>6.007250e+00</td>
        </tr>
        <tr>
          <th>CYTBD</th>
          <td>4.359899e+01</td>
          <td>4.359899e+01</td>
        </tr>
        <tr>
          <th>D_LACt2</th>
          <td>3.660315e-28</td>
          <td>4.140787e-15</td>
        </tr>
        <tr>
          <th>ENO</th>
          <td>1.471614e+01</td>
          <td>1.471614e+01</td>
        </tr>
        <tr>
          <th>ETOHt2r</th>
          <td>0.000000e+00</td>
          <td>5.729185e-15</td>
        </tr>
        <tr>
          <th>EX_ac_e</th>
          <td>-3.933509e-15</td>
          <td>0.000000e+00</td>
        </tr>
      </tbody>
    </table>
    </div>



Setting parameter fraction\_of\_optimium=0.90 would give the flux ranges
for reactions at 90% optimality.

.. code:: python

    fva_result = cobra.flux_analysis.flux_variability_analysis(model, model.reactions[:20], fraction_of_optimum=0.9)
    pandas.DataFrame.from_dict(fva_result).T




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>maximum</th>
          <th>minimum</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>ACALD</th>
          <td>9.466331e-29</td>
          <td>-2.542370</td>
        </tr>
        <tr>
          <th>ACALDt</th>
          <td>-6.310887e-29</td>
          <td>-2.542370</td>
        </tr>
        <tr>
          <th>ACKr</th>
          <td>-3.029226e-28</td>
          <td>-3.813556</td>
        </tr>
        <tr>
          <th>ACONTa</th>
          <td>8.894520e+00</td>
          <td>0.848587</td>
        </tr>
        <tr>
          <th>ACONTb</th>
          <td>8.894520e+00</td>
          <td>0.848587</td>
        </tr>
        <tr>
          <th>ACt2r</th>
          <td>3.407879e-28</td>
          <td>-3.813556</td>
        </tr>
        <tr>
          <th>ADK1</th>
          <td>1.716100e+01</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>AKGDH</th>
          <td>8.045934e+00</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>AKGt2r</th>
          <td>0.000000e+00</td>
          <td>-1.430083</td>
        </tr>
        <tr>
          <th>ALCD2x</th>
          <td>0.000000e+00</td>
          <td>-2.214323</td>
        </tr>
        <tr>
          <th>ATPM</th>
          <td>2.555100e+01</td>
          <td>8.390000</td>
        </tr>
        <tr>
          <th>ATPS4r</th>
          <td>5.938106e+01</td>
          <td>34.825618</td>
        </tr>
        <tr>
          <th>Biomass_Ecoli_core</th>
          <td>8.739215e-01</td>
          <td>0.786529</td>
        </tr>
        <tr>
          <th>CO2t</th>
          <td>-1.520653e+01</td>
          <td>-26.528850</td>
        </tr>
        <tr>
          <th>CS</th>
          <td>8.894520e+00</td>
          <td>0.848587</td>
        </tr>
        <tr>
          <th>CYTBD</th>
          <td>5.123909e+01</td>
          <td>35.984865</td>
        </tr>
        <tr>
          <th>D_LACt2</th>
          <td>0.000000e+00</td>
          <td>-2.145125</td>
        </tr>
        <tr>
          <th>ENO</th>
          <td>1.673252e+01</td>
          <td>8.686588</td>
        </tr>
        <tr>
          <th>ETOHt2r</th>
          <td>0.000000e+00</td>
          <td>-2.214323</td>
        </tr>
        <tr>
          <th>EX_ac_e</th>
          <td>3.813556e+00</td>
          <td>0.000000</td>
        </tr>
      </tbody>
    </table>
    </div>



Running pFBA
------------

Parsimonious FBA (often written pFBA) finds a flux distribution which
gives the optimal growth rate, but minimizes the total sum of flux. This
involves solving two sequential linear programs, but is handled
transparently by cobrapy. For more details on pFBA, please see `Lewis et
al. (2010) <http://dx.doi.org/10.1038/msb.2010.47>`__.

.. code:: python

    FBA_solution = model.optimize()
    pFBA_solution = cobra.flux_analysis.optimize_minimal_flux(model)

These functions should give approximately the same objective value

.. code:: python

    abs(FBA_solution.f - pFBA_solution.f)




.. parsed-literal::

    1.1102230246251565e-16


