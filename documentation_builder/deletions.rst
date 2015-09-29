
Simulating Deletions
====================

.. code:: python

    import pandas
    from time import time
    
    import cobra.test
    
    cobra_model = cobra.test.create_test_model("textbook")
    ecoli_model = cobra.test.create_test_model("ecoli")

Single Deletions
----------------

Perform all single gene deletions on a model

.. code:: python

    growth_rates, statuses = cobra.flux_analysis.single_gene_deletion(cobra_model)

These can also be done for only a subset of genes

.. code:: python

    growth_rates, statuses = cobra.flux_analysis.single_gene_deletion(cobra_model, cobra_model.genes[:20])
    pandas.DataFrame.from_dict({"growth_rates": growth_rates, "status": statuses})




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>growth_rates</th>
          <th>status</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>b0116</th>
          <td>0.782351</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>b0118</th>
          <td>0.873922</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>b0351</th>
          <td>0.873922</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>b0356</th>
          <td>0.873922</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>b0474</th>
          <td>0.873922</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>b0726</th>
          <td>0.858307</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>b0727</th>
          <td>0.858307</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>b1241</th>
          <td>0.873922</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>b1276</th>
          <td>0.873922</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>b1478</th>
          <td>0.873922</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>b1849</th>
          <td>0.873922</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>b2296</th>
          <td>0.873922</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>b2587</th>
          <td>0.873922</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>b3115</th>
          <td>0.873922</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>b3734</th>
          <td>0.374230</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>b3735</th>
          <td>0.374230</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>b3736</th>
          <td>0.374230</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>b3737</th>
          <td>0.374230</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>b3738</th>
          <td>0.374230</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>s0001</th>
          <td>0.211141</td>
          <td>optimal</td>
        </tr>
      </tbody>
    </table>
    </div>



This can also be done for reactions

.. code:: python

    growth_rates, statuses = cobra.flux_analysis.single_reaction_deletion(cobra_model, cobra_model.reactions[:20])
    pandas.DataFrame.from_dict({"growth_rates": growth_rates, "status": statuses})




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>growth_rates</th>
          <th>status</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>ACALD</th>
          <td>8.739215e-01</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>ACALDt</th>
          <td>8.739215e-01</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>ACKr</th>
          <td>8.739215e-01</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>ACONTa</th>
          <td>-3.963237e-27</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>ACONTb</th>
          <td>6.162976e-33</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>ACt2r</th>
          <td>8.739215e-01</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>ADK1</th>
          <td>8.739215e-01</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>AKGDH</th>
          <td>8.583074e-01</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>AKGt2r</th>
          <td>8.739215e-01</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>ALCD2x</th>
          <td>8.739215e-01</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>ATPM</th>
          <td>9.166475e-01</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>ATPS4r</th>
          <td>3.742299e-01</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>Biomass_Ecoli_core</th>
          <td>0.000000e+00</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>CO2t</th>
          <td>4.616696e-01</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>CS</th>
          <td>-5.916457e-30</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>CYTBD</th>
          <td>2.116629e-01</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>D_LACt2</th>
          <td>8.739215e-01</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>ENO</th>
          <td>-3.266892e-18</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>ETOHt2r</th>
          <td>8.739215e-01</td>
          <td>optimal</td>
        </tr>
        <tr>
          <th>EX_ac_e</th>
          <td>8.739215e-01</td>
          <td>optimal</td>
        </tr>
      </tbody>
    </table>
    </div>



Double Deletions
----------------

Double deletions run in a similar way. Passing in return\_frame=True
will cause them to format the results as a pandas Dataframe

.. code:: python

    cobra.flux_analysis.double_gene_deletion(cobra_model, cobra_model.genes[-10:], return_frame=True)




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>b0724</th>
          <th>b0723</th>
          <th>b0721</th>
          <th>b0729</th>
          <th>b0728</th>
          <th>b2464</th>
          <th>b0008</th>
          <th>b2935</th>
          <th>b2465</th>
          <th>b3919</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>b0724</th>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.704037</td>
        </tr>
        <tr>
          <th>b0723</th>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.704037</td>
        </tr>
        <tr>
          <th>b0721</th>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.704037</td>
        </tr>
        <tr>
          <th>b0729</th>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0.704037</td>
        </tr>
        <tr>
          <th>b0728</th>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0.704037</td>
        </tr>
        <tr>
          <th>b2464</th>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0.873922</td>
          <td>0.864759</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.704037</td>
        </tr>
        <tr>
          <th>b0008</th>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0.864759</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.704037</td>
        </tr>
        <tr>
          <th>b2935</th>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.000000</td>
          <td>0.704037</td>
        </tr>
        <tr>
          <th>b2465</th>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.814298</td>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.000000</td>
          <td>0.873922</td>
          <td>0.704037</td>
        </tr>
        <tr>
          <th>b3919</th>
          <td>0.704037</td>
          <td>0.704037</td>
          <td>0.704037</td>
          <td>0.704037</td>
          <td>0.704037</td>
          <td>0.704037</td>
          <td>0.704037</td>
          <td>0.704037</td>
          <td>0.704037</td>
          <td>0.704037</td>
        </tr>
      </tbody>
    </table>
    </div>



By default, the double deletion function will automatically use
multiprocessing, splitting the task over up to 4 cores if they are
available. The number of cores can be manually sepcified as well.
Setting use of a single core will disable use of the multiprocessing
library, which often aids debuggging.

.. code:: python

    start = time()  # start timer()
    cobra.flux_analysis.double_gene_deletion(ecoli_model, ecoli_model.genes[:100], number_of_processes=2)
    t1 = time() - start
    print("Double gene deletions for 100 genes completed in %.2f sec with 2 cores" % t1)
    
    start = time()  # start timer()
    cobra.flux_analysis.double_gene_deletion(ecoli_model, ecoli_model.genes[:100], number_of_processes=1)
    t2 = time() - start
    print("Double gene deletions for 100 genes completed in %.2f sec with 1 core" % t2)
    
    print("Speedup of %.2fx" % (t2/t1))


.. parsed-literal::

    Double gene deletions for 100 genes completed in 1.69 sec with 2 cores
    Double gene deletions for 100 genes completed in 2.02 sec with 1 core
    Speedup of 1.20x


Double deletions can also be run for reactions

.. code:: python

    cobra.flux_analysis.double_reaction_deletion(cobra_model, cobra_model.reactions[:10], return_frame=True)




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>ACALD</th>
          <th>ACALDt</th>
          <th>ACKr</th>
          <th>ACONTa</th>
          <th>ACONTb</th>
          <th>ACt2r</th>
          <th>ADK1</th>
          <th>AKGDH</th>
          <th>AKGt2r</th>
          <th>ALCD2x</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>ACALD</th>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0</td>
          <td>0</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.858307</td>
          <td>0.873922</td>
          <td>0.873922</td>
        </tr>
        <tr>
          <th>ACALDt</th>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0</td>
          <td>0</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.858307</td>
          <td>0.873922</td>
          <td>0.873922</td>
        </tr>
        <tr>
          <th>ACKr</th>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0</td>
          <td>0</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.858307</td>
          <td>0.873922</td>
          <td>0.873922</td>
        </tr>
        <tr>
          <th>ACONTa</th>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0</td>
          <td>0</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>ACONTb</th>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0</td>
          <td>0</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>ACt2r</th>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0</td>
          <td>0</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.858307</td>
          <td>0.873922</td>
          <td>0.873922</td>
        </tr>
        <tr>
          <th>ADK1</th>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0</td>
          <td>0</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.858307</td>
          <td>0.873922</td>
          <td>0.873922</td>
        </tr>
        <tr>
          <th>AKGDH</th>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0</td>
          <td>0</td>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0.858307</td>
          <td>0.858307</td>
        </tr>
        <tr>
          <th>AKGt2r</th>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0</td>
          <td>0</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.858307</td>
          <td>0.873922</td>
          <td>0.873922</td>
        </tr>
        <tr>
          <th>ALCD2x</th>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0</td>
          <td>0</td>
          <td>0.873922</td>
          <td>0.873922</td>
          <td>0.858307</td>
          <td>0.873922</td>
          <td>0.873922</td>
        </tr>
      </tbody>
    </table>
    </div>


