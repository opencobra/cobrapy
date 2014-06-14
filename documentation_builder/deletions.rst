
Simulating Deletions
====================

This example is available as an IPython
`notebook <http://nbviewer.ipython.org/github/opencobra/cobrapy/blob/master/documentation_builder/deletions.ipynb>`__.

.. code:: python

    from time import time
    
    
    from cobra.test import create_test_model, salmonella_pickle, ecoli_pickle
    from cobra.flux_analysis import single_deletion
    from cobra.flux_analysis import double_deletion
    
    
    cobra_model = create_test_model(salmonella_pickle)

Single Deletions
----------------

                Perform all single gene deletions on a model
                
.. code:: python

    start = time()  # start timer()
    growth_rates, statuses = single_deletion(cobra_model)
    print("All single gene deletions completed in %.2f seconds" % (time() - start))

.. parsed-literal::

    All single gene deletions completed in 4.15 seconds


These can also be done for only a subset of genes

.. code:: python

    single_deletion(cobra_model, element_list=cobra_model.genes[:100]);

Single deletions can also be run on reactions

.. code:: python

    start = time()  # start timer()
    growth_rates, statuses = single_deletion(cobra_model, element_type="reaction")
    print("All single reaction deletions completed in %.2f seconds" % (time() - start))

.. parsed-literal::

    All single reaction deletions completed in 7.62 seconds


Double Deletions
----------------

Double deletions run in a similar way

.. code:: python

    start = time()  # start timer()
    double_deletion(cobra_model, element_list_1=cobra_model.genes[:100])
    print("Double gene deletions for 100 genes completed in %.2f seconds" % (time() - start))

.. parsed-literal::

    Double gene deletions for 100 genes completed in 8.57 seconds


By default, the double deletion function will use multiprocessing,
splitting the task over up to 4 cores if they are available. The number
of cores can be manually sepcified as well. Setting use of a single core
will disable use of the multiprocessing library, which often aids
debuggging.

.. code:: python

    start = time()  # start timer()
    double_deletion(cobra_model, element_list_1=cobra_model.genes[:100], number_of_processes=4)
    print("Double gene deletions for 100 genes completed in %.2f seconds with 4 cores" % (time() - start))
    
    start = time()  # start timer()
    double_deletion(cobra_model, element_list_1=cobra_model.genes[:100], number_of_processes=1)
    print("Double gene deletions for 100 genes completed in %.2f seconds with 1 core" % (time() - start))

.. parsed-literal::

    Double gene deletions for 100 genes completed in 8.31 seconds with 4 cores
    Double gene deletions for 100 genes completed in 21.53 seconds with 1 core


Double deletions can also be run for reactions

.. code:: python

    start = time()
    double_deletion(cobra_model, element_list_1=cobra_model.reactions[:100], element_type="reaction")
    print("Double reaction deletions for 100 reactions completed in %.2f seconds" % (time() - start))

.. parsed-literal::

    Double reaction deletions for 100 reactions completed in 7.26 seconds


If pandas is installed, the results can be returned formatted as a
pandas.DataFrame

.. code:: python

    frame = double_deletion(cobra_model, element_list_1=cobra_model.reactions[300:308], element_type="reaction", return_frame=True)
    frame[frame < 1e-9] = 0.  # round small values to 0
    frame



.. raw:: html

    <div style="max-height:1000px;max-width:1500px;overflow:auto;">
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>ARBtex</th>
          <th>ARGAGMt7pp</th>
          <th>ARGDC</th>
          <th>ARGDCpp</th>
          <th>ARGORNt7pp</th>
          <th>ARGSL</th>
          <th>ARGSS</th>
          <th>ARGTRS</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>ARBtex</th>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0</td>
          <td> 0</td>
          <td> 0.380008</td>
        </tr>
        <tr>
          <th>ARGAGMt7pp</th>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0</td>
          <td> 0</td>
          <td> 0.380008</td>
        </tr>
        <tr>
          <th>ARGDC</th>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0</td>
          <td> 0</td>
          <td> 0.380008</td>
        </tr>
        <tr>
          <th>ARGDCpp</th>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0</td>
          <td> 0</td>
          <td> 0.380008</td>
        </tr>
        <tr>
          <th>ARGORNt7pp</th>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0</td>
          <td> 0</td>
          <td> 0.380008</td>
        </tr>
        <tr>
          <th>ARGSL</th>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0</td>
          <td> 0</td>
          <td> 0.000000</td>
        </tr>
        <tr>
          <th>ARGSS</th>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0.000000</td>
          <td> 0</td>
          <td> 0</td>
          <td> 0.000000</td>
        </tr>
        <tr>
          <th>ARGTRS</th>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0.380008</td>
          <td> 0</td>
          <td> 0</td>
          <td> 0.380008</td>
        </tr>
      </tbody>
    </table>
    </div>


