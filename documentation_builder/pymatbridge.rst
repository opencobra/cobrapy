
Using the COBRA toolbox with cobrapy
====================================

This example demonstrates using COBRA toolbox commands in MATLAB from
python through
`pymatbridge <http://arokem.github.io/python-matlab-bridge/>`__.

.. code:: python

    %load_ext pymatbridge

.. parsed-literal::

    Starting MATLAB on ZMQ socket ipc:///tmp/pymatbridge
    Send 'exit' command to kill the server
    .....MATLAB started and connected!


.. code:: python

    import cobra.test
    m = cobra.test.create_test_model()
The model\_to\_pymatbridge function will send the model to the workspace
with the given variable name.

.. code:: python

    from cobra.io.mat import model_to_pymatbridge
    model_to_pymatbridge(m, variable_name="model")
Now in the MATLAB workspace, the variable name 'model' holds a COBRA
toolbox struct encoding the model.

.. code:: python

    %%matlab
    model


.. parsed-literal::

    
    model = 
    
                rev: [2546x1 logical]
           metNames: {1802x1 cell}
                  b: [1802x1 double]
                  c: [2546x1 double]
             csense: [1802x1 char]
              genes: {1264x1 cell}
        metFormulas: {1802x1 cell}
               rxns: {2546x1 cell}
            grRules: {2546x1 cell}
           rxnNames: {2546x1 cell}
        description: [28x1 char]
                  S: [1802x2546 double]
                 ub: [2546x1 double]
                 lb: [2546x1 double]
               mets: {1802x1 cell}
         subSystems: {2546x1 cell}
    



First, we have to initialize the COBRA toolbox in MATLAB.

.. code:: python

    %%matlab --silent
    warning('off'); % this works around a pymatbridge bug
    addpath(genpath('~/cobratoolbox/'));
    initCobraToolbox();
Commands from the COBRA toolbox can now be run on the model

.. code:: python

    %%matlab
    optimizeCbModel(model)


.. parsed-literal::

    
    ans = 
    
               x: [2546x1 double]
               f: 0.3800
               y: [1801x1 double]
               w: [2546x1 double]
            stat: 1
        origStat: 5
          solver: 'glpk'
            time: 0.6857
    


