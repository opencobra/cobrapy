
Reading and Writing Models
==========================

This example is available as an IPython
`notebook <http://nbviewer.ipython.org/github/opencobra/cobrapy/blob/master/documentation_builder/io.ipynb>`__.

Functions for reading and writing models to various formats are included
with cobrapy. The package also ships with models of *E. coli* and
*Salmonella* in various formats for testing purposes. In this example,
we will use these functions to read models from these test files in
various formats.

.. code:: python

    import cobra.test
    
    print("E. coli test files: ")
    print(", ".join([i for i in dir(cobra.test) if i.startswith("ecoli")]))
    print("")
    print("Salmonella test files: ")
    print(", ".join([i for i in dir(cobra.test) if i.startswith("salmonella")]))

.. parsed-literal::

    E. coli test files: 
    ecoli_json, ecoli_mat, ecoli_pickle, ecoli_sbml
    
    Salmonella test files: 
    salmonella_fbc_sbml, salmonella_pickle, salmonella_reaction_p_values_pickle, salmonella_sbml


Pickle
------

Cobra models can be serialized using the python serialization format,
`pickle <https://docs.python.org/2/library/pickle.html>`__.

.. code:: python

    from cPickle import load, dump
    
    # read in the test models
    with open(cobra.test.ecoli_pickle, "rb") as infile:
        ecoli_model = load(infile)
    with open(cobra.test.salmonella_pickle, "rb") as infile:
        salmonella_model = load(infile)
    
    # output to a file
    with open("test_output.pickle", "wb") as outfile:
        dump(salmonella_model, outfile)

SBML
----

The `Systems Biology Markup Language <http://sbml.org>`__ is an
XML-based standard format for distributing models. Cobrapy can use
`libsbml <http://sbml.org/Software/libSBML>`__, which must be installed
separately (see installation instructions) to read and write SBML files.

Initially, the COBRA format for SBML files used the "notes" field in
SBML files. More recently, however, the `FBC
extension <http://sbml.org/Documents/Specifications/SBML_Level_3/Packages/Flux_Balance_Constraints_%28flux%29>`__
to SBML has come into existence, which defines its own fields.

Cobrapy can handle both formats (assuming libsbml has been installed
correctly). When reading in a model, it will automatically detect
whether fbc was used or not. When writing a model, the use\_fbc\_package
can be used.

.. code:: python

    cobra.io.read_sbml_model(cobra.test.salmonella_sbml)



.. parsed-literal::

    <Model Salmonella_consensus_build_1 at 0x3e24d90>



.. code:: python

    cobra.io.read_sbml_model(cobra.test.salmonella_fbc_sbml)



.. parsed-literal::

    <Model Salmonella_consensus_build_1 at 0xf9dda50>



.. code:: python

    cobra.io.write_sbml_model(salmonella_model, "test.xml",
                              use_fbc_package=False)
    cobra.io.write_sbml_model(salmonella_model, "test_fbc.xml",
                              use_fbc_package=True)

MATLAB
------

Often, models may be imported and exported soley for the purposes of
working with the same models in cobrapy and the `MATLAB cobra
toolbox <http://opencobra.github.io/cobratoolbox/>`__. MATLAB has its
own ".mat" format for storing variables. Reading and writing to these
mat files from python requires scipy, and is generally much faster than
using libsbml.

A mat file can contain multiple MATLAB variables. Therefore, the
variable name of the model in the MATLAB file can be passed into the
reading function:

.. code:: python

    cobra.io.load_matlab_model(cobra.test.ecoli_mat, variable_name="iJO1366")



.. parsed-literal::

    <Model iJO1366 at 0xf9dde90>



If the mat file contains only a single model, cobra can figure out which
variable to read from, and the variable\_name paramter is unnecessary.

.. code:: python

    cobra.io.load_matlab_model(cobra.test.ecoli_mat)



.. parsed-literal::

    <Model iJO1366 at 0xf9ddf90>



Saving models to mat files is also relatively straightforward

.. code:: python

    cobra.io.save_matlab_model(ecoli_model, "output_path.mat")