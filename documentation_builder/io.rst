
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
    
    salmonella_model = cobra.test.create_test_model("salmonella")

.. parsed-literal::

    E. coli test files: 
    ecoli_json, ecoli_mat, ecoli_pickle, ecoli_sbml
    
    Salmonella test files: 
    salmonella_fbc_sbml, salmonella_pickle, salmonella_sbml


JSON
----

cobrapy has a `JSON <https://en.wikipedia.org/wiki/JSON>`__ (JavaScript
Object Notation) representation. This is the ideal format for storing a
cobra model on a computer, or for interoperability with
`escher <https://escher.github.io>`__. Additional fields, however, will
not be saved.

.. code:: python

    cobra.io.load_json_model(cobra.test.ecoli_json)



.. parsed-literal::

    <Model iJO1366 at 0x7f8d2fa20150>



.. code:: python

    cobra.io.write_sbml_model(salmonella_model, "test.json")
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

    <Model Salmonella_consensus_build_1 at 0x7f8d480ad710>



.. code:: python

    cobra.io.read_sbml_model(cobra.test.salmonella_fbc_sbml)



.. parsed-literal::

    <Model Salmonella_consensus_build_1 at 0x7f8d480ad5d0>



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

    <Model iJO1366 at 0x7f8d48090b50>



If the mat file contains only a single model, cobra can figure out which
variable to read from, and the variable\_name paramter is unnecessary.

.. code:: python

    cobra.io.load_matlab_model(cobra.test.ecoli_mat)



.. parsed-literal::

    <Model iJO1366 at 0x7f8d2d85af10>



Saving models to mat files is also relatively straightforward

.. code:: python

    cobra.io.save_matlab_model(ecoli_model, "test_ecoli_model.mat")

::


    ---------------------------------------------------------------------------
    NameError                                 Traceback (most recent call last)

    <ipython-input-9-899c9ca38882> in <module>()
    ----> 1 cobra.io.save_matlab_model(ecoli_model, "test_ecoli_model.mat")
    

    NameError: name 'ecoli_model' is not defined


Pickle
------

Cobra models can be serialized using the python serialization format,
`pickle <https://docs.python.org/2/library/pickle.html>`__. While this
will save any extra fields which may have been created, it does not work
with any other tools and can break between cobrapy major versions. JSON
is generally the preferred format.

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
