
Reading and Writing Models
==========================

Cobrapy supports reading and writing models in SBML (with and without
FBC), JSON, MAT, and pickle formats. Generally, SBML with FBC is the
preferred format for general use.

The package also ships with test models in various formats for testing
purposes.

.. code:: python

    import cobra.test
    import os
    
    print("mini test files: ")
    print(", ".join([i for i in os.listdir(cobra.test.data_directory) if i.startswith("mini")]))
    
    textbook_model = cobra.test.create_test_model("textbook")
    ecoli_model = cobra.test.create_test_model("ecoli")
    salmonella_model = cobra.test.create_test_model("salmonella")


.. parsed-literal::

    mini test files: 
    mini.mat, mini_cobra.xml, mini.json, mini_fbc2.xml, mini_fbc1.xml, mini.pickle


SBML
----

The `Systems Biology Markup Language <http://sbml.org>`__ is an
XML-based standard format for distributing models which has support for
COBRA models through the `FBC
extension <http://sbml.org/Documents/Specifications/SBML_Level_3/Packages/Flux_Balance_Constraints_%28flux%29>`__
version 2.

Cobrapy has native support for reading and writing SBML with FBCv2.
Please note that all id's in the model must conform to the SBML SID
requirements in order to generate a valid SBML file.

.. code:: python

    cobra.io.read_sbml_model(os.path.join(cobra.test.data_directory, "mini_fbc2.xml"))




.. parsed-literal::

    <Model mini_textbook at 0x7f82fc4aec90>



.. code:: python

    cobra.io.write_sbml_model(textbook_model, "test_fbc2.xml")

There are other dialects of SBML prior to FBC 2 which have previously
been use to encode COBRA models. The primary ones is the "COBRA" dialect
which used the "notes" fields in SBML files.

Cobrapy can use `libsbml <http://sbml.org/Software/libSBML>`__, which
must be installed separately (see installation instructions) to read and
write these files. When reading in a model, it will automatically detect
whether fbc was used or not. When writing a model, the use\_fbc\_package
flag can be used can be used.

.. code:: python

    cobra.io.read_sbml_model(os.path.join(cobra.test.data_directory, "mini_cobra.xml"))




.. parsed-literal::

    <Model mini_textbook at 0x7f82cc744c10>



.. code:: python

    cobra.io.write_sbml_model(textbook_model, "test_cobra.xml", use_fbc_package=False)

JSON
----

cobrapy models have a `JSON <https://en.wikipedia.org/wiki/JSON>`__
(JavaScript Object Notation) representation. This format was crated for
interoperability with `escher <https://escher.github.io>`__. Additional
fields, however, will not be saved.

.. code:: python

    cobra.io.load_json_model(os.path.join(cobra.test.data_directory, "mini.json"))




.. parsed-literal::

    <Model mini_textbook at 0x7f82fc4d4f50>



.. code:: python

    cobra.io.write_sbml_model(textbook_model, "test.json")

MATLAB
------

Often, models may be imported and exported soley for the purposes of
working with the same models in cobrapy and the `MATLAB cobra
toolbox <http://opencobra.github.io/cobratoolbox/>`__. MATLAB has its
own ".mat" format for storing variables. Reading and writing to these
mat files from python requires scipy.

A mat file can contain multiple MATLAB variables. Therefore, the
variable name of the model in the MATLAB file can be passed into the
reading function:

.. code:: python

    cobra.io.load_matlab_model(os.path.join(cobra.test.data_directory, "mini.mat"),
                               variable_name="mini_textbook")




.. parsed-literal::

    <Model mini_textbook at 0x7f82cc97b110>



If the mat file contains only a single model, cobra can figure out which
variable to read from, and the variable\_name paramter is unnecessary.

.. code:: python

    cobra.io.load_matlab_model(os.path.join(cobra.test.data_directory, "mini.mat"))




.. parsed-literal::

    <Model mini_textbook at 0x7f82fc4d4450>



Saving models to mat files is also relatively straightforward

.. code:: python

    cobra.io.save_matlab_model(textbook_model, "test.mat")

Pickle
------

Cobra models can be serialized using the python serialization format,
`pickle <https://docs.python.org/2/library/pickle.html>`__. While this
will save any extra fields which may have been created, it does not work
with any other tools and can break between cobrapy major versions. JSON,
SBML, and MAT are generally the preferred format.

.. code:: python

    from pickle import load, dump
    
    # read in the test models
    with open(os.path.join(cobra.test.data_directory, "mini.pickle"), "rb") as infile:
        mini_model = load(infile)
    
    # output to a file
    with open("test.pickle", "wb") as outfile:
        dump(textbook_model, outfile)
