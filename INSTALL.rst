Installation of cobrapy
=======================

For installation help, please use the `Google
Group <http://groups.google.com/group/cobra-pie>`_. For usage
instructions, please see the
`documentation <https://cobrapy.readthedocs.org/en/latest/>`_.

All releases require Python 2.7+ or 3.4+ to be installed before
proceeding. Mac OS X (10.7+) and Ubuntu ship with Python. Windows users
without python can download and install python from the `python
website <https://www.python.org/ftp/python/2.7.9/python-2.7.9.amd64.msi>`_.
Please note that though Anaconda and other python distributions may work
with cobrapy, they are not explicitly supported (yet!).

Stable version installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

cobrapy can be installed with any recent installation of pip.
Instructions for several operating systems are below:

Mac OS X or Linux
-----------------

1. `install
   pip <http://pip.readthedocs.org/en/latest/installing.html>`_.
2. In a terminal, run ``sudo pip install cobra``

We highly recommend updating ``pip`` beforehand (``pip install pip --upgrade``).

Microsoft Windows
-----------------

The preferred installation method on Windows is also to use pip. The
latest Windows installers for Python 2.7 and 3.4 include pip, so if you
use those you will already have pip.

1. In a terminal, run ``C:\Python27\Scripts\pip.exe install cobra`` (you
   may need to adjust the path accordingly).

To install without pip, you will need to download and use the
appropriate installer for your version of python from the `python
package index <https://pypi.python.org/pypi/cobra/>`_.

Installation for development
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Get the `detailed contribution instructions <CONTRIBUTING.rst>`_ for contributing to cobrapy.

Installation of optional dependencies
=====================================

Optional dependencies
~~~~~~~~~~~~~~~~~~~~~

On windows, these can downloaded from [this site]
(http://www.lfd.uci.edu/~gohlke/pythonlibs/). On Mac/Linux, they can be
installed using pip, or from the OS package manager (e.g brew, apt,
yum).

1. `libsbml <http://sbml.org>`_ >= 5.10 to read/write SBML level 2
   files

   -  `Windows libsbml installer <http://www.lfd.uci.edu/~gohlke/pythonlibs/#libsbml>`_
   -  Use ``sudo pip install python-libsbml`` on Mac/Linux

2. `lxml <http://lxml.de/>`_ to speed up read/write of SBML level 3 files.

3. `scipy <http://scipy.org>`_ >= 0.11 for MOMA and saving to \*.mat files.
   -  `Windows scipy installer <http://www.lfd.uci.edu/~gohlke/pythonlibs/#scipy>`_

4. `pytest <http://docs.pytest.org/en/latest/>`_ and `pytest-benchmark
   <http://pytest-benchmark.readthedocs.io/en/latest/>`_ are required
   for testing

You can install all packages directly by

.. code:: shell

    pip install "cobra[all]"

Solvers
~~~~~~~~~~~~~

cobrapy uses `optlang <http://optlang.readthedocs.io>`_ to interface
the mathematical solvers used to optimize the created COBRA models,
which at the time of writing

-  ILOG/CPLEX (available with
   `Academic <https://www.ibm.com/developerworks/university/academicinitiative/>`_
   and
   `Commercial <http://www.ibm.com/software/integration/optimization/cplex-optimizer/>`_
   licenses).
-  `gurobi <http://gurobi.com>`_
-  `glpk <http://www.gnu.org/software/glpk/>`_

Testing your installation
=========================

While it is not a hard requirement for using cobrapy, you need pytest and
pytest-benchmark to run its tests. First do

.. code:: shell

    pip install pytest pytest-benchmark

or to install cobrapy directly with the test dependencies

.. code:: shell

   pip install "cobra[test]"

Then start python and type the following into the Python shell

.. code:: python

    from cobra.test import test_all
    test_all()

You should see some skipped tests and expected failures, and the
function should return ``0``. If you see a value other than ``0`` please file
an `issue report <CONTRIBUTING.rst>`_.
