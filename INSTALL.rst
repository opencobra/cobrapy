=======================
Installation of COBRApy
=======================

For installation help, please use the `Google Group
<http://groups.google.com/group/cobra-pie>`_. For usage instructions, please see
the `documentation <https://cobrapy.readthedocs.org/en/latest/>`_.

We only test against Python 3.7+, however, Python 3.4 or higher work mostly. 
For Windows users and possibly also Mac OS users, we recommend using the 
`Anaconda Python <https://www.anaconda.com/>`_ distribution.

Stable version installation
===========================

COBRApy can be installed with any recent installation of pip.  Instructions for
several operating systems are below:

Mac OS X or Linux
-----------------

1. We highly recommend that you create a `Python virtual environment
   <https://realpython.com/python-virtual-environments-a-primer>`_.
2. Install COBRApy when an environment is active by running ``pip install
   cobra`` in the terminal.

Microsoft Windows
-----------------

If you heed our recommendation to use Anaconda, you can open an Anaconda shell
and install COBRApy from the ``conda-forge`` channel.

.. code-block:: console

    conda install -c conda-forge cobra

Installation for development
============================

Get the `detailed contribution instructions <.github/CONTRIBUTING.rst>`_ for
contributing to COBRApy.

Solvers
=======

COBRApy uses `optlang <http://optlang.readthedocs.io>`_ to interface the
mathematical solvers used to optimize the created COBRA models.  At the time of
writing the supported solvers are:

- ILOG/CPLEX (available with `Academic
  <https://www.ibm.com/developerworks/university/academicinitiative/>`_ and
  `Commercial
  <http://www.ibm.com/software/integration/optimization/cplex-optimizer/>`_
  licenses)
- `Gurobi <http://gurobi.com>`_
- `GLPK <http://www.gnu.org/software/glpk/>`_ which is automatically installed
  as swiglpk
