cobrapy - Constraint-Based Reconstruction and Analysis in Python
================================================================

|Build Status| |Coverage Status| |Build status| |PyPI| |Gitter|

What is cobrapy?
~~~~~~~~~~~~~~~~

COBRA methods are widely used for genome-scale modeling of metabolic
networks in both prokaryotes and eukaryotes. **cobrapy** is a
constraint-based modeling package that is designed to accommodate the
biological complexity of the next generation of COBRA models and
provides access to commonly used COBRA methods, such as flux balance
analysis, flux variability analysis, and gene deletion analyses.

Our aim with cobrapy is to provide useful, efficient infrastructure
for:

- creating and managing metabolic models
- accessing popular solvers
- analyzing models with methods such as FVA, FBA, pFBA, MOMA etc.
- inspecting models and drawing conclusions on gene essentiality,
  testing consequences of knock-outs etc.

Our goal with cobrapy is for it to be useful on its own, and for it to
be the natural choice of infrastructure for developers that want to
build new COBRA related python packages for e.g. visualization,
strain-design and data driven analysis. By re-using the same classes
and design principles, we can make new methods both easier to
implement and easier to use, thereby bringing the power of COBRA to
more researchers.

The documentation is browseable online at
`readthedocs <https://cobrapy.readthedocs.org/en/stable/>`_ and can
also be
`downloaded <https://readthedocs.org/projects/cobrapy/downloads/>`_.

Please use the `Google
Group <http://groups.google.com/group/cobra-pie>`_ for help. By writing a well formulated question, with sufficient
detail, you are much more likely to quickly receive a good answer! Please refer to these `StackOverflow
guidelines <https://stackoverflow.com/help/how-to-ask>`_ on how to ask questions.
Alternatively, you can use `gitter.im <https://gitter.im/opencobra/cobrapy>`_ for quick questions
and discussions about cobrapy (faster response times). Please keep in mind that answers are provided on a volunteer
basis.

More information about opencobra is available at the
`website <http://opencobra.github.io/>`_.

If you use cobrapy in a scientific publication, please cite
`doi:10.1186/1752-0509-7-74 <http://dx.doi.org/doi:10.1186/1752-0509-7-74>`_

Installation
~~~~~~~~~~~~

Use pip to `install cobrapy from
PyPI <https://pypi.python.org/pypi/cobra>`_ (we recommend doing this
inside a `virtual
environment <http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_)::

    pip install cobra

In case you downloaded the source code, run::

    pip install -e .

In the ``cobrapy`` directory. For further information, please follow
the `detailed installation instructions <INSTALL.rst>`_.

Contributing
~~~~~~~~~~~~

Contributions are always welcome! Please read the `contributions
guideline <.github/CONTRIBUTING.rst>`_ to get started.

License
-------

The cobrapy source is released under both the GPL and LGPL licenses. You
may choose which license you choose to use the software under. However,
please note that binary packages which include GLPK (such as the binary
wheels distributed at https://pypi.python.org/pypi/cobra) will be bound
by its license as well.

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License or the Lesser GNU
General Public License as published by the Free Software Foundation,
either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

.. |Build Status| image:: https://travis-ci.org/opencobra/cobrapy.svg?branch=master
   :target: https://travis-ci.org/opencobra/cobrapy
.. |Coverage Status| image:: https://codecov.io/github/opencobra/cobrapy/coverage.svg?branch=master
   :target: https://codecov.io/github/opencobra/cobrapy
.. |Build status| image:: https://ci.appveyor.com/api/projects/status/qmqdy67jbwbc3ds7/branch/master?svg=true
   :target: https://ci.appveyor.com/project/cobrapy39491/cobrapy/branch/master
.. |PyPI| image:: https://badge.fury.io/py/cobra.svg
   :target: https://pypi.python.org/pypi/cobra
.. |Gitter| image:: https://badges.gitter.im/opencobra/cobrapy.svg
   :target: https://gitter.im/opencobra/cobrapy?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
