================================================================
COBRApy - Constraint-Based Reconstruction and Analysis in Python
================================================================

.. image:: https://img.shields.io/pypi/v/cobra.svg
   :target: https://pypi.org/project/cobra/
   :alt: Current PyPI Version

.. image:: https://img.shields.io/pypi/pyversions/cobra.svg
   :target: https://pypi.org/project/cobra/
   :alt: Supported Python Versions

.. image:: https://img.shields.io/pypi/l/cobra.svg
   :target: https://www.gnu.org/licenses/old-licenses/lgpl-2.0.html
   :alt: GNU Lesser General Public License 2 or later

.. image:: https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg
   :target: https://github.com/opencobra/cobrapy/blob/devel/.github/CODE_OF_CONDUCT.md
   :alt: Code of Conduct

.. image:: https://github.com/opencobra/cobrapy/workflows/CI-CD/badge.svg
   :target: https://github.com/opencobra/cobrapy/workflows/CI-CD
   :alt: GitHub Actions CI/CD Status

.. image:: https://codecov.io/gh/opencobra/cobrapy/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/opencobra/cobrapy
   :alt: Codecov

.. image:: https://readthedocs.org/projects/cobrapy/badge/?version=latest
   :target: https://cobrapy.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://badges.gitter.im/opencobra/cobrapy.svg
   :target: https://gitter.im/opencobra/cobrapy
   :alt: Gitter Chat Room

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/ambv/black
   :alt: Black

.. image:: https://zenodo.org/badge/6510063.svg
   :target: https://zenodo.org/badge/latestdoi/6510063
   :alt: Zenodo DOI

What is COBRApy?
================

COBRA methods are widely used for genome-scale modeling of metabolic networks in
both prokaryotes and eukaryotes. **COBRApy** is a constraint-based modeling
package that is designed to accommodate the biological complexity of the next
generation of COBRA models and provides access to commonly used COBRA methods,
such as flux balance analysis, flux variability analysis, and gene deletion
analyses.

Our aim with COBRApy is to provide useful, efficient infrastructure for:

- creating and managing metabolic models
- accessing popular solvers
- analyzing models with methods such as FVA, FBA, pFBA, MOMA etc.
- inspecting models and drawing conclusions on gene essentiality, testing
  consequences of knock-outs etc.

Our goal with COBRApy is for it to be useful on its own, and for it to be the
natural choice of infrastructure for developers that want to build new COBRA
related python packages for e.g. visualization, strain-design and data driven
analysis. By re-using the same classes and design principles, we can make new
methods both easier to implement and easier to use, thereby bringing the power
of COBRA to more researchers.

The documentation is browseable online at `readthedocs
<https://cobrapy.readthedocs.org/en/stable/>`_ and can also be `downloaded
<https://readthedocs.org/projects/cobrapy/downloads/>`_.

Please use the `Google Group <http://groups.google.com/group/cobra-pie>`_ for
help. By writing a well formulated question, with sufficient detail, you are
much more likely to quickly receive a good answer! Please refer to these
`StackOverflow guidelines <https://stackoverflow.com/help/how-to-ask>`_ on how
to ask questions.  Alternatively, you can use `gitter.im
<https://gitter.im/opencobra/cobrapy>`_ for quick questions and discussions
about COBRApy (faster response times). Please keep in mind that answers are
provided on a volunteer basis.

More information about opencobra is available at the `website
<http://opencobra.github.io/>`_.

If you use COBRApy in a scientific publication, please cite
`doi:10.1186/1752-0509-7-74 <http://dx.doi.org/doi:10.1186/1752-0509-7-74>`_

Installation
============

Use pip to `install COBRApy from PyPI <https://pypi.org/project/cobra/>`_ (we
recommend doing this inside a `virtual environment
<http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_)::

    pip install cobra

If you want to load MATLAB models, you will need additional dependencies. Please
install::

    pip install cobra[array]

For further information, please follow the `detailed installation instructions
<INSTALL.rst>`_.

Contributing
============

Contributions are always welcome! Please read the `contributing guidelines
<https://github.com/opencobra/cobrapy/blob/devel/.github/CONTRIBUTING.rst>`_ to
get started.

License
=======

The COBRApy source is released under both the GPL and LGPL licenses version 2 or
later. You may choose which license you choose to use the software under.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License or the GNU Lesser General Public
License as published by the Free Software Foundation, either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.
