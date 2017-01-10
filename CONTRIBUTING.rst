============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs using the `issue tracker <https://github.com/opencobra/cobrapy/issues>`__  

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub `issues <https://github.com/opencobra/cobrapy/issues>`__ for bugs. Anything tagged with "bug" and "help wanted" is open to whoever wants to
implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub `issues <https://github.com/opencobra/cobrapy/issues>`__ and `projects <https://github.com/opencobra/cobrapy/projects>`__ for features. Anything tagged with "enhancement" and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

cobrapy could always use more documentation, whether as part of the official cobrapy docs, in docstrings, or even on the web in blog posts, articles, and such - all contributions are welcome!

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an `issue <https://github.com/opencobra/cobrapy/issues>`__.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

If you like cobrapy please remember to 'star' our github page (click on the star at top right corner), that way we also have an idea of who is using cobrapy!

Get Started!
------------

Want to contribute a new feature or improvement? Consider starting by raising an issue and assign it to yourself to
describe what you want to achieve. This way, we reduce the risk of duplicated efforts and you may also get
suggestions on how to best proceed, e.g. there may be half-finished work in some branch that you could start with.

Here's how to set up `cobrapy` for local development to contribute smaller features or changes that you can implement yourself.

1. Fork the `cobrapy` repository on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/cobrapy.git

3. Install libglpk using your package manager. For macOS::

	$ brew install homebrew/science/glpk

   For Debian-based Linux systems (including Ubuntu and Mint)::

	$ sudo apt-get install libglpk-dev

4. If virtualenvwrapper is not installed, `follow the directions <https://virtualenvwrapper.readthedocs.io/en/latest/>`__
   to install virtualenvwrapper.

5. Install your local copy of cobrapy into a virtualenv with virtualenvwrapper::

    $ cd cobrapy
    $ mkvirtualenv cobrapy

   Use the ``--python`` option to select a specific version of Python for the virtualenv. Note on macOS, matplotlib
   requires Python be installed as a framework but virtualenv creates a non-framework build of Python.
   See the `matplotlib FAQ <http://matplotlib.org/1.5.3/faq/virtualenv_faq.html>`__ for details
   on a workaround.

6. Install the required packages for development in the virtualenv using pip install::

   (cobrapy)$ pip install --upgrade pip
   (cobrapy)$ pip install -r develop-requirements.txt

7. Setup cobrapy for development::

    (cobrapy)$ python setup.py develop

8. Create a branch for local development (see below for details on the branching model)::

    (cobrapy)$ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

9. When you are done making changes, check that your changes pass pep8
   and the tests with tox for the supported Python versions::

    (cobrapy)$ tox -e py27
    (cobrapy)$ tox -e py34
    (cobrapy)$ tox -e py35

10. Commit your changes and push your branch to GitHub::

    (cobrapy)$ git add .
    (cobrapy)$ git commit -m "Your detailed description of your changes."
    (cobrapy)$ git push origin name-of-your-bugfix-or-feature

11. Submit a pull request through the GitHub website.

For larger features that you want to work on collaboratively with other cobrapy team members, you may consider to first request to join the cobrapy developers team to get write access to the repository so that you can create a branch in the main repository (or simply ask the maintainer to create a branch for you). Once you have a new branch you can push your changes directly to the main repository and when finished, submit a pull request from that branch to ``devel``.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests in the ``cobra/test``
   directory. Except in rare circumstances, code coverage must
   not decrease (as reported by codecov which runs automatically when
   you submit your pull request)
2. If the pull request adds functionality, the docs should be
   updated. Put your new functionality into a function with a
   docstring and consider creating a notebook that demonstrates the
   usage in ``documentation_builder`` (documentation is written as
   jupyter notebooks in the ``documentation_builder`` directory, which
   are then converted to rst by the ``autodoc.sh`` script.)
3. The pull request should work for Python 2.7, 3.4 and 3.5. Check
   https://travis-ci.org/biosustain/cobrapy/pull_requests
   and make sure that the tests pass for all supported Python versions.
4. Assign a reviewer to your pull request. If in doubt, assign Henning
   Redestig. Your pull request must be approved by at least one
   reviewer before it can be merged.

Unit tests and benchmarks
-------------------------

cobrapy uses `pytest <http://docs.pytest.org/en/latest/>`_ for its
unit-tests and new features should in general always come with new
tests that make sure that the code runs as intended. Since COBRA
rapidly can become quite resource intensive fundamental methods such
as model manipulation, adding and removing reactions, metabolites etc
also must work efficiently. We use `pytest-benchmark
<https://pytest-benchmark.readthedocs.io/en/latest/>`_ to compare
different implementations to make sure that new code do not come with
unacceptable increased computation time. If you add benchmarked tests,
make sure to also include a test with and without the benchmark as we
do not want to slow down continuous integration by running benchmarks,
for examples, see e.g. ``test_add_metabolite`` in `test_model.py
<cobra/test/test_model.py>`_. ``test_add_metabolite`` is the main
test, ``test_add_metabolite_benchmark`` takes the special
``benchmark`` fixture that enables profiling the important code
snippet but is skipped when running::

    (cobrapy)$ pytest --benchmark-skip

When the test function itself is small and can safely be assumed to
not take many resources, we can directly profile the test as in
``test_subtract_metabolite_benchmark`` which calls
``benchmark(self.test_subtract_metabolite, model)``.

To run all tests and benchmarks do::

    (cobrapy)$ pytest

and to compare two implementations you may keep them in two branches
e.g. ``old`` and ``new`` and then do::

    (cobrapy)$ git checkout old
    (cobrapy)$ pytest --benchmark-save
    (cobrapy)$ git checkout new
    (cobrapy)$ pytest --benchmark-compare


Branching model
---------------

``devel``
    Is the branch all pull-requests should be based on.
``master``
    Is only touched by maintainers and is the branch with only tested, reviewed code that is released or ready for the
    next release.
``{fix, bugfix, doc, feature}/descriptive-name``
    Is the recommended naming scheme for smaller improvements, bugfixes, documentation improvement and new features respectively.

Please use concise descriptive commit messages and consider using ``git pull --rebase`` when you update your own fork to avoid merge commits.

1. Tests are in the ``cobra/test`` directory. They are automatically run
   through continuous integration services on both python 2 and python 3
   when pull requests are made.
2. Please write tests for new functions. Writing documentation as well
   would also be very helpful.
3. Ensure code will work with both python 2 and python 3. For example,
   instead of ``my_dict.iteritems()`` use ``six.iteritems(my_dict)``

Thank you very much for contributing to cobrapy!
