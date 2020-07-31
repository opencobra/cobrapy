============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Report Bugs
===========

Report bugs at https://github.com/opencobra/cobrapy/issues.

If you are reporting a bug, please follow the template guide lines. The more
detailed your report, the easier and thus faster we can help you.

Fix Bugs
========

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
==================

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
===================

cobrapy could always use more documentation, whether as part of the official
documentation, in docstrings, or even on the web in blog posts, articles, and
such.

Submit Feedback
===============

The best way to send feedback is to file an issue at
https://github.com/opencobra/cobrapy/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions are
  welcome :)

Get Started!
============

Ready to contribute? Here's how to set up cobrapy for local development.

1. Fork the https://github.com/opencobra/cobrapy repository on GitHub. If you
   have never done this before, `follow the official guide
   <https://guides.github.com/activities/forking/>`_
2. Clone your fork locally as described in the same guide.
3. Install your local copy into a a Python virtual environment.  You can `read
   this guide to learn more
   <https://realpython.com/python-virtual-environments-a-primer/>`_ about them
   and how to create one. Alternatively, particularly if you are a Windows or
   Mac user, you can also use `Anaconda <https://docs.anaconda.com/anaconda/>`_.
   Assuming you have virtualenvwrapper installed, this is how you set up your
   fork for local development

   .. code-block:: console

       mkvirtualenv my-env
       cd cobrapy/
       pip install -e ".[development]"

4. Create a branch for local development using the ``devel`` branch as a
   starting point. Use ``fix``, ``refactor``, or ``feat`` as a prefix

   .. code-block:: console

       git checkout devel
       git checkout -b fix-name-of-your-bugfix

   Now you can make your changes locally.

5. When making changes locally, it is helpful to ``git commit`` your work
   regularly. On one hand to save your work and on the other hand, the smaller
   the steps, the easier it is to review your work later. Please use `semantic
   commit messages
   <http://karma-runner.github.io/2.0/dev/git-commit-msg.html>`_.

   .. code-block:: console

       git add .
       git commit -m "fix: Your summary of changes"

6. When you're done making changes, check that your changes pass our test suite.
   This is all included with tox

   .. code-block:: console

       tox

   You can run all tests in parallel using detox. To get detox, just pip install
   it into your virtualenv.

7. Push your branch to GitHub.

   .. code-block:: console

       git push origin fix-name-of-your-bugfix

8. Open the link displayed in the message when pushing your new branch in order
   to submit a pull request. Please follow the template presented to you in the
   web interface to complete your pull request.

For larger features that you want to work on collaboratively with other cobrapy
team members, you may consider to first request to join the cobrapy developers
team to get write access to the repository so that you can create a branch in
the main repository (or simply ask the maintainer to create a branch for you).
Once you have a new branch you can push your changes directly to the main
repository and when finished, submit a pull request from that branch to
``devel``.

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
3. The pull request will be tested for several different Python versions.
4. Someone from the @opencobra/cobrapy-core team will review your work and guide
   you to a successful contribution.

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
``stable``
    Is only touched by maintainers and is the branch with only tested, reviewed
    code that is released or ready for the next release.
``{fix, bugfix, doc, feature}/descriptive-name``
    Is the recommended naming scheme for smaller improvements, bugfixes,
    documentation improvement and new features respectively.

Please use concise descriptive commit messages and consider using
``git pull --rebase`` when you update your own fork to avoid merge commits.

Thank you very much for contributing to cobrapy!

FAQs
----

Q1. Why do all of the tests that involve loading a pickled model fail on my branch?
	A: Pickling is the standard method for serializing objects in python,
	which is commonly done during operations like multiprocessing.
	Because of this, we need to maintain tests that run on pickled
	models, otherwise contributors may inadvertantly break
	multiprocessing features. If changes you made to cobrapy
	modify attributes of the ``cobra.Model`` class, the pickled
	models stored in the repository won't contain those changes
	and may fail tests that you add or modify. To resolve these
	errors, just run ``cobra/test/data/update_pickles.py`` on your
	branch, which will repickle the models.
