Contribution Guidelines
-----------------------

Generally, the following practices are recommended for making contributions to
cobrapy. These aren't all necessarily hard-and-fast rules, but should serve as
guidelines in most cases.

1. Please comment code.
2. All new python code should be pep8 compliant.
3. Please use git best practices, with a 50 line summary for each commit.
   Generally, separate features should be made in separate commits so
   they can be tested and merged independently. For example, adding a new
   solver would be a separate commit from fixing whitespace in cobra.core.
4. Documentation is written as IPython/jupyter notebooks in the
   ```documentation_builder``` directory, which are then converted to
   rst by the ```autodoc.sh``` script.
5. Tests are in the ```cobra/test``` directory. They are automatically run
   through continuous integration services on both python 2 and python 3
   when pull requests are made.
6. Please write tests for new functions. Writing documentation as well
   would also be very helpful.
7. Ensure code will work with both python 2 and python 3. For example,
   instead of ```my_dict.iteritems()``` use ```six.iteritems(my_dict)```

Thank you very much for contributing to cobrapy.
