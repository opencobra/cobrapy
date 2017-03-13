rm cobra.rst cobra.*.rst
sphinx-apidoc -o . ../cobra \
    ../cobra/test ../cobra/solvers/ \
    ../cobra/solvers/legacy.py
rm modules.rst
