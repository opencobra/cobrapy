rm cobra.rst cobra.*.rst
sphinx-apidoc -o . ../cobra ../cobra/external \
    ../cobra/test ../cobra/solvers/ ../cobra/test_all.py \
    ../cobra/version.py ../cobra/solvers/legacy.py
rm modules.rst
