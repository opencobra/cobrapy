rm cobra.rst cobra.*.rst
sphinx-apidoc -o . ../cobra ../cobra/oven ../cobra/external \
    ../cobra/test ../cobra/solvers/*_java.py ../cobra/test_all.py \
    ../cobra/version.py ../cobra/solvers/legacy.py
rm modules.rst

ipython nbconvert --to=rst *.ipynb

ipython nbconvert --to=python --template=pyscript.tpl *.ipynb
ls -1 *py | grep -v ^conf.py$ | xargs -Ipy mv py ../examples
