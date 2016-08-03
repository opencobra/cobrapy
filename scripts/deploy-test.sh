#!/bin/bash
echo -e " ... running twine to deploy ... "
echo -e "
[distutils]
index-servers=
    pypirepository

[pypirepository]
repository: https://testpypi.python.org/pypi
username: henred
password: pippi_langstrump
" > ~/.pypirc

pip install twine
twine upload --skip-existing ${TRAVIS_BUILD_DIR}/wheelhouse/* -r pypirepository
