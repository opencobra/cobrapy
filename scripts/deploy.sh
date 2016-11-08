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

if [[ "$TRAVIS_BRANCH" == "devel" ]]; then
	twine upload --skip-existing ${TRAVIS_BUILD_DIR}/wheelhouse/* -r pypirepository
elif [[ "$TRAVIS_BRANCH" == "master" ]]; then
	twine upload --skip-existing --username "${PYPI_USERNAME}" --password "${PYPI_PASSWORD}" ${TRAVIS_BUILD_DIR}/wheelhouse/*
fi
